#!/usr/bin/python
"""
With this module you can interact with the solvation free energy database.
The usage of the database is described in the README.

Please note:
This module was written for Python2.7. 
The compatibility with Python3 has not been tested yet.

Last Update 30.06.2018
@author: Christoph Hille
"""

import sys
import os
from os.path import join as pathjoin
import pybel
import numpy as np
import pandas as pd
from scipy import constants
import time
from subprocess import call
import pprint

try:
    from ase.io import sdf
except ImportError:
    pass

try:
    import bibtexparser
except ImportError:
    pass

class Database:
    """
    options:
        - path:         - str; path to the databasefile
        - name:         - str; databasename
                          default: 'database.sdf'
        - energy_unit:  - str; 'eV' or 'kcal'
                          in which unit the energy will be given
                          default: eV
        - temperature:  - float; temperature (for energy calculation);
                          default: 298.15
    """
    
    def __init__(self, path=None,
                       name='database.sdf',
                       energy_unit='eV',
                       temperature=298.15
                       ):
        
        if path:
            self.path = path
        else:
            self.path = os.path.dirname(__file__)

        self.__sdf_file = list(pybel.readfile("sdf", pathjoin(self.path, name)))	    
        self.solutes = self.__get_all_solutes()
        self.solvents = self.__get_all_solvents()
        self.energy_unit = energy_unit
        self.temperature = temperature
            
    def __get_all_solutes(self):
        """
        Returns dictionary of all solute-molecules in the database with solute ID as key and solute name as values
        """
        
        solutes = {}

        for mol in self.__sdf_file:
            solutes[mol.title] = mol.data['Name'].upper()
	
        return solutes

    def __get_all_solvents(self):
        """
        Returns dictionary of all solvent-molecules in the database with the solvent name as key and a tuple of the number of 
        solutes for which reference data exists and the dieelctric constant as values 
        
        Example:
            {'METHANOL': (104, 32.613), ...}
            There exists 104 reference values for Methanol as solvent and the dieelctric constant of Methanol is 32.613
        """
	
        solvents = {}
        eps = {}

        for mol in self.__sdf_file:
            data = mol.data
            for name in data:
                if 'logK' in name:
                    sol = name[6:-1].upper()
                elif 'WATER' in name:
                    sol = 'WATER'
                elif name == 'dielectric constant':
                    eps[data['Name']] = data[name]
                    sol = None
                else:
                    sol = None
                    continue
                if sol not in solvents:
                    solvents[sol] = 1
                else:
                    solvents[sol] += 1
        
        del solvents[None]
        solvents_out = {}
        for sol in solvents:
            if sol in eps:
                solvents_out[sol] = (solvents[sol], eps[sol])
            else:
                solvents_out[sol] = (solvents[sol], None)

        return solvents_out
	
    def __check_if_molecules_in_database(self, molecules_list, issolvent=False, disp=True):
        """
        Checks for every molecule in a given list if it has a geometry entry in the database
        
        I raises an error if a solute doesn't have an geometry entry (because thats means tha)
        """	
        if type(molecules_list) != list:
            raise TypeError('Molecules list has to be a list!')
        
        for sol in molecules_list:
            if sol not in self.solutes.keys() and sol.upper() not in self.solutes.values():
                if not issolvent:
                    raise KeyError('The given solute %s could not be identified as a solute in the database' %sol)
                else:
                    if sol.upper() not in self.solvents:
                        raise KeyError('The given solvent %s could not be identified as a solvent in the database' %sol)
                    if disp:
                        print 'There is no geometry entry for solvent %s in database' %sol
                    return 'no geometry'
	
    def __name_id_handler(self, molecules_list, direction='id', issolvent=False, disp=True):
        """ 
        Converts soluteID to solute name and vice versa
        (eg:    INPUT: solutes_list=[Argon], direction='id' -> OUTPUT: solutes_list=[002]
                INPUT: solutes_list=[081], direction='name'	-> OUTPUT: solutes_list=[METHANOL]
        """
        if self.__check_if_molecules_in_database(molecules_list, issolvent=issolvent, disp=disp) == 'no geometry':
            return ['no geometry']
        
        if direction == 'id':
            for i in range(len(molecules_list)):
                if molecules_list[i] in self.solutes.keys():
                    pass
                else:
                    molecules_list[i] = self.solutes.keys()[self.solutes.values().index(molecules_list[i].upper())]
                
        elif direction == 'name':
            for i in range(len(molecules_list)):
                if molecules_list[i].upper() in self.solutes.values():
                    molecules_list[i] = molecules_list[i].upper()
                else:
                    molecules_list[i] = self.solutes[molecules_list[i]]

        return molecules_list

    def __filter_solvent(self, solvent, solutes='all', disp=True):
        """
        Returns dictionary with solvent properties and a list of soluteIDs, solutes and logK values for the given solvent
        """
        if solutes == 'all':
            solutes = self.__sdf_file
            
        solvent = solvent.upper()
        solventname = self.__name_id_handler([solvent], direction='name', issolvent=True, disp=disp)[0]
        if solventname == 'no geometry':
            nogeo = True
        else:
            nogeo = False
            solvent = solventname
        
        if solvent not in self.solvents:
            print 'Solvent %s was not found in the database' %solvent
            return 

        if solvent == 'WATER':
            key = 'dgexp marzari (WATER)'
        elif solvent == '(\\XB1)-2-BUTANOL':
            key = 'logK ((\\xb1)-2-BUTANOL)'
        elif solvent == '(\\XB1)-1,2-PROPANEDIOL':
            key = 'logK ((\\xb1)-1,2-PROPANEDIOL)'
        else:
            key = 'logK (' + solvent + ')'
        solutes_list = []
        for mol in solutes:
            if key in mol.data.keys():
                if solvent == 'WATER':
                    logK = None
                    deltaGsolv = float(mol.data[key])
                    if self.energy_unit == 'eV':
                        deltaGsolv = round((deltaGsolv * constants.kilo * constants.calorie) / (constants.N_A * constants.eV), 4) #conversion kcal/mol to eV
                else:  
                    logK = float(mol.data[key])
                    deltaGsolv =  self.__calculate_G_solv(logK)
                solutes_list += [[mol.title, mol.data['Name'], logK, deltaGsolv]]

        if nogeo:
            props = {}
        else:
            props = self.__get_molecule_properties(solvent)
            
        return props, solutes_list

    def __one_mol_from_sdf(self, solute):
        """
        Returns the pybel mol-object for a given solute (name or id)
        """

        if type(solute) != str:
            raise TypeError('Solute has to be a single string')

        solute = self.__name_id_handler([solute])[0]

        return self.__sdf_file[int(solute)]

    def __filter_solute(self, solute):
        """
        Returns a list of solvents and logK values for a given solute 
        """
	
        mol = self.__one_mol_from_sdf(solute.upper())

        solvents_list = []
        for key in mol.data.keys():
            if key[0:4] == 'logK':
                solvents_list += [[key[6:-1 ], mol.data[key], self.__calculate_G_solv(eval(mol.data[key]))]]
        if solute.upper() == 'WATER':
            solvents_list +=[['WATER', None, mol.data['dgexp marzari (WATER)']]]
        props = self.__get_molecule_properties(solute)

        return props, solvents_list

    def __filter_solvent_solute(self, solvent, solute, disp=True):
        """
        Returns the logK and solvation free energy value for a given solvent, solute combinaton
        """

        solutemol = [self.__one_mol_from_sdf(solute)]
        
        solvent_solute = self.__filter_solvent(solvent, solutes=solutemol, disp=disp)[1]
        if len(solvent_solute) == 0:
            raise KeyError('No logK value for solute %s and solvent %s was found in the database' %(solute, solvent)) 
            
        logK, deltaGsolv = solvent_solute[0][2:]
        return {'logK': logK, 'deltaG_solv': deltaGsolv}

    def __calculate_G_solv(self, logK):
        """
        Calculates the solvation free energy for a given logK value in kcal/mol or eV.
        Formular: deltaGsolv = -ln(10)*RT*logK

        INPUT:
            - logK:     float; Partition Coefficient

            optional:
                - unit: str; 'kcal' or 'eV'; unit; default: 'eV'
                        (by initianlization of the class object one can also specify the energy unit with 'energy_unit') 
                - T:    float; temperature; default: 298.15
        """  

        unit = self.energy_unit
        T = self.temperature

        delta_G_solv = - constants.R * T * np.log(10) * logK #deltaGsolv in joule
        
        if unit == 'kcal/mol':
            delta_G_solv = round(delta_G_solv / (constants.calorie * constants.kilo), 2) #conversion joule to kcal/mol
        elif unit == 'eV':
            delta_G_solv = round(delta_G_solv / (constants.N_A * constants.eV), 5) # conversion joule to eV
        elif unit == 'J':
            delta_G_solv = round(delta_G_solv, 2)
        else:
            raise IOError("The unit %s is not supported as energy unit in this database. Please choose between 'kcal/mol', 'eV' or 'J'" %unit)

        return delta_G_solv

    def __get_molecule_properties(self, solvent):

        mol = self.__one_mol_from_sdf(solvent)

        eps = 'dielectric constant'
        eps_source = 'dielectric constant source'
        gamma = 'surface tension'
        polar = 'mean polarizability'
        dipole = 'dipole moment'
        polar_exp = 'experimental polarizability'

        props = {}

        props['molweight'] = round(float(mol.molwt),2)
        props['formula'] = mol.formula
        props['InChI'] = mol.data['InChI']
        props['vacuum energy'] = mol.data['vacuum energy']
        props['databaseID'] = mol.title

        if eps in mol.data:
            props[eps] =  float(mol.data[eps])
        if eps_source in mol.data:
            props[eps_source] =  mol.data[eps_source]
        elif eps in mol.data:
            print('There is no source for the dielectric constant value of %s' %solvent)
        if gamma in mol.data:
            props[gamma] = float(mol.data[gamma])
        if polar in mol.data:
            props[polar] = float(mol.data[polar])
        if dipole in mol.data:
            props[dipole] = float(mol.data[dipole])
        if polar_exp in mol.data:
            props[polar_exp] = float(mol.data[polar_exp])
        return props
    
    def create_sol_props_df(self):
        import pandas as pd

        solprops = {}
        for sol in self.solutes:
            solprops[int(sol)] = self.filtering(solute=sol)[0]
        solprops = pd.DataFrame(solprops).T
        solprops['solute_dipole'] = solprops['dipole moment']
        solprops['solute_polari'] = solprops['mean polarizability']
        solprops['solute_polari_exp'] = solprops['experimental polarizability']
        solprops = pd.DataFrame(data=solprops, columns=['solute_dipole', 'solute_polari', 'solute_polari_exp'])
        
        return solprops

    def filtering(self, solvent=None, solute=None, disp=True):

        if solute != None and solvent != None:
            return self.__filter_solvent_solute(solvent.upper(), solute.upper(), disp=disp)
        elif solute != None:
            return self.__filter_solute(solute.upper())
        elif solvent != None:
            return self.__filter_solvent(solvent.upper(), disp=disp)

    def add_dielectric_constant(self, molecule, value, source=None):
	
        self.add_property(molecule, 'dielectric constant', value)
        print('Added dielectric constant epsilon = %f for solute %s' %(value, molecule))

        if source is not None:
            self.add_property(molecule, 'dielectric constant source', source)
            print('Added %s as source for dielectric constant for solute %s' %(source, molecule))

    def add_property(self, molecule, prop_name, prop_value, name='database.sdf'):
	
        if type(molecule) != str:
            raise TypeError('Molecule has to be given as a string')

        molecule = self.__name_id_handler([molecule])[0]

        output = pybel.Outputfile('sdf', pathjoin(self.path, name), overwrite=True)

        for mol in self.__sdf_file:
            if molecule == mol.title:
                mol.data[prop_name] = prop_value
            output.write(mol)

        output.close()

    def seperate(self):
        """
        Seperates the whole SD file into a file per solute
        """
        time = time.strftime("%Y%m%d%H%M%S")
        folder = 'seperated_files_' + time
        os.mkdir(folder)
        os.chdir(folder)

        for mol in self.__sdf_file:
            singlefile = pybel.Outputfile("sdf", str(mol.title) + '.sdf', overwrite=True)
            singlefile.write(mol)
            singlefile.close()
    
    def get_reference(self, solvent):
        """
        Returns the reference of the partition coefficient for one solvent.
        Requires the installation of bibtexparser (http://bibtexparser.readthedocs.io/en/master/index.html)
        """

        if "bibtexparser" not in sys.modules:
            print(r"This feature requires the installation of bibtexparser." + "\n"  
                  r"Go to http://bibtexparser.readthedocs.io/en/master/index.html for more informations.")
        
        solvent = self.__name_id_handler([solvent], direction='name', issolvent=True, disp=False)[0]
        
        os.chdir(self.path)
        os.chdir('..')
        mainpath = os.getcwd()
        
        try:
            with open(pathjoin(mainpath, 'references', 'references.bib')) as bibtex_file:
                bib_database = bibtexparser.load(bibtex_file)
        except IOError:
            raise IOError("The bibtex file 'references.bib' could not be found." +
                          "Ensure the same folder structure as in the git repository")
        
        refs = bib_database.entries_dict
        
        try:
            ref_per_solvent = pd.read_csv(pathjoin(mainpath, 'references', 'references_per_solvent.csv'), index_col=0)
        except IOError:
            raise IOError("The file 'references_per_solvent.csv' could not be found in the folder references." +
                          "Ensure the same folder structure as in the git repository")
        ref_per_solvent = ref_per_solvent['reference']
        ref_per_solvent = ref_per_solvent.dropna()
        
        if solvent in ref_per_solvent:
            print ''
            keys = [u'title',  u'author',  u'journal', u'year', u'volume', u'pages',  u'doi',  u'url', 'ID']
            for key in  keys:
                try:
                    entry = refs[ref_per_solvent[solvent]][key]
                    print key + ': ' + entry
                except KeyError:
                    pass
        else:
            print('No reference for solvent %s deposited' %solvent)
        
    
    def draw(self, molecule, title=None, notebook=False, save=None):
        """
        Drawing a given molecules. 
        
        If using this method in jupyter notebooks, one can set 'notebook=True' for inline printing.
        
        Per default the database ID and the molecule name are printed as title.
        For printing another title use the key word 'title'.
        For printing no title, set 'title='' '.
        
        For saving the figure use 'save=<folder name>'. The folder will be created in the current directory.
        The picture will be saved as a png.
        """
        molecule = self.__name_id_handler([molecule])

        displays = []
        for mol in molecule:
            mol = self.__sdf_file[int(mol)]
            
            oldtitle = mol.title
            if title is not None:
                mol.title = title
            else:
                mol.title = mol.title + ' ' + self.solutes[mol.title]
            if notebook:
                from IPython.display import Image
                mol.draw(show=False, filename='temp.png')
                display = Image('temp.png')
                
                if save:
                    try:
                        os.mkdir(save)
                    except OSError:
                        pass
                    os.rename('temp.png', pathjoin(save, mol.title.replace(' ','_') + '.png'))
                    print 'Picture saved as ', pathjoin(save, mol.title.replace(' ','_') + '.png')
                else:
                    os.remove('temp.png')
                mol.title = oldtitle
                return display
        
            else:
                if save:
                    try:
                        os.mkdir(save)
                    except OSError:
                        pass
                    mol.draw(show=False, filename=pathjoin(save, mol.title.replace(' ','_') + '.png'))
                    print 'Picture saved as ', pathjoin(save, mol.title.replace(' ','_') + '.png')
                else:
                    mol.draw()

            mol.title = oldtitle
            
    def mol_to_ase(self, molecule):
        """
        Return the geometry of a molecule as ASE Atoms object. 
        ASE has to be installed (it is even enough to clone the git repo https://gitlab.com/ase/ase and set the Python path).
        More informations: https://wiki.fysik.dtu.dk/ase/
        """
        
        if "ase" not in sys.modules:
            print(r"You have to install ASE before you can use this feature." + "\n"  
                  r"Go to https://wiki.fysik.dtu.dk/ase/ for more informations.")
            
        molecule = self.__name_id_handler([molecule], disp=False)[0]    
        mol = self.__sdf_file[int(molecule)]
        
        output = pybel.Outputfile('sdf', ".tmp.sdf", overwrite=True)
        output.write(mol)
        
        ase_atoms = sdf.read_sdf(".tmp.sdf")
        os.remove(".tmp.sdf")
        
        return ase_atoms
        
    def d3_viewer(self, molecule, viewer='avogadro'):
        """
        Opens the geometry of the molecule in a 3d viewer. 
        Currently only avogadro is supported, but you can test other as well.
        """
        
        molecule = self.__name_id_handler([molecule], disp=False)[0]    
        mol = self.__sdf_file[int(molecule)]
        
        output = pybel.Outputfile('sdf', ".tmp.sdf", overwrite=True)
        output.write(mol)
        
        call([viewer, '.tmp.sdf'])
        os.remove(".tmp.sdf")
        
    def sol_has_ring(self, solute):
        """
        Checks if solute has ring system.
        
        Returns bool.
        """
        solute = self.__name_id_handler([solute], disp=False)[0]   
        
        mol = self.__one_mol_from_sdf(solute)
        
        for atom in mol.atoms:                    
            if atom.OBAtom.CountRingBonds() > 0:
                return True
        
        return False
        
        
        
    
        
