#!/usr/bin/python
"""
With this module you can interact with the Solv@TUM database.
The usage of this database interface is described in the README.

Please note:
This module was written for Python2.7.
The compatibility with Python3 has not been tested.

Last Update 18.08.2018
author: Christoph Hille (c.hille@tum.de)
"""

import sys
import os
from os.path import join as pathjoin
import time
from subprocess import call
import math
import pybel

try:
    from ase.io import sdf
except ImportError:
    pass

try:
    import bibtexparser
except ImportError:
    pass

try:
    import imolecule
except ImportError:
    pass

try:
    from scipy import constants
except ImportError:
    pass

try:
    import pandas as pd
except ImportError:
    pass


class Database(object):
    """
    By initialising of this class, the content of the SD-file will be loaded by using the Python module Pybel.
    Afterwards different methods can be applied on the database object.
    For more information have a look into the README.
    options:
        - path:             -   str; path to the sdf databasefile
                                default: same directory
        - name:             -   str; databasename
                                default: 'solvatum.sdf'
        - energy_unit:      -   str; 'eV', 'kcal' or 'J'
                                (unit for the energy)
                                default: eV
        - temperature:      -   float; temperature in K
                                (for energy calculation)
                                default: 298.15
    """
    def __init__(self,
                 path=None,
                 name='solvatum.sdf',
                 energy_unit='eV',
                 temperature=298.15):

        if path:
            self.path = path
        else:
            self.path = os.path.dirname(__file__)

        self.sdf_file = list(pybel.readfile("sdf", pathjoin(self.path, name)))
        self.solutes = self.__get_all_solutes()
        self.solvents = self.__get_all_solvents()
        self.energy_unit = energy_unit
        self.temperature = temperature

        if "imolecule" in sys.modules:
            pybel.ipython_3d = True

    def __get_all_solutes(self):
        """
        Returns a dictionary of all solute-molecules in the database
        with the unique ID as key and the solute name as value
        """

        solutes = {}

        for mol in self.sdf_file:
            solutes[mol.title] = mol.data['Name'].upper()

        return solutes

    def __get_all_solvents(self):
        """
        Returns dictionary of all solvent-molecules in the database
        with the solvent name as key and the number of solutes,
        for which reference data exists

        Example:
            {'METHANOL': 104, ...}
            There are 104 reference values for Methanol as solvent.
        """

        solvents = {}

        for mol in self.sdf_file:
            data = mol.data
            for name in data:
                if 'logK' in name:
                    sol = name[6:-1].upper()
                    if sol not in solvents:
                        solvents[sol] = 1
                    else:
                        solvents[sol] += 1

        return solvents

    def __check_if_molecules_in_database(self, molecule, issolvent=False, disp=True):
        """
        Checks for a given molecule if it has a geometry entry in the database.
        """
        if not isinstance(molecule, str):
            raise TypeError('Molecule has to be given as a string!')

        sol = molecule

        if sol not in self.solutes.keys() and sol.upper() not in self.solutes.values():
            if not issolvent:
                raise KeyError('The given solute %s could not be identified as a solute in the database' % sol)
            else:
                if sol.upper() not in self.solvents:
                    raise KeyError('The given solvent %s could not be identified as a solvent in the database' % sol)
                if disp:
                    print('There is no geometry entry for solvent %s in database' % sol)
                return 'no geometry'
        return 'has geometry'

    def __name_id_handler(self, molecule, direction='id', issolvent=False, disp=True):
        """
        Converts soluteID to solute name and vice versa
        (eg:    INPUT: molecule='Argon', direction='id' -> OUTPUT: '002'
                INPUT: molecule='081', direction='name'	-> OUTPUT: 'METHANOL'
        """
        if not isinstance(molecule, str):
            raise TypeError('Molecule has to be given as a string!')

        if direction == 'id':
            if issolvent:
                raise RuntimeError("Solvents have no database ID.")
            elif molecule not in self.solutes.keys():
                molecule = self.solutes.keys()[self.solutes.values().index(molecule.upper())]

        elif direction == 'name':
            if self.__check_if_molecules_in_database(molecule, issolvent=issolvent, disp=disp) == 'no geometry':
                molecule = molecule.upper()
            elif molecule.upper() in self.solutes.values():
                molecule = molecule.upper()
            else:
                molecule = self.solutes[molecule]

        return molecule

    def __filter_solvent(self, solvent, solutes='all', disp=True):
        """
        Returns dictionary with solvent properties and a list of soluteIDs,
        solutes and logK values for the given solvent.
        """
        solvent = solvent.upper()

        if solvent not in self.solvents:
            raise KeyError('Solvent %s was not found in the database' % solvent)

        nogeo = bool(self.__check_if_molecules_in_database(solvent, issolvent=True, disp=disp) == 'no geometry')

        # to avoid problem with unicode:
        if solvent == '(\\XB1)-2-BUTANOL':
            key = 'logK ((\\xb1)-2-BUTANOL)'
        elif solvent == '(\\XB1)-1,2-PROPANEDIOL':
            key = 'logK ((\\xb1)-1,2-PROPANEDIOL)'
        else:
            key = 'logK (' + solvent + ')'

        if solutes == 'all':
            solutes = self.sdf_file
        solutes_list = []
        for mol in solutes:
            if key in mol.data.keys():
                log_k = float(mol.data[key])
                delta_g_solv = calculate_g_solv(log_k, unit=self.energy_unit, temp=self.temperature)
                solutes_list += [[mol.title, mol.data['Name'], log_k, delta_g_solv]]

        if nogeo:
            props = {}
        else:
            props = self.get_molecule_properties(solvent)

        return props, solutes_list

    def one_mol_from_sdf(self, solute):
        """
        Returns the pybel mol-object for a given solute (name or id).
        """
        if not isinstance(solute, str):
            raise TypeError('Solute has to be given as a single string')

        solute = self.__name_id_handler(solute)

        return self.sdf_file[int(solute)]

    def __filter_solute(self, solute):
        """
        Returns a list of solvents and logK values for a given solute.
        """

        mol = self.one_mol_from_sdf(solute.upper())

        solvents_list = []
        for key in mol.data.keys():
            if key[0:4] == 'logK':
                solvents_list += [[key[6:-1],
                                   mol.data[key],
                                   calculate_g_solv(float(mol.data[key]),
                                                    unit=self.energy_unit,
                                                    temp=self.temperature)]]
        props = self.get_molecule_properties(solute)

        return props, solvents_list

    def __filter_solvent_solute(self, solvent, solute, disp=True):
        """
        Returns the logK and solvation free energy value for a given solvent, solute combinaton
        """

        solutemol = [self.one_mol_from_sdf(solute)]

        solvent_solute = self.__filter_solvent(solvent, solutes=solutemol, disp=disp)[1]
        if not solvent_solute:
            raise KeyError('No logK value for solute %s and solvent %s was found in the database' % (solute, solvent))

        log_k, delta_g_solv = solvent_solute[0][2:]
        return {'logK': log_k, 'deltaG_solv': delta_g_solv}

    def get_molecule_properties(self, molecule):
        """
        Returns dictionary with all data stored in the database for one molecule (except the logL values).
        (e.g.: SMILES, dipole moment, InChI, ...)
        """

        mol = self.one_mol_from_sdf(molecule)

        props = {}
        props['molweight'] = round(float(mol.molwt), 2)
        props['formula'] = mol.formula
        props['databaseID'] = mol.title

        for key in mol.data:
            if key[:4] != "logK":
                props[key] = mol.data[key]

        return props

    def create_sol_props_df(self):
        """
        Returns pandas data frame with solute dipole moment and solute polarizability.
        This requires the Python package pandas.
        """
        if "pandas" not in sys.modules:
            print(r"You have to install Pandas before you can use this method.")
            return None
        
        solprops = {}
        for sol in self.solutes:
            solprops[int(sol)] = self.filtering(solute=sol)[0]
        solprops = pd.DataFrame(solprops).T
        solprops['solute_dipole'] = solprops['dipole moment']
        solprops['solute_polari'] = solprops['mean polarizability']
        solprops = pd.DataFrame(data=solprops, columns=['solute_dipole', 'solute_polari'])

        return solprops

    def filtering(self, solvent=None, solute=None, disp=True):
        """
        Filtering the data for one solvent, solute or both.
        """
        if solute is not None and solvent is not None:
            if not isinstance(solute, str) or not isinstance(solvent, str):
                raise TypeError('Solute/ Solvent has to be given as a single string')
            return self.__filter_solvent_solute(solvent.upper(), solute.upper(), disp=disp)
        if solute is not None:
            if not isinstance(solute, str):
                raise TypeError('Solute has to be given as a single string')
            return self.__filter_solute(solute.upper())
        if solvent is not None:
            if not isinstance(solvent, str):
                raise TypeError('Solvent has to be given as a single string')
            return self.__filter_solvent(solvent.upper(), disp=disp)
        return None

    def add_property(self, molecule, prop_name, prop_value, name='solvatum.sdf'):
        """
        Adding further properties to the database.
        """
        if not isinstance(molecule, str):
            raise TypeError('Molecule has to be given as a string')

        molecule = self.__name_id_handler(molecule)

        output = pybel.Outputfile('sdf', pathjoin(self.path, name), overwrite=True)

        for mol in self.sdf_file:
            if molecule == mol.title:
                mol.data[prop_name] = prop_value
            output.write(mol)

        output.close()

    def seperate(self):
        """
        Seperates the whole SD file into a file per solute
        """
        timestamp = time.strftime("%Y%m%d%H%M%S")
        folder = 'seperated_files_' + timestamp
        os.mkdir(folder)
        os.chdir(folder)

        for mol in self.sdf_file:
            singlefile = pybel.Outputfile("sdf", str(mol.title) + '.sdf', overwrite=True)
            singlefile.write(mol)
            singlefile.close()

    def get_reference(self, solvent):
        """
        Returns the reference of the partition coefficient for one solvent.
        A proper usage of this method requires the installation of bibtexparser
        (http://bibtexparser.readthedocs.io/en/master/index.html)
        """

        # for solvents without a solute entry in the database the bibtex references have to be set here:
        further_refs = {'1-ETHYL-2-PYRROLIDINONE': 'Krummen2002',
                        '1-HEXADECENE': 'Abraham2012',
                        '1-METHYL-2-PIPERIDINONE': 'Abraham2009',
                        '1,5-DIMETHYL-2-PYRROLIDINONE': 'Krummen2002',
                        '5,8,11,14-TETRAOXAOCTADECANE': 'Hart2017',
                        'DIBENZYL ETHER': 'Park1987',
                        '1,2-DICHLOROETHANE': 'Sprunger2008a',
                        'DIETHYLDIGLYCOL': 'Park1987',
                        'DIETHYLENE GLYCOL DIBUTYL ETHER': 'Hart2017',
                        'ETHYLBENZOATE': 'Topphoff2000',
                        'ETHYLENE GLYCOL': 'Abraham2010',
                        'HEPTYLACETATE': 'Based on activity coefficient equal to unity',
                        'N-ETHYLACETAMIDE': 'Abraham2009',
                        'N-ETHYLFORMAMIDE': 'Abraham2009',
                        'N,N-DIBUTYLFORMAMIDE': 'Abraham2009',
                        'N,N-DIETHYLACETAMIDE': 'Abraham2009',
                        'PENTAN-3-OL': 'Miyano2005, Miyano2006',
                        'SULFOLANE': 'Stephens2011',
                        'TETRADECANE': 'Vrbka2002, Mokbel1998',
                        'TETRAETHYLENE GLYCOL DIMETHYL ETHER': 'Hart2017'}

        if "bibtexparser" in sys.modules:
            os.chdir(self.path)
            os.chdir('..')
            mainpath = os.getcwd()

            try:
                with open(pathjoin(mainpath, 'references', 'solvatum_references.bib')) as bibtex_file:
                    bib_database = bibtexparser.load(bibtex_file)
                    refs = bib_database.entries_dict
                    bibtex = True
            except IOError:
                print("WARNING: The bibtex file 'solvatum_references.bib' could not be found.\n\n")
                bibtex = False

        else:
            bibtex = False
            print(r"WARNING: A proper usage of this feature requires the installation of bibtexparser." + "\n"
                  r"Go to http://bibtexparser.readthedocs.io/en/master/index.html for more informations." + "\n"
                  r"If this is not installed, this method just prints the bibtex key " +
                  r"and the corresponding reference has to be manually looked up " +
                  r"in 'references/solvatum_references.bib'" + "\n\n")

        solvent = self.__name_id_handler(solvent, direction='name', issolvent=True, disp=False)

        if solvent in self.solutes.values():
            try:
                ref_per_solvent = self.get_molecule_properties(solvent)['reference']
            except KeyError:
                print("No reference for solvent %s deposited" % solvent)
                return
        else:
            try:
                ref_per_solvent = further_refs[solvent]
            except KeyError:
                print("No reference for solvent %s deposited" % solvent)
                return
        ref_per_solvent = ref_per_solvent.split(", ")

        print('------------------------------\n')
        for ref in ref_per_solvent:
            if ref == "Based on activity coefficient equal to unity":
                print(ref)
                continue
            if bibtex:
                ref = refs[ref]
                keys = [u'title', u'author', u'journal', u'year', u'volume', u'pages', u'doi', u'url', 'ID']
                for key in keys:
                    try:
                        entry = ref[key]
                        print(key + ': ' + entry)
                    except KeyError:
                        pass
            else:
                print(ref)
            print('------------------------------\n')

        return

    def draw(self, molecule, title=None, notebook=False, save=None):
        """
        Drawing the depiction of a given molecule (one can also give a list of molecules as input).

        If using this method in jupyter notebooks, one can set 'notebook=True' for inline printing.

        Per default the database ID and the molecule name are printed as title.
        For printing another title use the key word 'title'.
        For printing no title, set 'title='' '.

        For saving the figure use 'save=<folder name>'. The folder will be created in the current directory.
        The picture will be saved as a png.
        """
        if isinstance(molecule, str):
            molecule = [molecule]

        for mol in molecule:
            mol = self.one_mol_from_sdf(mol)

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
                    os.rename('temp.png', pathjoin(save, mol.title.replace(' ', '_') + '.png'))
                    print('Picture saved as ', pathjoin(save, mol.title.replace(' ', '_') + '.png'))
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
                    mol.draw(show=False, filename=pathjoin(save, mol.title.replace(' ', '_') + '.png'))
                    print('Picture saved as ', pathjoin(save, mol.title.replace(' ', '_') + '.png'))
                else:
                    mol.draw()

                mol.title = oldtitle
                return None

    def mol_to_ase(self, molecule):
        """
        Return the geometry of a molecule as ASE Atoms object.
        ASE has to be installed
        (it is even enough to clone the git repo https://gitlab.com/ase/ase and set the Python path).
        More informations: https://wiki.fysik.dtu.dk/ase/
        """

        if "ase" not in sys.modules:
            print(r"You have to install ASE before you can use this feature." + "\n"
                  r"Go to https://wiki.fysik.dtu.dk/ase/ for more informations.")
            return None

        mol = self.one_mol_from_sdf(molecule)

        output = pybel.Outputfile('sdf', ".tmp.sdf", overwrite=True)
        output.write(mol)

        ase_atoms = sdf.read_sdf(".tmp.sdf")
        os.remove(".tmp.sdf")

        return ase_atoms

    def d3_viewer(self, molecule, viewer='avogadro'):
        """
        Opens the geometry of the molecule in a 3d viewer.
        Currently only avogadro is supported, but you can test it with other programms as well.
        """

        mol = self.one_mol_from_sdf(molecule)

        output = pybel.Outputfile('sdf', ".tmp.sdf", overwrite=True)
        output.write(mol)

        call([viewer, '.tmp.sdf'])
        os.remove(".tmp.sdf")

    def sol_has_ring(self, solute):
        """
        Checks if solute has ring system.

        Returns bool.
        """
        mol = self.one_mol_from_sdf(solute)

        for atom in mol.atoms:
            if atom.OBAtom.CountRingBonds() > 0:
                return True

        return False

    def sol_is_aromatic(self, solute):
        """
        Checks if solute requires the necessary conditions for an aromatic compound.
        ATTENTION: it is not sufficient. The solute can also be an anti-aromatic compound.

        Returns bool.
        """
        if not self.sol_has_ring(solute):
            return False

        mol = self.one_mol_from_sdf(solute)

        for atom in mol.atoms:
            if atom.OBAtom.CountRingBonds() > 0:
                if atom.hyb != 2:
                    return False

        return True

    def atoms_in_sol(self, solute):
        """
        Returns the number of atoms in a given solute
        """
        mol = self.one_mol_from_sdf(solute)

        return len(mol.atoms)

    def sol_has_conj(self, solute):
        """
        Returns if solute has a conjugated system (an atom which is sp1 or sp2 hybridized)
        """
        mol = self.one_mol_from_sdf(solute)

        for atom in mol.atoms:
            if atom.hyb != 0 and atom.hyb != 3:
                return True

        return False


def calculate_g_solv(log_k, unit='eV', temp=298.15):
    """
    Calculates the solvation free energy for a given logK value in kcal/mol, eV or J.
    Formula: delta_g_solv = -ln(10)*RT*log_k

    INPUT:
        - log_k:     float; Partition Coefficient
        - unit:     string; energy unit (kcal/mol, eV or J)
        - temp:     float; temperature
    """

    if 'scipy' in sys.modules:
        con_r = constants.R
        con_cal = constants.calorie
        con_kilo = constants.kilo
        con_na = constants.N_A
        con_ev = constants.eV
    else:
        con_r = 8.3144598
        con_cal = 4.184
        con_kilo = 1000.0
        con_na = 6.022140857e+23
        con_ev = 1.6021766208e-19

    delta_g_solv = - con_r * temp * math.log(10) * log_k  # delta_g_solv in joule

    if unit == 'kcal/mol':
        # conversion joule to kcal/mol
        delta_g_solv = round(delta_g_solv / (con_cal * con_kilo), 2)
    elif unit == 'eV':
        # conversion joule to eV
        delta_g_solv = round(delta_g_solv / (con_na * con_ev), 5)
    elif unit == 'J':
        delta_g_solv = round(delta_g_solv, 2)
    else:
        raise IOError("The unit %s is not supported as energy unit." % unit +
                      "Please choose between 'kcal/mol', 'eV' or 'J'")

    return delta_g_solv
