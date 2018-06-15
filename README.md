This solvation free energy database contains partition coefficients (logK)
of 658 neutrale organic molekules and inorganic gases dissolved in organic solvents.
The values were provided by Prof. Dr. Acree Jr.
All publications of these values (36 in total) are listed in the [bib file](references/references.bib).
The [database](database/database.sdf) also contains in vacuum optimized geometries for each solute.

For a simple and fast access to the data one can use the database interface 
([databasehandler](database/databasehandler.py)), which was written in the Python programming language.

The installation and usage of this database interface is described in the following sections.

# Installation:

## Getting the repository:
    
    $ git clone https://gitlab.lrz.de/hille721/solvation-free-energy-database
    
## Installation of Pybel/OpenBabel

For using the database interface one need to install Pybel and also OpenBabel.
Installation guidelines can be found here: 

* [Pybel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html)
* [OpenBabel](http://openbabel.org/wiki/Main_Page) 

We recommend to install Pybel/OpenBabel with the conda package manager,
because then no further installations are necessary.
For this just type 
    
    $ conda install -c openbabel openbabel

in your command line.
(This requires of course that conda is installed)


For using the database interface from each directory, add the pythonpath to your bashrc:

```bash
export PYTHONPATH=$PYTHONPATH:<path_to_the_repository>/solvation-free-energy-database/database
``` 
  
# Usage of the database interface:
Open a Python prompt and type in:

```python
>>> from databasehandler import Database
>>> d = Database()
```

With this command the database gets initialized.
Now all solvents or solutes can be listed with:

```python
>>> print d.solvents
>>> print d.solutes
```

One can filter for one solvent, solute or both:

```python
>>> d.filtering(solvent='chloroform')
>>> d.filtering(solute='hexane')
>>> d.filtering(solvent='chloroform', solute='hexane')
{'deltaG_solv': -0.16979, 'logK': 2.87}

```
The output of the first lines is a dictionary with molecule specific propteries
and a list with all corresponding solutes and solvents, respectively.
The output of the latter, is a dictionary with the logK value 
and the solvation free energy calculated from it.
The unit of the solvation free energy is per default eV,
but one can switch to kcal/mol with

```python
>>> d.energy_unit='kcal/mol'
{'deltaG_solv': -3.92, 'logK': 2.87}
```

Adding new properties can simple done with

```python
>>> d.add_property(<molecule>, <property_name>, <property_value>)
```    
    
For adding dielectric constants, there exists a extra command:

```python
>>> d.add_dielectric_constant(<molecule>, <value>, <source>)
```

The database interface also provides the possibility to draw depictions of the molecules:

```python
>>> d.draw(<molecule>)
```

This method take the keywords 'save=True', for saving a png figure of the depiction and
'notebook=True' for showing the depiction inline in jupyter notebooks.

# Known Issues
There is a known issue in Pybel in the drawing functionality, which was already reported and is also fixed in the Pybel source code.
But it could be that with installing Pybel with conda or PIP and not directly from source the following error still appears:

```
ImportError: Tkinter or Python Imaging Library not found, but is required for image display. See installation instructions for more information.
```

For this we provided a fixed version of the pybel.py source code. Just put the Pythonpath to the directory fixed_pybel in your bashrc:

```bash
export PYTHONPATH=$PYTHONPATH:<path_to_the_repository>/solvation-free-energy-database/fixed_pybel
``` 

