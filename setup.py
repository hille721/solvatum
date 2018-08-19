#!/usr/bin/env python

import sys
from setuptools import setup

try:
    import pybel
except ImportError:
    sys.exit(u"You have to install the Open Babel chemistry library and the Python interface 'pybel' before you can install the Solv@TUM database interface. See the README for more informations.")

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='solvatum',
    version='1.0',
    packages=['solvatum'],
    include_package_data=True,
    author='Christoph Hille',
    author_email='c.hille@tum.de',
    description='Solvation Free Energy Database',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://mediatum.ub.tum.de/1452573',
    license='CC-BY-SA',
)
