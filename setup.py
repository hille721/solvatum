#!/usr/bin/env python
#
# Copyright (C) 2018 Christoph Hille

DESCRIPTION = "Solv@TUM - The Solvation Free Energy Database"
LONG_DESCRIPTION = """\
With this module you can interact with the Solv@TUM database. 
The usage of this database interface is described in the README. 

Please note: This module was written for Python=2.7. 
Currently there is no compatibility with Python=3.*.
"""

DISTNAME = 'solvatum'
MAINTAINER = 'Christoph Hille'
MAINTAINER_EMAIL = 'c.hille@tum.de'
LICENSE = 'CC-BY-SA'
DOWNLOAD_URL = 'https://mediatum.ub.tum.de/1452573'
VERSION = '1.0.0'

INSTALL_REQUIRES = ['openbabel>=2.4.0']

PACKAGES = ['solvatum']
        
import sys
from setuptools import setup

if __name__ ==  "__main__":
    try:
        import pybel
    except ImportError:
        sys.exit(u"You have to install the Open Babel chemistry library and the Python interface 'pybel' before you can install the Solv@TUM database interface. See the README for more informations.") 
        
    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        packages=PACKAGES,
        include_package_data=True,
        python_requires='>=2.7, !=3.*',
        )
