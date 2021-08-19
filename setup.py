#!/usr/bin/env python
#
# Copyright (C) 2018 Christoph Hille

DESCRIPTION = "Solv@TUM - The Solvation Free Energy Database"
LONG_DESCRIPTION = """\
With this module you can interact with the Solv@TUM database.
The usage of this database interface is described in the README.
"""

DISTNAME = 'solvatum'
MAINTAINER = 'Christoph Hille'
MAINTAINER_EMAIL = 'c.hille@tum.de'
LICENSE = 'CC-BY-SA'
URL = 'https://github.com/hille721/solvatum'
VERSION = '1.1.0'

PACKAGES = ['solvatum']

import sys
from setuptools import setup

if __name__ ==  "__main__":
    try:
        from openbabel import pybel
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
        url=URL,
        packages=PACKAGES,
        include_package_data=True,
        python_requires='>=3.6,<3.8.0a0',
        )
