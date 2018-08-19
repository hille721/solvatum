#!/usr/bin/env python

from distutils.core import setup

setup(name='solvatum',
      version='1.0',
      description='Solvation Free Energy Database',
      author='Christoph Hille',
      author_email='c.hille@tum.de',
      url='https://mediatum.ub.tum.de/1452573',
      license='cc-by-sa',
      packages=['solvatum'],
      data_files=[('solvatum', ['solvatum.sdf']),
                  ('references', ['solvatum_references.bib'])],
     )
