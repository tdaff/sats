#!/usr/bin/env python

from setuptools import setup
from glob import glob

scripts = glob('bin/*')

setup(name='sats',
      version='0.1',
      description='Automatic structure generation for forcefield training',
      long_description=open('README.rst').read(),
      author='Tom Daff',
      author_email='tdd20@cam.ac.uk',
      license='BSD',
      url='http://bitbucket.org/tdaff/sats/',
      packages=['sats'],
      scripts=scripts,
      requires=['ase', 'quippy', 'numpy'],
      classifiers=["Programming Language :: Python",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3",
                   "Development Status :: 3 - Alpha",
                   "Intended Audience :: Science/Research",
                   "Intended Audience :: System Administrators",
                   "License :: OSI Approved :: BSD License",
                   "Operating System :: OS Independent",
                   "Topic :: System :: Monitoring"])
