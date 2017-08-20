#!/usr/bin/env python
""" Distutils setup file for interference calculator. """
from setuptools import setup

# from distutils.core import setup
# import py2exe
#
# setup(console=['hello.py'])

# write test to check for either PyQt4 or PyQt5 (PySide not tested)


setup(
    # metadata
    name = 'interference calculator',
    version = '1.0',
    description = 'A tool to calculate mass interference based on a target material and a target mass. Useful for mass spectroscopy (SIMS, TIMS, ICP-MS) of elements and small molecules and their isotopes. The program can also display the standard ratio of all isotopes based on a given standard.',
    url = 'https://github.com/zanpeeters/interference_calculator',
    author = 'Zan Peeters',
    license = 'BSD 3-Clause Clear',
    classifiers = [
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering'
    ],
    keywords = 'interference, mass spectroscopy, isotope, element, standard ratio',

    # package info
    install_requires = ['pandas', 'PyQt4'],
    package_data = {'': ['periodic_table.csv']}
    entry_points = {
        'gui_scripts': ['interference_calculator = interference_calculator_ui:run']
    }
)
