#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Stand-alone script to fetch CIAAW/IUPAC data on isotope mass and abundance.
    This script will recreate the periodic_table database used by the rest of the
    interference_calculator package.
"""

import pandas as pd
import requests

header = '''# Periodic table
#
# The data in this file are published by the International Union of Pure and Applied Chemistry (IUPAC),
# Commission on Isotopic Abundances and Atomic Weights (CIAAW). The data were downloaded from the
# following sources. See links for more information and references.
#
# * atomic masses
#     publication: Wang et al, The Ame2012 atomic mass evaluation, Chinese Physics C, 2012, 36, 1603-2014
#     doi: https://doi.org/10.1088/1674-1137/36/12/003
#     data url: {}
#
# * isotopic abundances:
#     publication: Meija et al, Isotopic compositions of the elements 2013, Pure and Applied Chemistry, 2016, 88, 293-306 (table 1)
#     doi: https://doi.org/10.1515/pac-2015-0503
#     data url: {}
#
'''

mass_url = 'http://ciaaw.org/atomic-masses.htm'
abun_url = ' https://www.degruyter.com/table/j/pac.2016.88.issue-3/pac-2015-0503/pac-2015-0503.xml?id=j_pac-2015-0503_tab_001'
output = 'periodic_table.csv'

### Mass table

req = requests.get(mass_url)
mass = pd.read_html(req.text, encoding=req.encoding)[0]
mass = mass.drop(mass.index[-1])

# HTML table has rowspans, read_html does not handle it correctly.
# First 3 columns should be empty (NaN) for minor isotopes of the
# same parent element, but A and mass are in columns 0 and 1, resp.
# Split into two based on symbol == NaN, reorganize, concat back together.
mass.columns = ['atomic number', 'element', 'element name', 'atomic mass', 'mass']

partA = mass[mass['element name'].isnull()]
partA = partA[['element name', 'atomic mass', 'mass', 'atomic number', 'element']]
partA.columns = ['atomic number', 'element', 'element name', 'atomic mass', 'mass']

partB = mass[mass['element name'].notnull()]
mass = pd.concat([partA, partB]).sort_index()
mass = mass.fillna(method='pad')

mass['atomic number'] = pd.to_numeric(mass['atomic number'])
mass['atomic mass'] = pd.to_numeric(mass['atomic mass'].str.strip('*'))
# \xa0 is utf-8 encoded non-breaking space
mass['mass'] = pd.to_numeric(mass['mass'].str.split('(').str[0].str.replace('\xa0', ''))

### Abundance table

req = requests.get(abun_url)
abun = pd.read_html(req.text, encoding=req.encoding)[0]

abun.columns = ['atomic number', 'element', 'atomic mass', 'interval', 'annotation', 'abundance', 'reference', 'standard', 'interval2', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
abun = abun[['atomic number', 'element', 'atomic mass', 'abundance', 'standard']]

# No data for Po, At, Rn, Fr, Ra, Ac, also missing from mass table.
abun = abun.drop(abun[abun['element'].isin(['Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac'])].index)
abun.index = range(abun.shape[0])
# No data for Tc and Pm, but want to keep.
idx = abun[abun['element'] == 'Tc'].index
abun.loc[idx] = [43, 'Tc', 98, 0, '']
idx = abun[abun['element'] == 'Pm'].index
abun.loc[idx] = [61, 'Pm', 145, 0, '']

abun = abun.fillna(method='pad')
abun['atomic number'] = pd.to_numeric(abun['atomic number'])
abun['atomic mass'] = pd.to_numeric(abun['atomic mass'])
abun['abundance'] = pd.to_numeric(abun['abundance'].str.split('(').str[0].str.replace(' ', ''))
# \xe2\x80\x93 is utf-8 encoded en-dash
abun['standard'] = abun['standard'].str.strip('*').str.replace(b'\xe2\x80\x93'.decode('utf-8'), '')

# U233 missing, but listed in mass data, add.
u = abun.iloc[-1].copy()
u['atomic mass'] = 233
u['abundance'] = 0
abun = abun.append(u)
abun = abun.sort_values(['atomic number', 'atomic mass'])
abun.index = range(abun.shape[0])

### Merge

# Before merging, check that index, symbol, Z, and A are same.
if not mass.shape[0] == abun.shape[0]:
    raise ValueError('Mass and abundance tables have different length.')
if not (mass.index == abun.index).all():
    raise ValueError('Indices are not the same while merging mass and abundance tables.')
if not (mass['atomic number'] == abun['atomic number']).all():
    raise ValueError('Atomic number (Z) not same for all entries while merging mass and abundance tables.')
if not (mass['atomic mass'] == abun['atomic mass']).all():
    raise ValueError('Atomic mass (A) not same for all entries while merging mass and abundance tables.')
if not (mass['element'] == abun['element']).all():
    raise ValueError('Element symbols are not same for all entries while merging mass and abundance tables.')

mass['abundance'] = abun['abundance']
mass['standard'] = abun['standard']

with open(output, mode='wt', encoding='utf-8') as fh:
    fh.write(header.format(mass_url, abun_url))
    mass.to_csv(fh, index=False)
