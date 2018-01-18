#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Stand-alone script to fetch CIAAW/IUPAC data on isotope mass and abundance.
    This script will recreate the periodic_table database used by the rest of the
    interference_calculator package.
"""

import pandas as pd
import requests

mass_url = 'http://ciaaw.org/atomic-masses.htm'
abun_url = ' https://www.degruyter.com/table/j/pac.2016.88.issue-3/pac-2015-0503/pac-2015-0503.xml?id=j_pac-2015-0503_tab_001'
output = 'periodic_table.csv'

# Set to True to read from local files, useful for debugging.
_debug = False

### Mass table

if _debug:
    mass = pd.read_html('devel/mass_ciaaw.html', encoding='utf-8')[0]
else:
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

# Add isotope column
atomic_mass = mass['atomic mass'].values
element = mass['element'].values
isotope = [str(am) + el for am, el in zip(atomic_mass, element)]
mass['isotope'] = isotope

### Abundance table

if _debug:
    abun = pd.read_html('devel/abun_ciaaw.html', encoding='utf-8')[0]
else:
    req = requests.get(abun_url)
    abun = pd.read_html(req.text, encoding=req.encoding)[0]

abun.columns = ['atomic number', 'element', 'atomic mass', 'interval', 'annotation', 'abundance', 'reference', 'standard', 'interval2']
abun = abun[['atomic number', 'element', 'atomic mass', 'abundance', 'standard']]

# No data for Po, At, Rn, Fr, Ra, Ac, also missing from mass table.
abun = abun.drop(abun[abun['element'].isin(['Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac'])].index)
abun.index = range(abun.shape[0])
# No data for Tc and Pm, but want to keep.
idx = abun[abun['element'] == 'Tc'].index
abun.loc[idx] = [43, 'Tc', 98, '0.0', '']
idx = abun[abun['element'] == 'Pm'].index
abun.loc[idx] = [61, 'Pm', 145, '0.0', '']

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

### Major isotope
# For each element, determine the major isotope, the isotope with the highest abundance.
elements = mass['element'].unique()
major_isotope = []

for el in elements:
    el_slice = mass[mass['element'] == el]
    major_mass = el_slice.sort_values('abundance', ascending=False).iloc[0].loc['atomic mass']
    number_of_isotopes = el_slice.shape[0]
    major_isotope.extend([str(major_mass) + el] * number_of_isotopes)

mass['major isotope'] = major_isotope

# Reorder columns
mass = mass[['atomic number', 'element', 'element name', 'major isotope',
             'isotope', 'atomic mass', 'mass', 'abundance', 'standard']]

with open(output, mode='wt', encoding='utf-8') as fh:
    mass.to_csv(fh, index=False)
