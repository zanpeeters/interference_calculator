#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Calculate isotopic interference and standard ratios. """

import pandas as pd
import itertools
from .molecule import Molecule, mass_electron, periodic_table

__all__ = ['interference', 'standard_ratio']

def interference(atoms, mz, mzrange=0.3, maxsize=5, charge=[1],
                 chargesign='-', style='html'):
    """ For a list of atoms, calculate all molecular ions that can be formed
        from those atoms, including all isotopes, up to maxsize atoms,
        that have a mass-to-charge ratio within mz +/- mzrange.

        The mass-to-charge ratio (mz) can be given as a number or as a
        molecular formula. Molecular formulas are interpreted by Molecule().
        See Molecule() docstring for a detailed explanation on how to enter
        molecular formulas. If mz is None, no filtering will be done and
        all possible combinations of all isotopes up to maxsize length will
        be calculated.

        Charge is usually 1, irrespective of sign. Give charge = [1, 2, 3]
        to also include higher charged ions. Masses are adjusted for missing
        electrons (+ charge), extra electrons (- charge), or not adjusted
        (o charge, lower-case letter O). Setting charge=0 has the same
        effect as setting chargesign='o'.

        Molecular formulas are formatted in style (default is 'html').
        See Molecule() for more options.

        Returns a pandas.DataFrame with a column 'molecule' with molecular formula,
        a column 'charge', a column 'mass/charge' for the mass-to-charge ratio, a
        column 'mass/charge diff' for the mass/charge difference between this ion
        and the target mass/charge, a column 'MRP' which gives the mass-resolving
        power (mz/Δmz) needed to resolve this ion from the target ion, and a column
        'probability', which gives the combinatorial probability of encoutering this
        combination of isotopes, given the composition of the sample.
    """
    if isinstance(charge, int):
        charge = [charge]
    if chargesign not in ('+', '-', 'o', 'O', '0'):
        raise ValueError('chargesign must be either "+", "-", or "o".')

    picked_atoms = periodic_table[periodic_table['element'].isin(atoms)]

    # Mass-to-charge can be given as either a number, or as a molecule (string).
    # Calculate m/z from molecular formula.
    if mz:
        try:
            mz = float(mz)
        except ValueError:
            m = Molecule(mz)
            mz = m.mass
            if m.charge > 0:
                mz /= m.charge

    # Create a list with all possible combinations up to maxsize atoms.
    # Create same list for masses, combos are created in same order.
    isotope_combos = []
    mass_combos = []
    for size in range(1, maxsize + 1):
        i = itertools.combinations_with_replacement(picked_atoms['isotope'], size)
        m = itertools.combinations_with_replacement(picked_atoms['mass'], size)
        isotope_combos.extend(list(i))
        mass_combos.extend(list(m))

    masses = pd.DataFrame(mass_combos).sum(axis=1)

    # Using Molecule() to convert atom list to molecular formula and
    # calculate abundance is too slow for long list.
    # For 5 atoms, max size 5, 8567 combos: 80 s. Do for trimmed list later.
    molecules = [' '.join(m) for m in isotope_combos]

    data = pd.DataFrame({'molecule': molecules,
                         'mass/charge': masses})

    # ignore charge(s) for sign o
    if chargesign in ('o', 'O', '0'):
        data['charge'] = 0
    else:
        data_w_charge = []
        for ch in charge:
            d = data.copy()
            d['charge'] = ch
            if ch == 0:
                data_w_charge.append(d)
                continue
            elif ch == 1:
                charge_str = ' {}'.format(chargesign)
            else:
                charge_str = ' {}{}'.format(ch, chargesign)
            d['molecule'] += charge_str
            d['mass/charge'] /= ch
            if chargesign == '+':
                d['mass/charge'] -= mass_electron
            else:
                d['mass/charge'] += mass_electron
            data_w_charge.append(d)
        data = pd.concat(data_w_charge)

    if mz:
        data = data.loc[(data['mass/charge'] >= mz - mzrange)
                      & (data['mass/charge'] <= mz + mzrange)]
        data['mass/charge diff'] = data['mass/charge'] - mz
        data['MRP'] = mz/data['mass/charge diff'].abs()
    else:
        data['mass/charge diff'] = 0.0
        data['MRP'] = 0.0

    molec = []
    abun = []
    for molecule in data['molecule'].values:
        m = Molecule(molecule)
        abun.append(m.abundance)
        molec.append(m.formula(style=style))

    data['molecule'] = molec
    data['probability'] = abun
    return data[['molecule', 'charge', 'mass/charge',
                 'mass/charge diff', 'MRP', 'probability']]

def standard_ratio(atoms, style='html'):
    """ Give the stable isotopes and their standard abundance for the given element(s). """
    data = periodic_table[periodic_table['element'].isin(atoms)].copy()
    data['ratio'] = 1.0
    data['inverse ratio'] = 1.0
    for a in atoms:
        abun = data.loc[data['element'] == a, 'abundance'].copy()
        ratio = abun/abun.max()
        inv_ratio = 1/ratio
        data.loc[data['element'] == a, 'ratio'] = ratio
        data.loc[data['element'] == a, 'inverse ratio'] = inv_ratio

    pretty_isotopes = []
    for i in data['isotope'].values:
        m = Molecule(i)
        pretty_isotopes.append(m.formula(style=style, show_charge=False, all_isotopes=True))
    data['isotope'] = pretty_isotopes

    return data[['isotope', 'mass', 'abundance', 'ratio', 'inverse ratio', 'standard']]
