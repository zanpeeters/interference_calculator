# -*- coding: utf-8 -*-
""" Calculate isotopic interference and standard ratios. """

import pandas as pd
import itertools
from interference_calculator.molecule import Molecule, mass_electron, periodic_table

__all__ = ['interference', 'standard_ratio']

def interference(atoms, target, targetrange=0.3, maxsize=5, charge=[1],
                 chargesign='-', style='plain'):
    """ For a list of atoms (the composition of the sample),
        calculate all molecules that can be formed from a
        combination of those atoms (the interferences),
        including all stable isotopes, up to maxsize atoms,
        that have a mass-to-charge ratio within target ± targetrange.

        The target can be given as a mass-to-charge ratio or as a
        molecular formula. Molecular formulas are interpreted by Molecule().
        See Molecule() docstring for a detailed explanation on how to enter
        molecular formulas. If target is None, no filtering will be done and
        all possible combinations of all atoms and isotopes up to maxsize
        length will be calculated. Target information will be added to the
        output, unless target is None.

        Charge is usually 1, irrespective of sign. Give charge = [1, 2, 3]
        to also include interferences with higher charges. Masses are 
        adjusted for missing electrons (+ charge), extra electrons (- charge),
        or not adjusted (o charge, lower-case letter O). Setting charge=0
        has the same effect as setting chargesign='o'. The charge for the
        target ion, if target is specified as molecule instead of a number,
        can be different from the charge on the interferences. If no charge is
        specified for the target, the first charge and the chargesign of the
        interferences are used for the target.

        Molecular formulas are formatted in style (default is 'plain').
        See Molecule() for more options.

        Returns a pandas.DataFrame with a column 'molecule' with molecular formula,
        a column 'charge', a column 'mass/charge' for the mass-to-charge ratio, a
        column 'mass/charge diff' for the mass/charge difference between this ion
        and the target mass/charge, a column 'MRP' which gives the mass-resolving
        power (mz/Δmz) needed to resolve this ion from the target ion, a column
        'target', which indicates whether this row was specified as the target,
        and a column 'probability', which gives the combinatorial probability of
        encoutering this combination of isotopes, given the composition of the
        sample and the natural abundances of the isotopes.
    """
    if isinstance(charge, (int, float, str)):
        charge = tuple(int(charge))
    elif isinstance(charge, (tuple, list)):
        charge = tuple(int(c) for c in charge)
    else:
        raise ValueError('charge must be given as a number or a list of numbers.')

    if chargesign not in ('+', '-', 'o', '0'):
        raise ValueError('chargesign must be either "+", "-", "o", or "0".')

    # How to handle charge?
    # 1. charge for interferences
    #       - can be multiple values
    #       - specified by parameter
    # 2. charge for target
    #       - only one value
    #       - can be different from 1
    #       - must be specified in target formula
    #       - if unspecified, take sign and first value from 1
    if target:
        try:
            target_mz = float(target)
            target = str(target)
            target_charge = 0
            target_chargesign = 'o'
            target_abun = 1
        except ValueError:
            m = Molecule(target)
            if m.chargesign:
                target_chargesign = m.chargesign
            else:
                target_chargesign = chargesign
            if m.charge:
                target_charge = m.charge
            else:
                target_charge = charge[0]
            target_mz = m.mass
            target_abun = m.abundance
            if m.charge > 0:
                # mass correction done in Molecule.parse()
                target_mz /= m.charge
    else:
        target_mz = 0
        target_charge = 0
        target_chargesign = '0'
        target_abun = 0

    # Retrieve info from perioic table for all atoms in sample.
    # Create a list with all possible combinations up to maxsize atoms.
    # Create same list for masses, combos are created in same order.
    picked_atoms = periodic_table[periodic_table['element'].isin(atoms)]
    isotope_combos = []
    mass_combos = []
    for size in range(1, maxsize + 1):
        i = itertools.combinations_with_replacement(picked_atoms['isotope'], size)
        m = itertools.combinations_with_replacement(picked_atoms['mass'], size)
        isotope_combos.extend(list(i))
        mass_combos.extend(list(m))

    masses = pd.DataFrame(mass_combos).sum(axis=1)
    molecules = [' '.join(m) for m in isotope_combos]
    data = pd.DataFrame({'molecule': molecules,
                         'mass/charge': masses})

    # ignore charge(s) for sign o
    if chargesign in ('o', '0'):
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

    if target:
        data = data.loc[(data['mass/charge'] >= target_mz - targetrange)
                      & (data['mass/charge'] <= target_mz + targetrange)]
        data['mass/charge diff'] = data['mass/charge'] - target_mz
        data['MRP'] = target_mz/data['mass/charge diff'].abs()
    else:
        data['mass/charge diff'] = 0.0
        data['MRP'] = pd.np.inf

    molec = []
    abun = []
    for molecule in data['molecule'].values:
        m = Molecule(molecule)
        abun.append(m.abundance)
        molec.append(m.formula(style=style))

    data['molecule'] = molec
    data['probability'] = abun
    data['target'] = False
    target_data = {
        'molecule': target,
        'charge': target_charge,
        'mass/charge': target_mz,
        'mass/charge diff': 0,
        'MRP': pd.np.inf,
        'probability': target_abun,
        'target': True
    }
    data = data.append(target_data, ignore_index=True)
    return data[['molecule', 'charge', 'mass/charge',
                 'mass/charge diff', 'MRP', 'probability', 'target']]

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
