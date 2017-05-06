#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Calculate isotopic interference and standard ratios. """
import interference_calculator as ic
from pandas import concat

__all__ = ['interference', 'standard_ratio']

def interference(atoms, mz, mzrange=0.3, maxsize=5, charge=[1],
                chargesign='-', formula_style='html'):
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
        electrons (+ charge) or extra electrons (- charge).

        Molecular formulas are represented in the style in formula_style (default
        is 'html'). See Molecule() for more options.

        Returns a pandas.DataFrame with a column 'molecule' with molecular formula,
        a column 'mass/charge' for the mass-to-charge ratio, a column 'mass/charge diff'
        for the mass/charge difference between this ion and the target ion, a column
        'MRP' which gives the mass-resolving power (mz/Δmz) needed to resolve the
        target ion from this ion, a column 'charge', and a column 'probability',
        which for now is simply the product of the abundances of the isotopes in
        the given molecular ion.
    """
    picked_atoms = []
    for a in atoms:
        picked_atoms.append(ic.isotope_data[ic.isotope_data.element == a])

    picked_atoms = concat(picked_atoms).query('index not in ("H", "D", "T")')

    # Mass-to-charge can be given as either a number, or as a molecule (string).
    # Calculate m/z from molecular formula.
    if mz:
        try:
            mz = float(mz)
        except ValueError:
            molec = ic.Molecule(mz)
            ch = molec.charge
            if ch == 0:
                ch = 1
            mz = molec.mass / ch

    # Create a list with all possible combinations up to maxsize atoms
    # Create same list for masses and abundances; combos are created in same order.
    isotope_combos = []
    mass_combos = []
    abundance_combos = []
    for size in range(1, maxsize + 1):
        i = itertools.combinations_with_replacement(picked_atoms.index, size)
        m = itertools.combinations_with_replacement(picked_atoms.isotope_mass, size)
        a = itertools.combinations_with_replacement(picked_atoms.abundance, size)
        isotope_combos += list(i)
        mass_combos += list(m)
        abundance_combos += list(a)

    # Calculate total mass (sum), probability (product, below)
    masses = DataFrame(mass_combos).sum(axis=1)
    abundances = DataFrame(abundance_combos)

    # There can be missing values because: (A) missing from NIST
    # table, (B) because combos were of different length, and turning
    # them into a square array fills up each row with NaNs.
    # Replaced A-type NaN with -1 or minabundance earlier, replace
    # B-type NaN with 1 here (will be eaten by product), then
    # -1 back to NaN (nothing should be < 0 if minabundance was given).
    print(abundances)
    # abundances = abundances.fillna(value=1)
    # abundances = abundances.where(abundances >= 0)
    abundances = abundances.prod(axis=1, skipna=False)
    print(abundances)

    # Using Molecule() to convert atom list to molecular formula
    # is too slow for long list. Simple join for now.
    molecules = [' '.join(m) for m in isotope_combos]

    data = DataFrame({'molecule': molecules,
                      'mass/charge': masses,
                      'probability': abundances})

    if len(charge) == 1:
        if charge[0] == 1:
            charge_str = [''] * len(data.index)
        else:
            charge_str = [str(charge[0])] * len(data.index)
        data['charge'] = charge[0]
    else:
        cc = []
        charge_str = []
        for c in charge:
            cc.extend([c] * len(data.index))
            if c == 1:
                charge_str.extend([''] * len(data.index))
            else:
                charge_str.extend([str(c)] * len(data.index))

        # repeat data once for each charge
        data = concat([data] * len(charge))
        data['charge'] = cc

    charge_str = [' [' + c + chargesign + ']' for c in charge_str]
    data['molecule'] += charge_str

    # Correct mass for charge * electron mass
    if chargesign == '+':
        data['mass/charge'] -= data['charge'] * mass_electron
    elif chargesign == '-':
        data['mass/charge'] += data['charge'] * mass_electron
    else:
        raise ValueError('chargesign must be either "+" or "-".')

    data['mass/charge'] /= data['charge']
    if mz:
        data['mass/charge diff'] = data['mass/charge'] - mz
        data['MRP'] = mz/data['mass/charge diff'].abs()
    else:
        data['mass/charge diff'] = 0
        data['MRP'] = 0

    if mz:
        results = data.loc[(data['mass/charge'] >= mz - mzrange)
                         & (data['mass/charge'] <= mz + mzrange)]
    else:
        results = data

    def convert(s):
        m = ic.Molecule(s)
        return m.formula(style=formula_style)

    # Next line gives SettingWithCopyWarning no matter how I do it.
    results.loc[:, 'molecule'] = results['molecule'].apply(convert)
    results = results.ix[:, ['molecule',
                             'charge',
                             'mass/charge',
                             'mass/charge diff',
                             'MRP',
                             'probability']]

    return results.sort_values('mass/charge')


def standard_ratio(atoms, formula_style='html'):
    """ Give the stable isotope and solar abundance for the given element(s). """
    data = []
    for a in atoms:
        eldata = ic.isotope_data[ic.isotope_data['element'] == a]
        highest = eldata['abundance'].max()
        eldata['inverse'] = highest/eldata['abundance']
        data.append(eldata)

    results = concat(data)

    def convert(s):
        m = ic.Molecule(s)
        return m.formula(style=formula_style, show_charge=False)

    results.index = [convert(i) for i in results.index]
    results = results.ix[:, ['isotope_mass', 'abundance', 'inverse']]

    return results

if __name__ == '__main__':
    print('Possible interferences for MgO given as sample consisting of H, C, O, N, and Si:')
    print(interference(['H', 'C', 'O', 'N', 'Si'], 'MgO'))
    print('Standard ratios of C and S:')
    print(standard_ratios(['C', 'S']))
