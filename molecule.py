#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Molecule() is a class that takes an input string of a chemical formula,
    parses the string into atomic units, and stores relevant molecular data.
    The chemical formula can be output in a number of ways, including custom
    formatting using simple templates.
"""
import interference_calculator as ic
from collections import Counter
from numpy import prod
from math import factorial
# from scipy.misc import factorial
import re

__all__ = ['Molecule', 'templates', 'html_template', 'latex_template', 'isotope_template']

# How to determine the input notation?
# If it contains a separation marker (space or anything not letter, numbers, or []+-),
# then it is isotope notation, otherwise empirical.
# However, treat a single element (one capital in string) as isotope notation.
_isotope_notation_rx = re.compile(r'(\d*)([A-Z][a-z]{0,2}|[+-])(\d*)')
_empirical_notation_rx = re.compile(r'(?:\[(\d*)\])?([A-Z][a-z]{0,2}|\[\d*[+-]\])(\d*)')
_is_empirical_rx = re.compile(r'^[A-Za-z\d\[\]+-]+$')
_is_single_rx = re.compile(r'[A-Z]')

templates = ['html_template', 'latex_template', 'isotope_template']

html_template = {
    'isotope': '<sup>{}</sup>',
    'element': '{}',
    'stoich': '<sub>{}</sub>',
    'charge': '<sup>{}</sup>',
    'minorjoin': '',
    'majorjoin': ''
}

latex_template = {
    'isotope': '{{}}^{{{}}}',
    'element': '{{{}}}',
    'stoich': '_{{{}}}',
    'charge': '{{}}^{{{}}}',
    'minorjoin': '',
    'majorjoin': ''
}

isotope_template = {
    'isotope': '{}',
    'element': '{}',
    'stoich': '{}',
    'charge': '{}',
    'minorjoin': '',
    'majorjoin': ' '
}

class FormatError(Exception):
    """ Raised when an error in the formatting of molecular formula is detected. """
    pass


class Molecule(object):
    """ Represents a molecule or molecular ion. """
    def __init__(self, molecule):
        """ Parses a chemical formula string and returns an object that
            holds properties of the molecule or molecular ion.

            Two forms of input string are supported: isotope notation and
            empirical formula notation. These names and notations are used
            for input and output.

            Isotope notation is a separated list of elements, where each
            element is an isotope of the form NXxxn, with N the isotope
            number, Xxx is the element name, and n is the stoichiometric
            number. Any character except A-Z, 0-9, +, -, [, or ] may be
            used to separate the elements. If no isotope number is specified,
            the most common isotope is assumed (e.g. C -> 12C). A charge
            may optionally be given as the last element. This form is
            useful for inputting many unusual isotopes.

            Isotope notation: '12C2 15N O3 2-'

            Emperical formula notation is a form of shorthand. It contains
            no spaces and no isotope numbers, only stoichiometric numbers.
            Isotopes can optionally be given, surrounded by []. A charge
            may optionally be given at the end, also surrounded by [].
            This form is useful for inputting larger molecules with few
            unusual isotopes.

            Empirical formula notation: 'HCOOCH2[15]NH3[2-]'

            D and T are accepted as aliases for 2H and 3H, but are internally
            converted to 2H and 3H. See Molecule.formula() for output options.

            After parsing the input formula, the compositional information
            is stored under Molecule.composition as a pandas.DataFrame. This
            is a slice (copy) of the MassTable, with an extra column containing
            the stoichiometry. Other information stored in the Molecule class
            is Molecule.mass, Molecule.charge, and Molecule.chargesign.
        """
        self.input = molecule
        self.mass = 0.0
        self.abundance = 1.0
        self.charge = 0
        self.chargesign = ''

        self.parse(molecule)

    def __str__(self):
        return self.input + ' --> ' + self.formula()

    def parse(self, molecule, _mass_table=None):
        """ Parse input, retrieve elements. """
        # TODO: convert to using pyparsing, easier for future maintenance.
        molecule = molecule.strip()

        # Determine notation
        if (len(re.findall(_is_single_rx, molecule)) == 1 or
            not re.match(_is_empirical_rx, molecule)):
                units = re.findall(_isotope_notation_rx, molecule)
        else:
            units = re.findall(_empirical_notation_rx, molecule)

        # Check for charge and sign in input
        if '[' in units[-1][1]:
            chsgn = units.pop(-1)[1].strip(' []')
            charge = chsgn.strip('+-')
        elif '+' in units[-1][1] or '-' in units[-1][1]:
            charge, chsgn, skip = units.pop(-1)
        else:
            charge = 0
            chsgn = ''

        if charge == '':
            self.charge = 0
        else:
            self.charge = int(charge)

        if '+' in chsgn:
            self.chargesign = '+'
        else:
            self.chargesign = '-'

        # Fix certain shorthand notations
        iunits = []
        for u in units:
            # u = ['12', 'C', '2']
            u = list(u)

            # D -> 2H
            if u[1] == 'H' and not u[0]:
                u[0] = '1'
            elif u[1] == 'D':
                u[0] = '2'
                u[1] = 'H'
            elif u[1] == 'T':
                u[0] = '3'
                u[1] = 'H'

            # C -> 12C
            if u[0] == '':
                idx = ic.elements.index(u[1])
                mi = ic.main_isotopes[idx]
                u[0] = mi.strip(u[1])

            # Stoichiometry
            try:
                u[2] = int(u[2])
            except ValueError:
                if u[2] == '':
                    u[2] = 1
                else:
                    msg = 'Stoichiometry not a number: {}'
                    raise FormatError(msg.format(''.join(u)))

            # Isotope from element and isotope number.
            # Repeat each isotope stoich. times, so Counter can count all.
            iunits += [str(u[0]) + u[1]] * u[2]

        c = Counter(iunits)
        self.isotopes = list(c.keys())
        # TODO: make Molecule sortable, use self.sort()
        self.isotopes.sort()
        self.isotopes = tuple(self.isotopes)

        # Don't use c.values, isot/keys are sorted
        self.stoichiometry = tuple(c[i] for i in self.isotopes)
        self.indices = tuple(ic.isotopes.index(i) for i in self.isotopes)
        self.elements = tuple(ic.elements[i] for i in self.indices)
        self.atomic_numbers = tuple(ic.atomic_numbers[i] for i in self.indices)
        self.unit_masses = tuple(ic.unit_masses[i] for i in self.indices)
        self.isotope_masses = tuple(ic.isotope_masses[i] for i in self.indices)
        self.abundances = tuple(ic.abundances[i] for i in self.indices)

        for m, s in zip(self.isotope_masses, self.stoichiometry):
            self.mass += m * s

        # call abun func here


    def sort(self):
        """ Sort elements from heavy to light, but keep same elements together. """
        raise NotImplementedError('not yet')

    def formula(self, style='', HtoD=True, show_charge=True, template={}):
        """ Return the molecular formula as a string.

            The molecular formula can be formatted as html
            (style='html'), LaTeX (style='latex'), plain
            text isotope notation (style='isotope', default),
            plain text empirical formula (style='empirical'),
            or in a custom format (style='custom'), see below.

            1H, 2H, 3H will be converted to H, D, and T; set
            HtoD=False to output as 1H, 2H, and 3H instead.

            Charge and sign will be automatically added, unless
            show_charge is set to False.

            If style='custom', a custom template can be used to
            format the molecular formula. The template must be
            a dict containing 6 keys: isotope, element, stoich,
            charge, minorjoin, and majorjoin. The isotope,
            element, stoich, and charge keys should refer to a
            template string containing a curly bracket pair,
            which will be replaced using string.format(). The
            minorjoin string will be used to connect all the
            parts into one unit, while the majorjoin string
            will be used to connect all the units into the
            final output string.
        """
        if show_charge:
            if self.charge == 0:
                charge = ''
            elif self.charge == 1:
                charge = self.chargesign
            else:
                charge = str(self.charge) + self.chargesign
        else:
            charge = ''

        elem = list(self.elements)
        umass = [str(u) for u in self.unit_masses]
        stoich = [str(s) if s > 1 else '' for s in self.stoichiometry]

        if HtoD:
            for n, x in enumerate(zip(umass, elem)):
                um, el = x
                if el == 'H':
                    if um == '1':
                        umass[n] = ''
                    elif um == '2':
                        umass[n] = ''
                        elem[n] = 'D'
                    elif um == '3':
                        umass[n] = ''
                        elem[n] = 'T'

        if style == 'html':
            templ = html_template
        elif style == 'latex':
            templ = latex_template
        elif style == 'empirical':
            raise NotImplementedError('not yet')
        elif style == 'custom':
            templ = template
        else:
            templ = isotope_template

        molecule = []
        for um, el, st in zip(umass, elem, stoich):
            if um:
                um_str = templ['isotope'].format(um)
            else:
                um_str = ''
            el_str = templ['element'].format(el)
            if st:
                st_str = templ['stoich'].format(st)
            else:
                st_str = ''
            m = templ['minorjoin'].join((um_str, el_str, st_str))
            molecule.append(m)
        if charge:
            molecule.append(templ['charge'].format(charge))
        return templ['majorjoin'].join(molecule)

    def relative_abundance(isotopes, stoichometry):
        """ Given a list of isotopes and a list of stoichiometric numbers
            calculate relative abundance for entire molecule. """
        # multiple isotopes e.g. 28Si (92.2%) 29Si (4.7%) 30Si (3.1%)
        # In this type of mass spectorscopy we only look at total mass of molecule,
        # not position of isotope. Therefore Si4-29Si has 5 isobaric structures:
        #   Si4-29Si, Si3-29Si-Si, Si2-29Si-Si2, Si-29Si-Si3, 29Si-Si4
        #
        # Same problem as drawing 3 green + 2 red balls from a bag of multi-coloured
        # balls. Calculate the probability mass function of multinomial distribution
        # with replacement. Replacement condition is acceptable, since pool from
        # which isotopes are pulled from which molecule is made (nature, sample
        # in mass spectrometer during sputtering) is large compared to number of 
        # molecules being formed and therefore composition of pool does not change
        # significantly (assuming homogenous distribution).
        # 
        #   f(xi, n, pi) = n!/(x1!*x2!*...xk!) * (p1**x1 * p2**x2 * ... pk**xk)
        #   for i = 1..k
        # with
        #   n = total number of all isotopes from the same parent element
        #     = sum(xi) for i=1..k
        #   k = number of different isotopes in molecule
        #   xi = number of isotope i = stoichiometry
        #   pi = probability of isotope i = natural abundance
        #
        # Example: molecule 12C 16O2 18O
        # C is independent of O
        # there are 3 O in the molecule, n = 3
        # there are 2 O isotopes in this molecule, k = 2
        # for 16O: xi = 2, for 18O: xi = 1
        # for 16O: pi = 0.9976 for 18O: pi = 0.002 (and 0.0004 for 17O)

        data = ic.isotope_data.loc[isotopes]
        data['stoich'] = stoichiometry

        parents = data['main_isotope'].value_counts().to_dict()
        abun_per_el = []
        for el, k in parents.items():
            d = data[data['main_isotope'] == el]
            n = d['stoich'].sum()

            if k == 1:
                # Simple case of single isotope, even if it occurs n times
                abun = d['abundance'].iat[0] ** n
            else:
                abun = factorial(n)/factorial(d['stoich']).prod() * (d['abundance'] ** d['stoich']).prod()

            abun_per_el.append(abun)
        return prod(abun_per_el)


if __name__ == '__main__':
    print(Molecule('12C2 15N H4 2+').formula(style='html'))