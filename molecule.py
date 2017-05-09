#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Molecule() is a class that takes an input string of a chemical formula,
    parses the string into atomic units, and stores relevant molecular data.
    The chemical formula can be output in a number of ways, including custom
    formatting using simple templates.
"""
import pandas as pd
from collections import Counter
from numpy import prod
from scipy.misc import factorial
import re

__all__ = ['Molecule', 'mass_electron', 'periodic_table',
           'templates', 'html_template', 'latex_template', 'isotope_template']

periodic_table = pd.read_csv('periodic_table.csv', comment='#')

# CODATA 2014, http://physics.nist.gov/cgi-bin/cuu/Value?me
mass_electron = 0.0005485799090

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
    'count': '<sub>{}</sub>',
    'charge': '<sup>{}</sup>',
    'minorjoin': '',
    'majorjoin': ''
}

latex_template = {
    'isotope': '{{}}^{{{}}}',
    'element': '{{{}}}',
    'count': '_{{{}}}',
    'charge': '{{}}^{{{}}}',
    'minorjoin': '',
    'majorjoin': ''
}

isotope_template = {
    'isotope': '{}',
    'element': '{}',
    'count': '{}',
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
            number, Xxx is the element name, and n is the count
            number (subscript). Any character except A-Z, 0-9, +, -, [, or ]
            may be used to separate the elements. If no isotope number is
            specified, the most common isotope is assumed (e.g. C -> 12C).
            A charge may optionally be given as the last element. This form
            is useful for inputting many unusual isotopes.

            Isotope notation: '12C2 15N O3 2-'

            Emperical formula notation is a form of shorthand. It contains
            no spaces and no isotope numbers, only count numbers.
            Isotopes can optionally be given, surrounded by []. A charge
            may optionally be given at the end, also surrounded by [].
            This form is useful for inputting larger molecules with few
            unusual isotopes.

            Empirical formula notation: 'HCOOCH2[15]NH3[2-]'

            D is an accepted alias for 2H, but is internally converted
            to 2H. See Molecule.formula() for output options.

            After parsing, relevant information about the molecule, such as
            total mass, abundance, charge, and sign, as well as lists of
            atoms, isotopes, atomic masses, and a few others are stored in
            the Molecule() object.
        """
        self.input = molecule
        self.mass = 0.0
        self.abundance = 1.0
        self.charge = 0
        self.chargesign = ''

        self.elements = []
        self.counts = []
        self.atomic_numbers = []
        self.atomic_masses = []
        self.masses = []
        self.abundances = []

        self.parse()
        self.relative_abundance()
        self.molecular_formula = self.formula()

    def __str__(self):
        return self.input + ' --> ' + self.molecular_formula

    def parse(self):
        """ Parse input, retrieve elements from periodic table,
            calculate mass and abundance.
        """
        # TODO: convert to using pyparsing, easier for future maintenance.
        self.input = self.input.strip()

        # Determine notation
        if (len(re.findall(_is_single_rx, self.input)) == 1 or
            not re.match(_is_empirical_rx, self.input)):
                units = re.findall(_isotope_notation_rx, self.input)
        else:
            units = re.findall(_empirical_notation_rx, self.input)

        # Check for charge and sign in input
        if '[' in units[-1][1]:
            chsgn = units.pop(-1)[1].strip(' []')
            charge = chsgn.strip('+-')
        elif '+' in units[-1][1] or '-' in units[-1][1]:
            charge, chsgn, skip = units.pop(-1)
        else:
            charge = 0
            chsgn = ''

        if chsgn in ('+', '-', 'o', 'O', '0', ''):
            self.chargesign = chsgn
        else:
            raise FormatError('Unknown chargesign: {}'.format(chsgn))

        if charge == '':
            if chsgn:
                self.charge = 1
            else:
                self.charge = 0
        else:
            self.charge = int(charge)

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
                mi = periodic_table[periodic_table['element'] == u[1]].iloc[0].loc['major isotope']
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
            # Repeat each isotope count times, so Counter can count all.
            iunits += [str(u[0]) + u[1]] * u[2]

        c = Counter(iunits)
        self.isotopes = sorted(c.keys())
        self.counts = [c[i] for i in self.isotopes]

        for i in self.isotopes:
            isotope = periodic_table[periodic_table['isotope'] == i].iloc[0]
            self.elements.append(isotope['element'])
            self.atomic_numbers.append(isotope['atomic number'])
            self.atomic_masses.append(isotope['atomic mass'])
            self.masses.append(isotope['mass'])
            self.abundances.append(isotope['abundance'])

        for m, s in zip(self.masses, self.counts):
            self.mass += m * s

        # adjust mass for extra or missing electrons (charge)
        if self.chargesign == '+':
            self.mass -= mass_electron * self.charge
        elif self.chargesign == '-':
            self.mass += mass_electron * self.charge

    def relative_abundance(self):
        """ Given a list of isotopes and a list of count numbers
            (subscripts) calculate relative abundance for entire molecule.
        """
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
        #   xi = number of isotope i = count
        #   pi = probability of isotope i = natural abundance
        #
        # Example: molecule 12C 16O2 18O
        # C is independent of O
        # there are 3 O in the molecule, n = 3
        # there are 2 O isotopes in this molecule, k = 2
        # for 16O: xi = 2, for 18O: xi = 1
        # for 16O: pi = 0.9976 for 18O: pi = 0.002 (and 0.0004 for 17O)

        data = periodic_table[periodic_table['isotope'].isin(self.isotopes)].copy()
        data['count'] = self.counts

        parents = data['major isotope'].value_counts().to_dict()
        abun_per_el = []
        for el, k in parents.items():
            d = data[data['major isotope'] == el]
            n = d['count'].sum()

            if k == 1:
                # Simple case of single isotope, even if it occurs n times
                abun = d['abundance'].iat[0] ** n
            else:
                abun = factorial(n)/factorial(d['count']).prod() * (d['abundance'] ** d['count']).prod()

            abun_per_el.append(abun)
        self.abundance = prod(abun_per_el)

    def sort(self):
        """ Sort elements from heavy to light, but keep same elements together. """
        raise NotImplementedError('not yet')

    def formula(self, style='plain', HtoD=True, show_charge=True, template={}):
        """ Return the molecular formula as a string.

            The molecular formula can be formatted as html
            (style='html'), LaTeX (style='latex'), plain
            text isotope notation (style='isotope', default),
            plain text empirical formula (style='empirical'),
            or in a custom format (style='custom'), see below.

            1H and 2H will be converted to H and D; set
            HtoD=False to output as 1H and 2H instead.

            Charge and sign will be automatically added, unless
            show_charge is set to False.

            If style='custom', a custom template can be used to
            format the molecular formula. The template must be
            a dict containing 6 keys: isotope, element, count,
            charge, minorjoin, and majorjoin. The isotope,
            element, count, and charge keys should refer to a
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

        elem = self.elements
        amass = [str(u) for u in self.atomic_masses]
        count = [str(c) if c > 1 else '' for c in self.counts]

        if HtoD:
            for n, x in enumerate(zip(amass, elem)):
                am, el = x
                if el == 'H':
                    if am == '1':
                        amass[n] = ''
                    elif am == '2':
                        amass[n] = ''
                        elem[n] = 'D'

        if style == 'html':
            templ = html_template
        elif style == 'latex':
            templ = latex_template
        elif style == 'empirical':
            raise NotImplementedError('not yet')
        elif style == 'custom':
            if not template:
                raise ValueError('If you select style="custom", you must supply a custom template.')
            templ = template
        elif style in ('plain', 'isotope'):
            templ = isotope_template
        else:
            raise ValueError('style must be one of "html", "latex", "plain", "isotope", or "custom".')

        molecule = []
        for am, el, ct in zip(amass, elem, count):
            if am:
                am_str = templ['isotope'].format(am)
            else:
                am_str = ''
            el_str = templ['element'].format(el)
            if ct:
                ct_str = templ['count'].format(ct)
            else:
                ct_str = ''
            m = templ['minorjoin'].join((am_str, el_str, ct_str))
            molecule.append(m)
        if charge:
            molecule.append(templ['charge'].format(charge))
        return templ['majorjoin'].join(molecule)


if __name__ == '__main__':
    m = Molecule('12C2 15N H4 2+')
    print(m.formula(style='html'))
    print(m.mass, m.abundance)
