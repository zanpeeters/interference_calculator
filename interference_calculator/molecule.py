# -*- coding: utf-8 -*-
""" Molecule() is a class that takes an input string of a chemical formula,
    parses the string into atomic units, and stores relevant molecular data.
    The chemical formula can be output in a number of ways, including custom
    formatting using simple templates.
"""
import pandas as pd
import pyparsing as pp
import pkg_resources

from numpy import prod
from scipy.misc import factorial

_periodic_table_file = pkg_resources.resource_filename(__name__, 'periodic_table.csv')
periodic_table = pd.read_csv(_periodic_table_file, comment='#')

# CODATA 2014, http://physics.nist.gov/cgi-bin/cuu/Value?me
mass_electron = 0.0005485799090

# parser elements used by all forms
_opt_int = pp.Optional(pp.Word(pp.nums))
_element = pp.Combine(pp.Word(pp.alphas.upper(), exact=1) + pp.Optional(pp.Word(pp.alphas.lower(), max=2)))
_neutral = pp.oneOf('o 0')
_charged = pp.oneOf('+ -')

### isotope notation in Backus-Naur form (-ish)
# example: 12C2 18O -
#
# element     ::= capital + [lowercase letter] + [lowercase letter]
# atomic mass ::= integer
# count       ::= integer
# delimiter   ::= one or more character not A-Z, a-z, 0-9, +, -
# unit        ::= [atomic mass] + element + [multiplier] + [delimiter]
# charge      ::= ("o"|"0") | ([integer] + ("+"|"-"))
# molecule    ::= one or more units + [charge]
#

_in_delimiter = pp.CharsNotIn(pp.alphanums + '+-').setParseAction(pp.replaceWith(','))
_in_comma = pp.Optional(pp.Suppress(','))
_in_unit = pp.OneOrMore(pp.Group(
                _opt_int('atomic_mass') + _element('element') + _opt_int('count') + _in_comma
              ))
_in_charge = pp.Optional(
                _neutral('charge_sign') |
                _opt_int('charge_count') + _charged('charge_sign')
             )
_in_molecule = _in_unit('units') + _in_charge

### molecular notation in Backus-Naur form (-ish)
# example: C2H5COOCH[15]NH3[+]
#
# element     ::= capital + [lowercase letter] + [lowercase letter]
# atomic mass ::= "[" + integer + "]"
# count       ::= integer
# charge      ::= "[" + (("o"|"0") | ([integer] + ("+"|"-"))) + "]"
# unit        ::= [atomic mass] + element + [multiplier]
# molecule    ::= one or more units + [charge]

_mn_atomic_mass = pp.Optional(pp.Combine(
                    pp.Suppress('[') + pp.Word(pp.nums) + pp.Suppress(']')
                  ))
_mn_unit = pp.OneOrMore(pp.Group(
                _mn_atomic_mass('atomic_mass') + _element('element') + _opt_int('count')
              ))
_mn_charge = pp.Optional(
                pp.Suppress('[') + (
                    _neutral('charge_sign') |
                    _opt_int('charge_count') + _charged('charge_sign')
                ) + pp.Suppress(']')
            )
_mn_molecule = _mn_unit('units') + _mn_charge

# Just a list of names to choose from, for convenience.
templates = ['html', 'latex', 'mhchem', 'isotope', 'plain', 'molecular']

html_template = {
    'begin': '',
    'atomic_mass': '<sup>{}</sup>',
    'element': '{}',
    'count': '<sub>{}</sub>',
    'charge': '<sup>{}</sup>',
    'minorjoin': '',
    'majorjoin': '',
    'minus': '&ndash;',
    'end': ''
}

latex_template = {
    'begin': '$\mathrm{',
    'atomic_mass': '{{}}^{{{}}}',
    'element': '{{{}}}',
    'count': '_{{{}}}',
    'charge': '{{}}^{{{}}}',
    'minorjoin': '',
    'majorjoin': '',
    'minus': '-',
    'end': '}$'
}

mhchem_template = {
    'begin': '\ce{',
    'atomic_mass': '^{{{}}}',
    'element': '{}',
    'count': '{}',
    'charge': '^{}',
    'minorjoin': '',
    'majorjoin': '',
    'minus': '',
    'end': '}'
}

isotope_template = {
    'begin': '',
    'atomic_mass': '{}',
    'element': '{}',
    'count': '{}',
    'charge': '{}',
    'minorjoin': '',
    'majorjoin': ' ',
    'minus': '',
    'end': ''
}

molecular_template = {
    'begin': '',
    'atomic_mass': '[{}]',
    'element': '{}',
    'count': '{}',
    'charge': '[{}]',
    'minorjoin': '',
    'majorjoin': '',
    'minus': '',
    'end': ''
}

class Molecule(object):
    """ Represents a molecule or molecular ion. """

    def __init__(self, molecule):
        """ Parses a chemical formula string and returns an object that
            holds properties of the molecule or molecular ion.

            Two forms of input string are supported: isotope notation and
            molecular formula notation. These names and notations are used
            for input and output.

            Isotope notation is a list of units, where each unit is of
            the form NXxxn, where N the atomic mass, Xxx is the element,
            and n is the count (subscript). Any character except A-Z, 0-9,
            +, -, [, or ] may be used to separate the units in the list,
            space is most common. If no atomic mass is specified, the most
            common isotope is assumed (e.g. C -> 12C). A charge may optionally
            be given as the last element. This notation is useful for inputting
            many unusual isotopes.

            Isotope notation: '12C2 15N O3 2-'

            Molecular formula notation is a form of shorthand. It contains
            no spaces and no atomic masses, only count numbers. If an atomic
            mass needs to be given for an isotope, it must be surrounded by [].
            A charge may optionally be given at the end, also surrounded by [].
            This form is useful for inputting larger molecules with few
            unusual isotopes.

            Molecular formula notation: 'HCOOCH2[15]NH3[2-]'

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
        self.isotopes = []
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
        if not self.input:
            return
        self.input = self.input.strip()

        # Parse input string into pyparsing.ParseResult objects
        try:
            molec = _mn_molecule.parseString(self.input, parseAll=True)
        except pp.ParseException:
            delim_string = _in_delimiter.transformString(self.input)
            molec = _in_molecule.parseString(delim_string, parseAll=True)

        # Collect data from ParseResult objects,
        # merge mulitple occurances of same element.
        data = {}
        for unit in molec.units:
            label = unit.atomic_mass + unit.element
            if label not in data.keys():
                data[label] = {
                    'atomic_mass': unit.atomic_mass,
                    'element': unit.element,
                    'count': int(unit.get('count', 1))
                }
            else:
                data[label]['count'] += int(unit.get('count', 1))

        # Sort and split data into lists.
        for k in sorted(data.keys()):
            am = data[k]['atomic_mass']
            el = data[k]['element']
            if el == 'D':
                # special case
                el = 'H'
                am = 2
            elif am:
                am = int(am)
            else:
                # no atomic mass given, find major isotope, e.g. C -> 12C
                am = periodic_table[periodic_table['element'] == el].iloc[0].loc['major isotope']
                am = int(am.strip(el))
            self.atomic_masses.append(am)
            self.elements.append(el)
            self.isotopes.append(str(am) + el)
            self.counts.append(data[k]['count'])

        # Retrieve additional information from periodic table
        for i in self.isotopes:
            isotope = periodic_table[periodic_table['isotope'] == i].iloc[0]
            self.atomic_numbers.append(isotope['atomic number'])
            self.masses.append(isotope['mass'])
            self.abundances.append(isotope['abundance'])

        # Calculate total mass of molecule
        for m, c in zip(self.masses, self.counts):
            self.mass += m * c

        # Find charge and sign
        self.chargesign = molec.get('charge_sign', '')
        if self.chargesign in ('o', '0', ''):
            self.charge = 0
        else:
            self.charge = int(molec.get('charge_count', 1))

        # Adjust mass for extra or missing electrons (charge)
        if self.chargesign == '+':
            self.mass -= mass_electron * self.charge
        elif self.chargesign == '-':
            self.mass += mass_electron * self.charge

    def relative_abundance(self):
        """ Given a list of isotopes and a list of count numbers
            (subscripts) calculate relative abundance for entire molecule.
        """
        # multiple isotopes e.g. 28Si (92.2%) 29Si (4.7%) 30Si (3.1%)
        # In this type of mass spectrometry we only look at total mass of molecule,
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

    def formula(self, style='plain', HtoD=True, show_charge=True, all_isotopes=False, template={}):
        """ Return the molecular formula as a string.

            The molecular formula can be formatted as html
            (style='html'), LaTeX (style='latex'), LaTeX with
            mhchem package (style='mhchem'), plain
            text isotope notation (style='isotope' or
            'plain', default), molecular formula notation
            (style='molecular'), or in a custom format
            (style='custom'), see below.

            1H and 2H will be converted to H and D; set
            HtoD=False to output as 1H and 2H instead.

            Charge and sign will be automatically added, unless
            show_charge is set to False.

            Only atomic masses of minor isotopes will be given
            in the output; set all_isotopes=True to explicitly
            give the atomic mass for each isotope. HtoD takes
            precedence over all_isotopes: H and D will be used
            even if all_isotopes=True, unless HtoD=False.

            If style='custom', a custom template can be used to
            format the molecular formula. The template must be
            a dict containing 9 keys: begin, atomic_mass, element,
            count, charge, minorjoin, majorjoin, minus, and end.
            A curly brace pair in the atomic_mass, element, count,
            and charge template strings will be replaced using
            string.format(). The minorjoin string will be used
            to join the atomic mass, element, and count
            into a unit, and the majorjoin string connects all
            the units into the output string. The begin and end
            strings are added to the beginning and end of the final
            output string, respectively.
        """
        # Force copy of list without using list.copy (python 2)
        elem = list(self.elements)
        amass = [str(u) for u in self.atomic_masses]
        count = [str(c) if c > 1 else '' for c in self.counts]

        if HtoD:
            for n, (am, el) in enumerate(zip(amass, elem)):
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
        elif style == 'mhchem':
            templ = mhchem_template
        elif style == 'molecular':
            templ = molecular_template
        elif style in ('plain', 'isotope'):
            templ = isotope_template
        elif style == 'custom':
            if not template:
                raise ValueError('If you select style="custom", you must supply a custom template.')
            templ = template
        else:
            msg = 'style must be one of "html", "latex", "mhchem", '
            msg += '"plain", "isotope", "molecular", or "custom".'
            raise ValueError(msg)

        if show_charge:
            if self.chargesign == '-' and templ['minus']:
                chargesign = templ['minus']
            else:
                chargesign = self.chargesign

            if self.charge == 0:
                charge = ''
            elif self.charge == 1:
                charge = chargesign
            else:
                charge = str(self.charge) + chargesign
        else:
            charge = ''

        molecule = []
        for am, el, ct in zip(amass, elem, count):
            if am:
                if (not all_isotopes and
                    (periodic_table['major isotope'] == am + el).any()):
                        am_str = ''
                else:
                    am_str = templ['atomic_mass'].format(am)
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

        return templ['begin'] + templ['majorjoin'].join(molecule) + templ['end']
