#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Create parser for elements, isotopes, and molecules based on pyparser. """

import pyparsing as pp

_opt_int = pp.Optional(pp.Word(pp.nums))
_element = pp.Combine(pp.Word(pp.alphas.upper(), exact=1) + pp.Optional(pp.Word(pp.alphas.lower(), max=2)))
_neutral = pp.oneOf('o 0')
_charged = pp.oneOf('+ -')

### isotope formula in Backus-Naur form (-ish)
# example: 12C2 18O -
#
# element     ::= capital + [lowercase letter] + [lowercase letter]
# atomic mass ::= integer
# multiplier  ::= integer
# delimiter   ::= one or more character not A-Z, a-z, 0-9, +, -
# isotope     ::= [atomic mass] + element + [multiplier] + [delimiter]
# charge      ::= ("o"|"0") | ([integer] + ("+"|"-"))
# molecule    ::= one or more isotope + [charge]
#

_if_delimiter = pp.CharsNotIn(pp.alphanums + '+-').setParseAction(pp.replaceWith(','))
_if_comma = pp.Optional(pp.Suppress(','))
_if_isotope = pp.OneOrMore(pp.Group(
                _opt_int('atomic_mass') + _element('element') + _opt_int('multiplier') + _if_comma
              ))
_if_charge = pp.Optional(pp.Group(
                _neutral('charge_sign') | 
                _opt_int('charge_count') + _charged('charge_sign')
             ))
_if_molecule = _if_isotope('isotopes') + _if_charge('charge')

def parse_isotope_formula(input_string):
    """ Parses a string in the 'isotope formula' format, returns a
        pyparsing.ParseResult object.

        For the returned molecule, molecule.isotopes gives a list
        of all elements with atomic mass and multipliers, while charge
        information is given in molecule.charge.

        The resulting lists are not uniq'ed, multiple occurences of
        the same elements/isotope are not combined.
    """
    clean_string = _if_delimiter.transformString(input_string)
    return _if_molecule.parseString(clean_string, parseAll=True)

### molecular formula in Backus-Naur form (-ish)
# example: C2H5COOCH[15]NH3[+]
#
# element     ::= capital + [lowercase letter] + [lowercase letter]
# atomic mass ::= "[" + integer + "]"
# multiplier  ::= integer
# charge      ::= "[" + (("o"|"0") | ([integer] + ("+"|"-"))) + "]"
# isotope     ::= [atomic mass] + element + [multiplier]
# molecule    ::= one or more isotope + [charge]

_mf_allowed = pp.CharsNotIn(pp.alphanums + '[]+-')
_mf_atomic_mass = pp.Optional(pp.Combine(
                    pp.Suppress('[') + pp.Word(pp.nums) + pp.Suppress(']')
                  ))
_mf_isotope = pp.OneOrMore(pp.Group(
                _mf_atomic_mass('atomic_mass') + _element('element') + _opt_int('multiplier')
              ))
_mf_charge = pp.Optional(
                pp.Suppress('[') + (
                    _neutral('charge_sign') | 
                    _opt_int('charge_count') + _charged('charge_sign')
                ) + pp.Suppress(']')
            )
_mf_molecule = _mf_isotope('isotopes') + _mf_charge('charge')

def parse_molecular_formula(input_string):
    """ Parses a string in the 'molecular formula' format, returns a
        pyparsing.ParseResult object.

        For the returned molecule, molecule.isotopes gives a list
        of all elements with atomic mass and multipliers, while charge
        information is given in molecule.charge.

        The resulting lists are not uniq'ed, multiple occurences of
        the same elements/isotope are not combined.
    """
    return _mf_molecule.parseString(input_string, parseAll=True)

def is_molecular_formula(input_string):
    """ Parses a string to determine if it is in the 'molecular format'.

        Returns True if the input string is in the 'molecular format',
        False otherwise.
    """
    return _mf_allowed.matches(input_string)

def parse(input_string):
    """ Parses an input string and retrieves information about the molecule
        it represents.
    """
    if is_molecular_formula(input_string):
        molec = parse_molecular_formula(input_string)
    else:
        molec = parse_isotope_formula(input_string)

    data = {}
    for m in molec.isotopes:
        label = m.atomic_mass + m.element
        if label not in data.keys():
            data[label] = {
                'atomic_mass': int(m.atomic_mass),
                'element': m.element,
                'multiplier': int(m.multiplier)
            }
        else:
            x = m.get('multiplier', 1)
            data[label]['multiplier'] += int(x)

    atomic_masses = []
    elements = []
    multipliers = []
    for k in sorted(data.keys()):
        atomic_masses.append(data[k]['atomic_mass'])
        elements.append(data[k]['element'])
        multipliers.append(data[k]['multiplier'])
    
    charge = molec.charge.charge_number
    charge_sign = molec.charge.charge_sign
    print(atomic_masses)
    print(elements)
    print(multipliers)
    print(charge, charge_sign)

