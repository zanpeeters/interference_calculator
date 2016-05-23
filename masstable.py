#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Manage table of accurate isotope masses and abundances, downloaded from NIST. """

from pandas.io.sql import read_sql
from pandas import Series, DataFrame, concat
from numpy import nan, prod
from scipy.misc import factorial
from html.parser import HTMLParser
from collections import Counter
import sqlite3
import urllib.request
import datetime, dateutil
import sys, os, itertools, re

__all__ = ['Molecule', 'MassTable',  # classes
           'mass_interference', 'get_masstable',  # methods
           'standard_ratios', 'relative_abundance',
           'masstable',  # main table
           'isotopes', 'atomic_numbers', 'elements',  # lists
           'unit_masses', 'isotope_masses', 'abundances',
           'main_isotopes'
]

_default_filename = 'masstable.db'
_default_path = os.path.dirname(__file__)
_default_location = os.path.join(_default_path, _default_filename)
_mass_url = 'http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all&ascii=ascii2&isotyp=some'
# _url = 'http://ciaaw.org/atomic-masses.htm'

# Precompile regexes
# 
# How to determine the notation?
# If it contains a separation marker (space or anything not letter, numbers, or []+-),
# then it is isotope notation, otherwise empirical.
# However, treat a single element (one capital in string) as isotope notation.
_isotope_notation_rx = re.compile(r'(\d*)([A-Z][a-z]{0,2}|[+-])(\d*)')
_empirical_notation_rx = re.compile(r'(?:\[(\d*)\])?([A-Z][a-z]{0,2}|\[\d*[+-]\])(\d*)')
_is_empirical_rx = re.compile(r'^[A-Za-z\d\[\]+-]+$')
_is_single_rx = re.compile(r'[A-Z]')

# Molecular formula formatting templates
_html_template = {
    'isotope': '<sup>{}</sup>',
    'element': '{}',
    'stoich': '<sub>{}</sub>',
    'charge': '<sup>{}</sup>',
    'minorjoin': '',
    'majorjoin': ''
}

_latex_template = {
    'isotope': '{{}}^{{{}}}',
    'element': '{{{}}}',
    'stoich': '_{{{}}}',
    'charge': '{{}}^{{{}}}',
    'minorjoin': '',
    'majorjoin': ''
}

_isotope_template = {
    'isotope': '{}',
    'element': '{}',
    'stoich': '{}',
    'charge': '{}',
    'minorjoin': '',
    'majorjoin': ' '
}


class HTMLPreParser(HTMLParser):
    """ Extracts <pre>...</pre> blocks from a HTML page.
        All pre-blocks are concatenated into a single text string.
    """
    def __init__(self):
        super().__init__()
        self._in_pre = False
        self.text = ''

    def handle_starttag(self, tag, attrs):
        """ Handle html opening tags. """
        if tag == 'pre':
            self._in_pre = True

    def handle_data(self, data):
        """ Handle text inside html tags. """
        if self._in_pre:
            self.text += data

    def handle_endtag(self, tag):
        """  Handle html closing tag. """
        if tag == 'pre':
            self._in_pre = False


class MassTable(object):
    """ MassTable class, maintains a local database of isotope data. """
    # Support 2 forms of object instantiation:
    # 1. a = A()
    # 2. with A() as a:
    #
    # In case 1, __init__ is called on instantiation, __del__ on destruction,
    # even if there was an error. Also when e.g. a = None.
    # In case 2:
    #   instantiation: __init__ then __enter__
    #   destruction: __exit__ then __del__
    # __exit__ is send an error object, or None.
    #
    # To support both cases, let __init__ and __del__ do all the work.
    # __enter__ only needs to return the self object.
    #
    def __init__(self, filename=_default_location, new=False, _debug=False):
        """ Retrieve or update NIST isotope data, store locally in a SQLite database.

            Usage: mt = MassTable(filename='masstable.db', new=False)

            On first run, if 'filename' is not found, a new database will be
            created, by default in the installation directory of this module.
            Set 'new' to True to delete the old file and start from scratch,
            e.g. when a database gets corrupt for some reason.
        """
        # With _debug set, reads local file 'test.html' instead of fetching from NIST website.
        self._debug = _debug
        self.statements = {
            'isotopes_exists': '''SELECT name FROM sqlite_master
                                 WHERE type="table" AND name="isotopes";''',
            'isotopes_create': '''CREATE TABLE isotopes
                            ({} int not null,
                             {} text not null,
                             {} int not null,
                             {} text unique not null,
                             {} real not null,
                             {} real);''',
            'isotopes_fill': '''INSERT INTO isotopes
                          ({}, {}, {}, {}, {}, {})
                          VALUES (?,?,?,?,?,?);''',
            'isotopes_drop': '''DROP TABLE IF EXISTS isotopes;''',
            'metadata_create': '''CREATE TABLE metadata
                                 (var text unique not null,
                                  value text);''',
            'metadata_fill':  '''INSERT INTO metadata
                               (var, value)
                               VALUES (?,?);''',
            'metadata_drop': '''DROP TABLE IF EXISTS metadata;''',
            'metadata_getvar': '''SELECT value FROM metadata WHERE var=?;''',
            'abundance_create': '''CREATE TABLE abundance
                                   ({} text unique not null,
                                    {} real,
                                    {} text);''',
            'abundance_fill': '''INSERT INTO abundance
                                ({}, {}, {}) VALUES (?,?,?);''',
            'abundance_drop': '''DROP TABLE IF EXIST abundance;'''
        }
        self.metadata = {
            'last_update': None,
            'mass_url': _mass_url,
            'abundance_url': _abundance_url
        }

        if new:
            try:
                os.unlink(filename)
            except FileNotFoundError:
                pass

        # sqlite3.connect will create empty file if it doesn't exist.
        db_exists = os.path.exists(filename)

        self.db = sqlite3.connect(filename)
        self.cursor = self.db.cursor()

        if not db_exists:
            self.update()

    def __enter__(self):
        """ MassTable supports the with statement. """
        return self

    def __del__(self):
        """ Destructor for class, make sure SQLite connection is closed.
            Rollback if there was an error, commit otherwise.
        """
        if self.db:
            if sys.exc_info()[1]:
                self.db.rollback()
            else:
                self.db.commit()
            self.db.close()

    def __exit__(self, errortype, error, traceback):
        """ Finalizes the with-statement. Defer to __del__. """
        pass

    def needs_update(self, period=365):
        """ Test whether data in database is older than period days.
            Default period is 1 year (365 days). Returns True or False.
        """
        query = self.cursor.execute(self.statements['metadata_getvar'], ('last_update',))
        lastupdate = query.fetchone()[0]
        lastupdate = dateutil.parser.parse(lastupdate)
        period = datetime.timedelta(days=period)
        now = datetime.datetime.now()
        if now > lastupdate + period:
            return True
        else:
            return False

    def update(self):
        """ Update the database from the NIST website.
            Change the download URL in self.metadata['url'] before calling update().
        """
        tbl = HTMLPreParser()
        if self._debug:
            f = open('test.html', encoding='ascii', mode='rt')
            tbl.feed(f.read())
        else:
            with urllib.request.urlopen(self.metadata['url']) as response:
                info = response.info()
                charset = info.get_content_charset()
                if not charset:
                    charset = 'utf-8'
                tbl.feed(response.read().decode(charset))

        # Data is partitioned in blocks separated by empty lines
        # First and last block are links.
        blocks = tbl.text.split('\n\n')[1:-1]

        table = []
        isotopes = []
        for block in blocks:
            data = []
            for line in block.split('\n'):
                data.append(line.split('=')[1].strip())

            atomic_number = int(data[0])
            element = data[1]
            unit_mass = int(data[2])
            isotope = '{}{}'.format(unit_mass, element)
            isotope_mass = float(data[3].split('(')[0])
            try:
                abundance = float(data[4].split('(')[0])
            except ValueError:
                abundance = None

            if element in ('H', 'D', 'T'):
                if element == 'H':
                    alias = '1H'
                    isotope = 'H'
                else:
                    alias = '{}H'.format(unit_mass)
                    isotope = element
                    element = 'H'
                table.append([atomic_number, element, unit_mass,
                              isotope_mass, abundance, ''])
                isotopes.append(alias)

            table.append([atomic_number, element, unit_mass,
                          isotope_mass, abundance, ''])
            isotopes.append(isotope)

        headers = ['atomic_number', 'element', 'unit_mass',
                   'isotope_mass', 'abundance', 'main_isotope']
        d = DataFrame(data=table, columns=headers, index=isotopes)

        # Find the main isotope for each element.
        # Determined by (1) maximum abundance, or (2) lowest unit mass
        # if no abundance value is available.
        for element in d.element.unique():
            sel = d.element == element
            if all(d.loc[sel, 'abundance'].isnull()):
                main_isotope = d.loc[sel, 'unit_mass'].argmin()
            else:
                main_isotope = d.loc[sel, 'abundance'].argmax()
            d.loc[sel, 'main_isotope'] = main_isotope

        d.to_sql('isotopes', self.db, if_exists='replace', index_label='isotope')

        self.metadata['last_update'] = datetime.datetime.now()
        self.cursor.execute(self.statements['metadata_drop'])
        self.cursor.execute(self.statements['metadata_create'])
        self.cursor.executemany(self.statements['metadata_fill'], self.metadata.items())
        self.db.commit()


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
        molecule = molecule.strip()

        # Determine notation
        if (len(re.findall(_is_single_rx, molecule)) == 1 or
            not re.match(_is_empirical_rx, molecule)):
                units = re.findall(_isotope_notation_rx, molecule)
        else:
            units = re.findall(_empirical_notation_rx, molecule)

        # Check for charge an sign in input
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
                idx = elements.index(u[1])
                mi = main_isotopes[idx]
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
        self.counter = c  # debug
        self.isotopes = list(c.keys())
        # TODO: make Molecule sortable, use self.sort()
        self.isotopes.sort()
        self.isotopes = tuple(self.isotopes)

        # Don't use c.values, isot/keys are sorted
        self.stoichiometry = tuple(c[i] for i in self.isotopes)
        self.indices = tuple(isotopes.index(i) for i in self.isotopes)
        self.elements = tuple(elements[i] for i in self.indices)
        self.atomic_numbers = tuple(atomic_numbers[i] for i in self.indices)
        self.unit_masses = tuple(unit_masses[i] for i in self.indices)
        self.isotope_masses = tuple(isotope_masses[i] for i in self.indices)
        self.abundances = tuple(abundances[i] for i in self.indices)

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
            templ = _html_template
        elif style == 'latex':
            templ = _latex_template
        elif style == 'empirical':
            raise NotImplementedError('not yet')
        elif style == 'custom':
            templ = template
        else:
            templ = _isotope_template

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


def get_masstable():
    """ Read mass table from database, returns pandas.DataFrame with isotope names as index. """
    with MassTable() as mt:
        if mt.needs_update():
            warnings.warn('Data in table is outdated, run MassTable.update().')
        return read_sql('SELECT * FROM isotopes', mt.db, index_col='isotope')


def mass_interference(atoms, mz, mzrange=0.3, maxsize=5, charge=[1], chargesign='-',
                      minabundance=1e-8, formula_style='html'):
    """ For a list of atoms, calculate all molecular ions that can be formed
        from those atoms, including all isotopes, up to maxsize atoms,
        that have a mass-to-charge ratio within mz +/- mzrange.

        The mass-to-charge ratio (mz) can be given as a number, or as a
        molecular formula. Molecular formulas are interpreted by Molecule().
        See Molecule() docstring for a detailed explanation of how to enter
        molecular formulas. If mz is None, no filtering will be done and
        all possible combinations of all isotopes up to maxsize length will
        be calculated.

        Charge is usually 1, irrespective of sign. Give charge = [1, 2, 3]
        to also include higher charged ions. Chargesign is only used for display.

        Not all isotopes have a known natural abundance. In the case where an
        abundance value is missing, use minabundance as a replacement. Setting
        this to None, will use NaN instead, and all ions with an isotope for
        which the abundance is not known, will not be calculated (i.e. NaN is
        respected).

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
        picked_atoms.append(masstable[masstable.element == a])

    picked_atoms = concat(picked_atoms).query('index not in ("H", "D", "T")')

    if not minabundance:
        # Signalling value, place back later
        minabundance = -1
    picked_atoms.loc[:, 'abundance'] = picked_atoms.loc[:, 'abundance'].fillna(value=minabundance)

    # Mass-to-charge can be given as either a number, or as a molecule (string).
    # Calculate m/z from molecular formula.
    if mz:
        try:
            mz = float(mz)
        except ValueError:
            molec = Molecule(mz)
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
    # print(abundances)
    abundances = abundances.fillna(value=1)
    abundances = abundances.where(abundances >= 0)
    abundances = abundances.prod(axis=1, skipna=False)

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
        m = Molecule(s)
        return m.formula(style=formula_style)

    # Next line gives SettingWithCopyWarning no matter how I do it.
    results.loc[:, 'molecule'] = results['molecule'].apply(convert)
    results = results.ix[:, ['molecule',
                             'charge',
                             'mass/charge',
                             'mass/charge diff',
                             'MRP',
                             'probability']]

    return results.sort(columns='mass/charge')


def standard_ratios(atoms, formula_style='html'):
    """ Give the stable isotopes and solar abundances for the given elements. """
    # Which cols: isotope, mass, std ratio, std ratio %
    data = []
    for a in atoms:
        eldata = masstable[masstable['element'] == a]
        highest = eldata['abundance'].max()
        eldata['inverse'] = highest/eldata['abundance']
        # eldata['percent'] = 100 * eldata['abundance']
        data.append(eldata)

    results = concat(data)

    def convert(s):
        m = Molecule(s)
        return m.formula(style=formula_style, show_charge=False)

    results.index = [convert(i) for i in results.index]
    results = results.ix[:, ['isotope_mass', 'abundance', 'inverse']] #, 'percent']]

    return results


def relative_abundance(isot, stoich):
    """ Given a list of isotopes and a list of stoichiometric numbers
        calculate relative abundance for entire molecule. """
    # multiple isotopes e.g. 28Si (92.2%) 29Si (4.7%) 30Si (3.1%)
    # In this type of mass spectorscopy we only look at mass,
    # not position of isotope. Therefore Si4-29Si has 5 equivalent structures:
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
    #   k = number of isotopes in molecule
    #   xi = number of isotope i, stoichiometry
    #   pi = probability of isotope i, natural abundance
    #
    # Example: molecule 12C 16O2 18O
    # C is independent of O
    # there are 3 O in the molecule, n = 3
    # there are 2 O isotopes in this molecule, k = 2
    # for 16O: xi = 2, for 18O: xi = 1
    # for 16O: pi = 0.9976 for 18O: pi = 0.002 (and 0.0004 for 17O)

    data = masstable.loc[isot]
    data['stoich'] = stoich

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


# Load masstable only once, into namespace of module
# Split into separate python lists for faster lookup.
masstable = get_masstable()

isotopes = tuple(masstable.index.tolist())
atomic_numbers = tuple(masstable.get('atomic_number').tolist())
elements = tuple(masstable.get('element').tolist())
unit_masses = tuple(masstable.get('unit_mass').tolist())
isotope_masses = tuple(masstable.get('isotope_mass').tolist())
abundances = tuple(masstable.get('abundance').tolist())
main_isotopes = tuple(masstable.get('main_isotope').tolist())

if __name__ == '__main__':
    print(standard_ratios(['C', 'S']))
