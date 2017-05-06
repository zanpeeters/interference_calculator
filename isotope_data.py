#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Manage table of accurate isotope masses and abundances. """

from pandas.io.sql import read_sql
from pandas import DataFrame
from html.parser import HTMLParser
import sqlite3
import urllib.request
import datetime, dateutil
import sys, os, warnings

__all__ = ['IsotopeData', 'update', 'get_isotope_data',
           'relative_abundance',
           'isotope_data', 'isotopes', 'atomic_numbers', 'elements',
           'unit_masses', 'isotope_masses', 'abundances',
           'main_isotopes', 'mass_electron'
]

_default_filename = 'isotope_data.db'
_default_path = os.path.dirname(__file__)
_default_location = os.path.join(_default_path, _default_filename)
# _mass_url = 'http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all&ascii=ascii2&isotyp=some'
_mass_url = 'http://ciaaw.org/atomic-masses.htm'
_abundance_url = ' https://www.degruyter.com/table/j/pac.2016.88.issue-3/pac-2015-0503/pac-2015-0503.xml?id=j_pac-2015-0503_tab_001'

# CODATA 2014, http://physics.nist.gov/cgi-bin/cuu/Value?me
mass_electron = 0.0005485799090

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


class IsotopeData(object):
    """ Maintain a local database of isotope data. """

    def __init__(self, filename=_default_location, new=False, _debug=False):
        """ Retrieve or update isotope data, store locally in a SQLite database.

            Usage: data = IsotopeData(filename='isotope_data.db', new=False)

            On first run, if 'filename' is not found, a new database will be
            created, by default in the installation directory of this module.
            Set 'new' to True to delete the old file and start from scratch,
            e.g. when a database gets corrupt for some reason.
        """
        # With _debug set, reads local file 'test.html' instead of fetching from website.
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
        """ IsotopeData supports the with statement. """
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
        """ Update the database from the website.
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


def get_isotope_data():
    """ Read isotope data from database, returns pandas.DataFrame with isotope names as index. """
    with IsotopeData() as idb:
        if idb.needs_update():
            warnings.warn('Database is outdated, run isotope_data.update().')
        return read_sql('SELECT * FROM isotopes', idb.db, index_col='isotope')

def update(period=365, force=False):
    """ Update database from website.
    
        Checks last update and updates database if data in database is older
        than period days (365 days by default). Force update now by setting
        force=True.
    """
    with IsotopeData() as idb:
        if idb.needs_update(period=period) or force:
            idb.update()
        else:
            print('Database is already up to date.')

# Load isotope_data only once, into namespace of module
# Split into separate python lists for faster lookup.
isotope_data = get_isotope_data()

isotopes = tuple(isotope_data.index.tolist())
atomic_numbers = tuple(isotope_data.get('atomic_number').tolist())
elements = tuple(isotope_data.get('element').tolist())
unit_masses = tuple(isotope_data.get('unit_mass').tolist())
isotope_masses = tuple(isotope_data.get('isotope_mass').tolist())
abundances = tuple(isotope_data.get('abundance').tolist())
main_isotopes = tuple(isotope_data.get('main_isotope').tolist())
