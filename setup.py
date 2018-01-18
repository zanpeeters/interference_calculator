#!/usr/bin/env python
""" Setuptools setup file for interference calculator. """
from setuptools import setup, find_packages
import os

with open(os.path.join('interference_calculator', '__init__.py'), mode='rt', encoding='utf-8') as fh:
    script = []
    for l in fh.readlines():
        l = l.strip()
        # don't import anything, just get metadata
        if l.startswith('from') or l.startswith('import'):
            continue
        script.append(l)

exec('\n'.join(script))

with open('README.md', mode='rt', encoding='utf-8') as fh:
    __long_description__ = fh.read()

try:
    import PyQt5
    pyqtdep = 'PyQt5'
except ImportError:
    try:
        import PyQt4
        pyqtdep = 'PyQt4'
    except ImportError:
        pyqtdep = 'PyQt5'

setup(
    name = __name__,
    version = __version__,
    description = __description__,
    long_description = __long_description__,
    url = __url__,
    author = __author__,
    author_email = 'me@example.com',
    license = __license__,
    classifiers = [
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering'
    ],
    keywords = 'interference mass spectrometry isotope element standard ratio',

    install_requires = [
        'matplotlib',
        'pandas',
        'pyparsing',
        pyqtdep,
        'scipy'
    ],
    entry_points = {
        'gui_scripts': ['interference_calculator=interference_calculator.ui:run']
    },

    packages = find_packages(),
    package_data = {'interference_calculator': [
        'periodic_table.csv',
        'icon.svg',
        'display_button_icon.svg',
        'help_button_icon.svg']
    },
    zip_safe = False
)
