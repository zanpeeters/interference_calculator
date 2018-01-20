.. image:: interference_calculator/icon.svg
    :width: 128px
    :height: 128px
    :align: left
    :alt: icon.svg

***********************
Interference calculator
***********************

Interference calculator calculates all molecules that can be formed (the interferences) from a combination of a list of atoms (sample composition), given a target mass and range. The calculation considers all isotopes of the sample atoms and build molecules up to a given size. The results are displayed in both a table and a mass spectrum.

The program can also display the standard ratios of the isotopes for any given element. The results include the natural abundance, standard ratio, inverse ratio, and the standard material in which the isotopic ratios where measured.

Installation
============

To use interference_calculator, first you need to have Python installed. Download Python `here <https://www.python.org>`_. Once Python is installed, open a terminal window (command window on Windows) and type:
.. codeblock:: shell

    $ pip install interference_calculator

on the command line to install interference calculator.

Running
=======

To start the program, simply run the ui.py script.
.. codeblock:: shell

    $ python interference_calculator/ui.py

Use in IPython
==============

Interference_calculator can also be used from an interactive interpreter or in another Python script. For example, to calculate the mass interference around iron (Fe, mass 56), given a sample that consists of Si, Ca, O, and H, use it like this. ::

    >>> import interference_calculator as ic
    >>> ic.interference(['Ca', 'O', 'H', 'Si'], 'Fe')
             molecule  charge  mass/charge  mass/charge diff          MRP  \
    0          O Ca -       1    55.958054          0.023118  2419.529988
    1           Si2 -       1    55.954402          0.019466  2873.520086
    2        18O3 D -       1    56.012129          0.077193   724.609657
    3       18O3 H2 -       1    56.013677          0.078741   710.361726
    4  17O 18O2 H D -       1    56.019926          0.084990   658.132608
    5       D4 48Ca -       1    56.009478          0.074542   750.376520
    6     O 18O2 D2 -       1    56.021986          0.087050   642.561143
    7   17O2 18O D2 -       1    56.026175          0.091239   613.057735
    8            Fe -       1    55.934936          0.000000          inf

        probability  target
    0  9.671034e-01   False
    1  8.506314e-01   False
    2  7.556442e-15   False
    3  3.999731e-06   False
    4  1.103164e-13   False
    5  1.100126e-18   False
    6  3.721482e-09   False
    7  1.103556e-16   False
    8  9.175400e-01    True

    >>> ic.standard_ratio(['Ca', 'O'])
       isotope       mass  abundance     ratio  inverse ratio      standard
    13     16O  15.994915   0.997621  1.000000       1.000000         VSMOW
    14     17O  16.999132   0.000379  0.000380    2632.244327         VSMOW
    15     18O  17.999160   0.002000  0.002005     498.710558         VSMOW
    41    40Ca  39.962591   0.969410  1.000000       1.000000  NIST SRM 915
    42    42Ca  41.958618   0.006470  0.006674     149.831530  NIST SRM 915
    43    43Ca  42.958766   0.001350  0.001393     718.081481  NIST SRM 915
    44    44Ca  43.955482   0.020860  0.021518      46.472196  NIST SRM 915
    45    46Ca  45.953690   0.000040  0.000041   24235.250000  NIST SRM 915
    46    48Ca  47.952523   0.001870  0.001929     518.401070  NIST SRM 915

There is also a class ``Molecule``, that can parse strings with molecular formulas. After parsing, it holds information about the molecule, such as mass, elements, isotopes, and relative abundance. ``Molecule.formula()`` can be used to typeset the molecular formula in various ways.::

    >>> m = ic.Molecule('C2 15N O3 2+')
    >>> m.elements
    ['N', 'C', 'O']

    >>> m.isotopes
    ['15N', '12C', '16O']

    >>> m.masses
    [15.000108899, 12.0, 15.99491462]

    >>> m.mass
    86.985949918818

    >>> m.abundances
    [0.003663, 0.988922, 0.9976206]

    >>> m.abundance
    1.317443808955884e-05

    >>> m.formula(style='latex')
    '$\\mathrm{{}^{15}{N}{C}_{2}{O}_{3}{}^{2-}}$'

See the docstrings for detailed help and options.
