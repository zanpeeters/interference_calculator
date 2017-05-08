#!/usr/bin/env python
# -*- encoding: utf-8 -*-
""" Load the interference_calculator module.

    IsotopeData uses a local database to store isotope data fetched from
    NIST and other websites. This file is stored in the same directory
    as where this module is installed. The installation path may require
    admin rights for write access. If you wish to update the database, or
    create one on first installation of this module, run the following
    as root/Administrator/God:

    $ sudo python -M interference_calculator.isotope_data -c \
         "interference_calculator.isotope_data.update()"
"""
from __future__ import absolute_import
from .molecule import *
from .interference_calculator import *
