#!/usr/bin/env python
# -*- encoding: utf-8 -*-
""" Load the masstable module.
    
    Masstable uses a local database to store isotope data fetched from
    NIST and other wedsites. This file is stored in the same directory
    as where this module is installed. The installation path may require
    admin rights for write access. If you wish to update the database, or
    create one on first installation of this module, run the following
    as root/Administrator/God:
    
    on first install
    $ sudo python -M masstable -c "masstable.MassTable()"
    
    to update
    $ sudo python -M masstable -c "mt = masstable.MassTable(); mt.update()"
"""
from __future__ import absolute_import
from .masstable import *
