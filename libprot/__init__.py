# -*- coding: utf-8 -*-
import os
from pathlib import Path

"""Top-level package for libprot."""
from prody import confProDy

__author__ = """Nate Guerin"""
__email__ = 'nate@textpla.in'
__version__ = '0.1.2'

confProDy(verbosity='none')

REDUCE_BIN = Path(os.path.dirname(__file__), 'bin', 'reduce', 'reduce')
REDUCE_DB = Path(os.path.dirname(__file__), 'bin', 'reduce', 'reduce_wwPDB_het_dict.txt')
