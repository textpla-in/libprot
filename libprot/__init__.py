# -*- coding: utf-8 -*-
import os
from pathlib import Path
from platform import system

"""Top-level package for libprot."""
from prody import confProDy

__author__ = """Nate Guerin"""
__email__ = 'nate@textpla.in'
__version__ = '0.1.4'

confProDy(verbosity='none')

REDUCE_BIN = Path(os.path.dirname(__file__), 'bin', 'reduce', system(), 'reduce')
REDUCE_DB = Path(os.path.dirname(__file__), 'bin', 'reduce', 'reduce_wwPDB_het_dict.txt')

if not REDUCE_BIN.exists():
    raise Exception(f'Platform {system()} is not supported. `reduce` was not found at {REDUCE_BIN.absolute()}')
REDUCE_BIN = str(REDUCE_BIN.absolute())
REDUCE_DB = str(REDUCE_DB.absolute())
