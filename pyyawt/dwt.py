# -*- coding: utf-8 -*-

# Copyright (c) 2015 Holger Nahrstaedt
# See COPYING for license details.

"""
Helper function for wavelet denoising
"""

from __future__ import division, print_function, absolute_import
import numpy as np
import sys as sys
from ._pyyawt import *

__all__ = ['dbwavf', 'coifwavf']


def dbwavf(wname):
    ret = _wavelet_parser(wname)
    if (ret[0] != 1):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_dbwavf_length(wname),dtype=np.float64)
    _dbwavf(wname,lowPass)
    return lowPass


def coifwavf(wname):
    ret = _wavelet_parser(wname)
    if (ret[0] != 2):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_coifwavf_length(wname),dtype=np.float64)
    _coifwavf(wname,lowPass)
    return lowPass
