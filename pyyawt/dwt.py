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
    """
    dbwavf is an utility function for obtaining scaling filter of daubechies wavelet.

    Parameters
    ----------
    wname: str
         wavelet name, 'db1' to 'db36'

    Returns
    -------
    F: array_like
         scaling filter

    Examples
    --------
    F = dbwavf("db2")
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != 1):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_dbwavf_length(wname.encode()),dtype=np.float64)
    _dbwavf(wname.encode(),lowPass)
    return lowPass


def coifwavf(wname):
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != 2):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_coifwavf_length(wname.encode()),dtype=np.float64)
    _coifwavf(wname.encode(),lowPass)
    return lowPass
