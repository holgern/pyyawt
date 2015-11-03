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

__all__ = ['dbwavf', 'coifwavf', 'symwavf', 'legdwavf']


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
    if (ret[0] != PYYAWT_DAUBECHIES):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_dbwavf_length(wname.encode()),dtype=np.float64)
    _dbwavf(wname.encode(),lowPass)
    return lowPass


def coifwavf(wname):
    """
    coifwavf is an utility function for obtaining scaling filter of coiflets wavelet.

    Parameters
    ----------
    wname: str
         wavelet name, 'coif1' to 'coif5'

    Returns
    -------
    F: array_like
         scaling filter

    Examples
    --------
    F = coifwavf('coif3')
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != PYYAWT_COIFLETS):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_coifwavf_length(wname.encode()),dtype=np.float64)
    _coifwavf(wname.encode(),lowPass)
    return lowPass


def symwavf(wname):
    """
    symwavf is an utility function for obtaining scaling filter of symlets wavelet.

    Parameters
    ----------
    wname: str
         wavelet name, 'sym2' to 'sym20'

    Returns
    -------
    F: array_like
         scaling filter

    Examples
    --------
    F = symwavf('sym7')
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != PYYAWT_SYMLETS):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_symwavf_length(wname.encode()),dtype=np.float64)
    _symwavf(wname.encode(),lowPass)
    return lowPass


def legdwavf(wname):
    """
    legdwavf is an utility function for obtaining scaling filter of legendre wavelet.

    Parameters
    ----------
    wname: str
         wavelet name, 'legd1' to 'legd9'

    Returns
    -------
    F: array_like
         scaling filter

    Examples
    --------
    F = legdwavf('sym7')
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != PYYAWT_LEGENDRE):
        raise Exception("Wrong wavelet name!")
    lowPass = np.zeros(_legdwavf_length(wname.encode()),dtype=np.float64)
    _legdwavf(wname.encode(),lowPass)
    return lowPass
