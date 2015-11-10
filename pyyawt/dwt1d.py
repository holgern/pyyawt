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
from .dwt import *

__all__ = ['dwt']


def dwt(x,wname):
    # flow 1
    ret = _wavelet_parser(wname.encode())
    filterLength = _wfilters_length(wname.encode())
    stride, val = _wave_len_validate(x.shape[0],filterLength)
    if (val == 0):
        raise Exception("Input signal is not valid for selected decompostion level and wavelets!")

    m3 = 1
    m4 = 1
    n3 = np.floor((x.shape[0] + filterLength - 1)/2).astype(int)
    if (dwtmode("status","nodisp") == 'per'):
        n3 = np.ceil(((x.shape[0]))/2.0).astype(int)
    n4 = n3
    Lo_D, Hi_D = wfilters(wname,'d')
    out1 = np.zeros(n3*m3,dtype=np.float64)
    out2 = np.zeros(n4*m4,dtype=np.float64)

    _dwt_neo(x,Lo_D,Hi_D,out1,out2)
    return out1, out2
