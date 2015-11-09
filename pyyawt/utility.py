# -*- coding: utf-8 -*-

# Copyright (c) 2015 Holger Nahrstaedt
# See COPYING for license details.

"""
Helper function for pyyawt
"""

from __future__ import division, print_function, absolute_import
import numpy as np
import sys as sys
from ._pyyawt import *

__all__ = ['conv']


def conv(a,b):
    m1,n1 = a.shape
    m2,n2 = b.shape
    m3 = 1
    n3 = m1 * n1 + m2 * n2 - 1
    Y = np.zeros((m3,n3),dtype=np.float64)
    _conv(a,b,Y)
    return Y
