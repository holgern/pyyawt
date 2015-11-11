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

__all__ = ['conv', 'iconv', 'wrev', 'qmf', 'dyaddown', 'dyadup', 'wkeep', 'wextend', 'wcodemat',
           'mat3Dtran', 'wrev3', 'wrev2', 'wnorm', 'waveletfamilies']


def conv(a,b):
    m1,n1 = a.shape
    m2,n2 = b.shape
    m3 = 1
    n3 = m1 * n1 + m2 * n2 - 1
    Y = np.zeros((m3,n3),dtype=np.float64)
    _conv(a,b,Y)
    return Y


def iconv(*args):
    raise Exception("Not yet implemented!!")


def wrev(*args):
    raise Exception("Not yet implemented!!")


def qmf(*args):
    raise Exception("Not yet implemented!!")


def dyaddown(*args):
    raise Exception("Not yet implemented!!")


def dyadup(*args):
    raise Exception("Not yet implemented!!")


def wkeep(*args):
    raise Exception("Not yet implemented!!")


def wextend(*args):
    raise Exception("Not yet implemented!!")


def wcodemat(*args):
    raise Exception("Not yet implemented!!")


def mat3Dtran(*args):
    raise Exception("Not yet implemented!!")


def wrev3(*args):
    raise Exception("Not yet implemented!!")


def wrev2(*args):
    raise Exception("Not yet implemented!!")


def wnorm(*args):
    raise Exception("Not yet implemented!!")


def waveletfamilies(*args):
    raise Exception("Not yet implemented!!")
