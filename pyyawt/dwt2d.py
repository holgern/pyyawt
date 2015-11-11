# -*- coding: utf-8 -*-

# Copyright (c) 2015 Holger Nahrstaedt
# See COPYING for license details.

"""
dwt2d function for wavelet denoising
"""

from __future__ import division, print_function, absolute_import
import numpy as np
import sys as sys
from ._pyyawt import *
from .dwt import *

__all__ = ['dwt2', 'idwt2', 'wavedec2', 'waverec2', 'wenergy2', 'detcoef2', 'appcoef2', 'wrcoef2'
           'upcoef2', 'upwlev2']


def dwt2(x,*args):
    raise Exception("Not yet implemented!!")


def idwt2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wavedec2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def waverec2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wrcoef2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def appcoef2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def detcoef2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wenergy2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upcoef2(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upwlev2(cA, cD, *args):
    raise Exception("Not yet implemented!!")
