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

__all__ = ['dwt', 'idwt', 'wavedec', 'waverec', 'wrcoef', 'appcoef', 'detcoef', 'wenergy',
           'upcoef', 'upwlev']


def dwt(x,*args):
    """
    dwt is for discrete fast wavelet transform with the signal extension method optional argument. Available wavelets include haar, daubechies (db1 to db20), coiflets (coif1 to coif5), symlets (sym2 to sym20), legendre (leg1 to leg9), bathlets, dmey, beyklin, vaidyanathan, biorthogonal B-spline wavelets (bior1.1 to bior6.8).

    Calling Sequence
    ----------------

    [cA,cD]=dwt(x,wname,['mode',extMethod])
    [cA,cD]=dwt(x,Lo_D,Hi_D,['mode',extMethod])

    Parameters
    ----------
    wname: str
         wavelet name,  wavelet name, haar( "haar"), daubechies ("db1" to "db20"), coiflets ("coif1" to "coif5"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    x : array_like
         input vector
    Lo_D : array_like
         lowpass analysis filter
    Hi_D : array_like
         highpass analysis filter
    extMethod : str
         extension mode, 'zpd' for example
    Returns
    -------
    cA: array_like
         approximation coefficent
    cD: array_like
         detail coefficent
    Examples
    --------
    cA,cD=dwt(x,'db2','mode','asymh')
    """
    x = x.flatten()
    if (len(args) == 1 or (len(args) == 3 and args[1] == 'mode')):
        # flow 1
        wname = args[0]
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
        if (len(args) == 3 and args[1] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[2])
        _dwt_neo(x,Lo_D,Hi_D,out1,out2)
        if (len(args) == 3 and args[1] == 'mode'):
            dwtmode(ST)
        return out1, out2
    elif (len(args) == 2 or (len(args) == 3 and args[2] == 'mode')):
        # flow 1
        Lo_D = args[0]
        Hi_D = args[1]
        filterLength = args[0].shape[0]
        stride, val = _wave_len_validate(x.shape[0],filterLength)
        if (val == 0):
            raise Exception("Input signal is not valid for selected decompostion level and wavelets!")

        m3 = 1
        m4 = 1
        n3 = np.floor((x.shape[0] + filterLength - 1)/2).astype(int)
        if (dwtmode("status","nodisp") == 'per'):
            n3 = np.ceil(((x.shape[0]))/2.0).astype(int)
        n4 = n3
        out1 = np.zeros(n3*m3,dtype=np.float64)
        out2 = np.zeros(n4*m4,dtype=np.float64)
        if (len(args) == 3 and args[1] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[2])
        _dwt_neo(x,Lo_D,Hi_D,out1,out2)
        if (len(args) == 3 and args[1] == 'mode'):
            dwtmode(ST)
        return out1, out2
    else:
        raise Exception("Wrong input!")


def idwt(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wavedec(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def waverec(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wrcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def appcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def detcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def wenergy(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upwlev(cA, cD, *args):
    raise Exception("Not yet implemented!!")
