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


def wavedec(x, N, *args):
    """
    Multiple level 1-D discrete fast wavelet decomposition
    Calling Sequence
    [C,L]=wavedec(X,N,wname)
    [C,L]=wavedec(X,N,Lo_D,Hi_D)
    Parameters
    wname : wavelet name, haar( "haar"), daubechies ("db1" to "db36"), coiflets ("coif1" to "coif17"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    X : signal vector
    N : decompostion level
    Lo_D : lowpass analysis filter
    Hi_D : highpass analysis filter
    C : coefficient vector
    L : length vector
    Description
    wavedec can be used for multiple-level 1-D discrete fast wavelet
    decompostion using a specific wavelet name wname or wavelet decompostion
    filters Lo_D and Hi_D. Such filters can be generated using wfilters.

    The global extension mode which can be change using dwtmode is used.

    The coefficient vector C contains the approximation coefficient at level N
    and all detail coefficient from level 1 to N

    The first entry of L is the length of the approximation coefficent,
    then the length of the detail coefficients are stored and the last
    value of L is the length of the signal vector.

    The approximation coefficient can be extracted with C(1:L(1)).
    The detail coefficients can be obtained with C(L(1):sum(L(1:2))),
    C(sum(L(1:2)):sum(L(1:3))),.... until C(sum(L(1:length(L)-2)):sum(L(1:length(L)-1))).
    Examples
    X = wnoise(4,10,0.5); //doppler with N=1024
    [C,L]=wavedec(X,3,'db2')
    """
    x = x.flatten()
    m1 = 1
    n1 = x.shape[0]
    if (len(args) == 1 and isinstance(args[0], str)):
        wname = args[0]
        ret = _wavelet_parser(wname.encode())
        filterLength = _wfilters_length(wname.encode())
        Lo_D, Hi_D = wfilters(wname,'d')
    elif(len(args) == 2):
        Lo_D = args[0]
        Hi_D = args[1]
        filterLength = Lo_D.shape[0]
    else:
        raise Exception("Wrong input!!")

    stride, val = _wave_len_validate(x.shape[0],filterLength)
    if (val == 0 or stride < N):
        raise Exception("Input signal is not valid for selected decompostion level and wavelets!")
    m4 = 1
    m5 = 1
    n4 = 0
    calLen = n1 * m1
    for count in np.arange(N):
        calLen += filterLength - 1
        temLen = np.floor(calLen/2).astype(int)
        n4 += temLen
        calLen = temLen

    n4 += temLen
    if (dwtmode("status","nodisp") == 'per'):
        n4 = 0
        calLen = n1 * m1
        for count in np.arange(N):
            # calLen += m3*n3 - 1;
            calLen = np.ceil(calLen/2.0).astype(int)
            temLen = calLen
            n4 += temLen
            # calLen = temLen;
        n4 += temLen
    n5 = N + 2

    output1 = np.zeros(n4*m4,dtype=np.float64)
    output2 = np.zeros(n5*m5,dtype=np.int32)
    _wave_dec_len_cal(filterLength, m1*n1, N, output2)
    _wavedec(x, output1, Lo_D, Hi_D, output2, N)
    return output1, output2


def waverec(C, L, *args):
    """
    Multiple level 1-D inverse discrete fast wavelet reconstruction
    Calling Sequence
    x0=waverec(C,L,wname)
    x0=waverec(C,L,Lo_R,Hi_R)
    Parameters
    wname : wavelet name, haar( "haar"), daubechies ("db1" to "db36"), coiflets ("coif1" to "coif17"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    x0 : reconstructed vector
    Lo_R : lowpass synthesis filter
    Hi_R : highpass synthesis filter
    C : coefficent array
    L : length array
    Description
    waverec can be used for multiple-level 1-D inverse discrete fast wavelet
    reconstruction.

    waverec supports only orthogonal or biorthogonal wavelets.
    Examples
    X = wnoise(4,10,0.5); //doppler with N=1024
    [C,L]=wavedec(X,3,'db2');
    x0=waverec(C,L,'db2');
    err = sum(abs(X-x0))
    """
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
