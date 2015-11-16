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
    ----------------
    Calling Sequence
    ----------------
    [cA,cD]=dwt(x,wname,['mode',extMethod])
    [cA,cD]=dwt(x,Lo_D,Hi_D,['mode',extMethod])
    ----------
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
    -------
    Returns
    -------
    cA: array_like
         approximation coefficent
    cD: array_like
         detail coefficent
    --------
    Examples
    --------
    cA,cD=dwt(x,'db2','mode','asymh')
    """
    x = x.flatten()
    m1 = 1
    n1 = x.shape[0]
    ST = None
    if ((len(args) == 1 or (len(args) == 3 and isinstance(args[-1], str) and args[-2] == 'mode')) and isinstance(args[0], str)):
        wname = args[0]
        ret = _wavelet_parser(wname.encode())
        filterLength = _wfilters_length(wname.encode())
        Lo_D, Hi_D = wfilters(wname,'d')
        if (len(args) == 3 and isinstance(args[-1], str) and args[-2] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[-1])
    elif((len(args) == 2 or (len(args) == 4 and isinstance(args[-1], str) and args[-2] == 'mode')) and not isinstance(args[0], str)):
        Lo_D = args[0]
        Hi_D = args[1]
        filterLength = Lo_D.shape[0]
        if (len(args) == 4 and isinstance(args[-1], str) and args[-2] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[-1])
    else:
        raise Exception("Wrong input!!")

    stride, val = _wave_len_validate(x.shape[0],filterLength)
    if (val == 0):
        if (ST is not None):
            dwtmode(ST)
        raise Exception("Input signal is not valid for selected decompostion level and wavelets!")
    m3 = 1
    m4 = 1
    n3 = np.floor((m1*n1 + filterLength - 1)/2).astype(int)
    if (dwtmode("status","nodisp") == 'per'):
        n3 = np.ceil(((x.shape[0]))/2.0).astype(int)
    n4 = n3
    out1 = np.zeros(n3*m3,dtype=np.float64)
    out2 = np.zeros(n4*m4,dtype=np.float64)
    _dwt_neo(x,Lo_D,Hi_D,out1,out2)
    if (ST is not None):
        dwtmode(ST)
    return out1, out2


def idwt(cA, cD, *args):
    """
    Inverse Discrete Fast Wavelet Transform

    Calling Sequence
    ----------------
    X=idwt(cA,cD,wname,[L],['mode',extMethod])
    X=idwt(cA,cD,Lo_R,Hi_R,[L],['mode',extMethod])

    Parameters
    ----------
    wname: wavelet name, haar( "haar"), daubechies ("db1" to "db36"), coiflets ("coif1" to "coif17"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    x : reconstructed vector
    Lo_R: lowpass synthesis filter
    Hi_R: highpass syntheis filter
    L : restruction length
    cA: approximation coefficent
    cD: detail coefficent

    Description
    -----------
    idwt is for inverse discrete fast wavelet transform. Coefficent could be void vector as '[]' for cA or cD.

    Examples
    --------
    x=np.random.rand(1,100)
    [cA,cD]=dwt(x,'db2','mode','asymh')
    x0=idwt(cA,cD,'db2',100)
    """
    cA = cA.flatten()
    m1 = 1
    n1 = cA.shape[0]
    cD = cD.flatten()
    m2 = 1
    n2 = cD.shape[0]
    L = None
    ST = None
    if (len(args) >= 1 and (len(args) < 3 or (len(args) >= 3 and isinstance(args[-1], str) and args[-2] == 'mode')) and isinstance(args[0], str)):
        wname = args[0]
        ret = _wavelet_parser(wname.encode())
        filterLength = _wfilters_length(wname.encode())
        Lo_R, Hi_R = wfilters(wname,'r')
        if (len(args) > 1 and not isinstance(args[1],str)):
            L = args[1]
        if (len(args) >= 3 and isinstance(args[-1], str) and args[-2] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[-1])
    elif(len(args) >= 2 and (len(args) < 4 or (len(args) >= 4 and isinstance(args[-1], str) and args[-2] == 'mode')) and not isinstance(args[0], str)):
        Lo_R = args[0]
        Hi_R = args[1]
        filterLength = Lo_R.shape[0]
        if (len(args) > 2 and not isinstance(args[2],str)):
            L = args[2]
        if (len(args) >= 4 and isinstance(args[-1], str) and args[-2] == 'mode'):
            ST = dwtmode("status","nodisp")
            dwtmode(args[-1])
    else:
        raise Exception("Wrong input!!")

    stride, val = _wave_len_validate(np.max((m1*n1, m2*n2)), filterLength)
    if (val == 0):
        raise Exception("Input signal is not valid for selected decompostion level and wavelets!")
    m4 = 1
    if (L is None):
        if ((m1*n1 != 0)):  # and (dwtmode("status","nodisp") != "per")):
            n4 = m1*n1*2 - filterLength + 2
        if ((m1*n1 != 0) and (dwtmode("status","nodisp") == "per")):
            n4 = m1*n1*2
        if ((m1*n1 == 0)):  # and (dwtmode("status","nodisp") != "per")):
            n4 = m2*n2*2 - pWaveStruct.length + 2
        if ((m1*n1 == 0) and (dwtmode("status","nodisp") == "per")):
            n4 = m2*n2*2
    else:
        n4 = L
        if (L > np.max((m1*n1, m2*n2))*2 + filterLength):
            if (ST is not None):
                dwtmode(ST)
            raise Exception("Length Parameter is not valid for input vector!!")
    output1 = np.zeros(n4*m4,dtype=np.float64)
    if ((m1*n1 == 0) and (m2*n2 != 0)):
        _idwt_detail_neo(cD, Hi_R, output1)
    elif ((m1*n1 != 0) and (m2*n2 == 0)):
        _idwt_approx_neo(cA, Lo_R, output1)
    elif ((m1*n1 != 0) and (m2*n2 != 0)):
        _idwt_neo(cA, cD, Lo_R, Hi_R, output1)
    else:
        if (ST is not None):
            dwtmode(ST)
        raise Exception("Wrong input!!")
    if (ST is not None):
        dwtmode(ST)
    return output1


def wavedec(x, N, *args):
    """
    Multiple level 1-D discrete fast wavelet decomposition

    Calling Sequence
    ----------------
    [C,L]=wavedec(X,N,wname)
    [C,L]=wavedec(X,N,Lo_D,Hi_D)

    Parameters
    ----------
    wname : wavelet name, haar( "haar"), daubechies ("db1" to "db36"), coiflets ("coif1" to "coif17"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    X : signal vector
    N : decompostion level
    Lo_D : lowpass analysis filter
    Hi_D : highpass analysis filter
    C : coefficient vector
    L : length vector

    Description
    -----------
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
    C(sum(L(1:2)):sum(L(1:3))),.... until C(sum(L(1:length(L)-2)):sum(L(1:length(L)-1)))

    Examples
    --------
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
    ----------------
    x0=waverec(C,L,wname)
    x0=waverec(C,L,Lo_R,Hi_R)

    Parameters
    ----------
    wname : wavelet name, haar( "haar"), daubechies ("db1" to "db36"), coiflets ("coif1" to "coif17"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"
    x0 : reconstructed vector
    Lo_R : lowpass synthesis filter
    Hi_R : highpass synthesis filter
    C : coefficent array
    L : length array

    Description
    -----------
    waverec can be used for multiple-level 1-D inverse discrete fast wavelet
    reconstruction.

    waverec supports only orthogonal or biorthogonal wavelets.

    Examples
    --------
    X = wnoise(4,10,0.5); //doppler with N=1024
    [C,L]=wavedec(X,3,'db2');
    x0=waverec(C,L,'db2');
    err = sum(abs(X-x0))
    """
    C = C.flatten()
    m1 = 1
    n1 = C.shape[0]
    L = L.flatten()
    m2 = 1
    n2 = L.shape[0]
    L_summed_len = 0
    for count in np.arange(m2 * n2 - 1):
        L_summed_len += L[count]
    if (L_summed_len != m1*n1):
        raise Exception("Inputs are not coef and length array!!!")
    val = 0
    for count in np.arange(m2 * n2 - 1):
        if (L[count] > L[count+1]):
            val = 1
            break
    if (val != 0):
        raise Exception("Inputs are not coef and length array!!!")

    if (len(args) == 1 and isinstance(args[0], str)):
        wname = args[0]
        ret = _wavelet_parser(wname.encode())
        filterLength = _wfilters_length(wname.encode())
        Lo_R, Hi_R = wfilters(wname,'r')
    elif(len(args) == 2):
        Lo_R = args[0]
        Hi_R = args[1]
        filterLength = Lo_R.shape[0]
    else:
        raise Exception("Wrong input!!")
    if (L[0] < filterLength):
        raise Exception("Input signal is not valid for selected decompostion level and wavelets!\n")
    m4 = 1
    n4 = L[m2*n2-1]
    output1 = np.zeros(n4*m4,dtype=np.float64)
    _waverec(C, output1, Lo_R, Hi_R, L, m2*n2-2)
    return output1


def wrcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def appcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def detcoef(C, L, N=None):
    """
    1-D detail coefficients extraction

    Calling Sequence
    ----------------
    D=detcoef(C,L,[N])

    Parameters
    ----------
    D : reconstructed detail coefficient
    C : coefficent array
    L : length array
    N : restruction level with N<=length(L)-2
    Description
    -----------
    detcoef is for extraction of detail coeffient at different level
    after a multiple level decompostion. Extension mode is stored as
    a global variable and could be changed with dwtmode. If N is omitted,
    the detail coefficients will extract at the  maximum level (length(L)-2).

    The length of D depends on the level N.

    C and L can be generated using wavedec.

    Examples
    --------
    X = wnoise(4,10,0.5); //doppler with N=1024
    [C,L]=wavedec(X,3,'db2');
    D2=detcoef(C,L,2)
    """
    C = C.flatten()
    m1 = 1
    n1 = C.shape[0]
    L = L.flatten()
    m2 = 1
    n2 = L.shape[0]

    L_summed_len = 0
    for count in np.arange(m2 * n2 - 1):
        L_summed_len += L[count]
    if (L_summed_len != m1*n1):
        raise Exception("Inputs are not coef and length array!!!")
    val = 0
    for count in np.arange(m2 * n2 - 1):
        if (L[count] > L[count+1]):
            val = 1
            break
    if (val != 0):
        raise Exception("Inputs are not coef and length array!!!")
    if (N is None):
        m4 = 1
        n4 = L[0]
        N = m2*n2 - 2
    else:
        if ((N > m2*n2 - 2) or N < 1):
            raise Exception("Level Parameter is not valid for input vector!!!")
        m4 = 1
        n4 = L[n2*m2 - N - 1]
    output1 = np.zeros(n4*m4,dtype=np.float64)
    _detcoef(C,L,output1,m2*n2-2,N)
    return output1


def wenergy(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upcoef(cA, cD, *args):
    raise Exception("Not yet implemented!!")


def upwlev(cA, cD, *args):
    raise Exception("Not yet implemented!!")
