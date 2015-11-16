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

__all__ = ['orthfilt', 'biorfilt', 'dbwavf', 'coifwavf', 'symwavf', 'legdwavf', 'biorwavf', 'rbiorwavf', 'wfilters', 'wmaxlev', 'dwtmode']


def orthfilt(w):
    """
    orthfilt is an utility function for obtaining analysis and synthesis filter set of given orthogonal wavelets including haar, daubechies, coiflets and symlets

    Parameters
    ----------
    w: array_like
         scaling filter

    Returns
    -------
    Lo_D: array_like
         lowpass analysis filter
    Hi_D: array_like
         highpass analysis filter
    Lo_R: array_like
         lowpass synthesis filter
    Hi_R: array_like
         highpass synthesis filter

    Examples
    --------
    F = dbwavf("db2")
    [lo_d,hi_d,lo_r,hi_r]=orthfilt(F)
    """
    m1 = 1
    n1 = w.shape[0]
    Lo_D = np.zeros(n1,dtype=np.float64)
    Hi_D = np.zeros(n1,dtype=np.float64)
    Lo_R = np.zeros(n1,dtype=np.float64)
    Hi_R = np.zeros(n1,dtype=np.float64)
    _orthfilt(w,Lo_D,Hi_D,Lo_R,Hi_R)
    return Lo_D,Hi_D,Lo_R,Hi_R


def biorfilt(df,rf):
    """
    biorfilt is an utility function for obtaining analysis and synthesis filter set of given bi-orthogonal spline wavelets. DF and RF should be output of biorfilt result with the same length.

    Parameters
    ----------
    df: array_like
         analysis scaling filter
    rf: array_like
         synthesis scaling filter

    Returns
    -------
    Lo_D: array_like
         lowpass analysis filter
    Hi_D: array_like
         highpass analysis filter
    Lo_R: array_like
         lowpass synthesis filter
    Hi_R: array_like
         highpass synthesis filter

    Examples
    --------
    RF,DF = biorwavf('bior3.3')
    [lo_d,hi_d,lo_r,hi_r]=biorfilt(DF,RF)
    """
    m1 = 1
    n1 = df.shape[0]
    Lo_D = np.zeros(n1,dtype=np.float64)
    Hi_D = np.zeros(n1,dtype=np.float64)
    Lo_R = np.zeros(n1,dtype=np.float64)
    Hi_R = np.zeros(n1,dtype=np.float64)
    _biorfilt(df,rf,Lo_D,Hi_D,Lo_R,Hi_R)
    return Lo_D,Hi_D,Lo_R,Hi_R


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
    lowPass = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    _dbwavf(wname.encode(),b'Lo_R',lowPass)
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
    lowPass = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    _coifwavf(wname.encode(),b'Lo_R',lowPass)
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
    lowPass = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    _symwavf(wname.encode(),b'Lo_R',lowPass)
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
    lowPass = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    _legdwavf(wname.encode(),b'Lo_R',lowPass)
    return lowPass


def biorwavf(wname):
    """
    biorwavf is an utility function for obtaining twin scaling filters of bi-orthogonal spline wavelet including bior1.1, bior1.3, bior1.5, bior2.2, bior2.4, bior2.6, bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4, bior5.5 and bior6.8. Although the twin filters have different length, zeros has been fed to keep two filters the same length.

    Parameters
    ----------
    wname: str
         wavelet name,  'bior1.1' to 'bior6.8'

    Returns
    -------
    RF: array_like
         synthesis scaling filter
    DF: array_like
         analysis scaling filter

    Examples
    --------
    RF,DF = biorwavf('bior3.3');
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != PYYAWT_SPLINE_BIORTH):
        raise Exception("Wrong wavelet name!")
    filterLength = _wfilters_length(wname.encode())
    RF = np.zeros(filterLength, dtype=np.float64)
    DF = np.zeros(filterLength, dtype=np.float64)
    _biorwavf(wname.encode(),b'Lo_R',False,RF)
    _biorwavf(wname.encode(),b'Lo_D',True,DF)
    return RF,DF


def rbiorwavf(wname):
    """
    rbiorwavf is an utility function for obtaining twin scaling filters of bi-orthogonal spline wavelet including bior1.1, bior1.3, bior1.5, bior2.2, bior2.4, bior2.6, bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4, bior5.5 and bior6.8. Although the twin filters have different length, zeros has been fed to keep two filters the same length. rbiorwavf is reversing the results of biorwavf.

    Parameters
    ----------
    wname: str
         wavelet name,  'rbior1.1' to 'rbior6.8'

    Returns
    -------
    RF: array_like
         synthesis scaling filter
    DF: array_like
         analysis scaling filter

    Examples
    --------
    [RF,DF]=rbiorwavf('rbior3.3')
    """
    ret = _wavelet_parser(wname.encode())
    if (ret[0] != PYYAWT_SPLINE_RBIORTH):
        raise Exception("Wrong wavelet name!")
    RF = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    DF = np.zeros(_wfilters_length(wname.encode()),dtype=np.float64)
    _rbiorwavf(wname.encode(),b'Lo_R',False,RF)
    _rbiorwavf(wname.encode(),b'Lo_D',True,DF)
    return RF,DF


def wfilters(wname,filterType=None):
    """wfilters is an utility function for obtaining analysis and synthesis filter set.

    Calling Sequence
    ----------------

    [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname)
    [Lo_D,Hi_D]=wfilters(wname,'d')
    [Lo_R,Hi_R]=wfilters(wname,'r')
    [Lo_D,Lo_R]=wfilters(wname,'l')
    [Hi_D,Hi_R]=wfilters(wname,'h')

    Parameters
    ----------
    wname: str
         wavelet name,  wavelet name, haar( "haar"), daubechies ("db1" to "db20"), coiflets ("coif1" to "coif5"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"

    Returns
    -------
    Lo_D: array_like
         lowpass analysis filter
    Hi_D: array_like
         highpass analysis filter
    Lo_R: array_like
         lowpass synthesis filter
    Hi_R: array_like
         highpass synthesis filter

    Examples
    --------
    [lo_d,hi_d,lo_r,hi_r]=wfilters('db2')
    """
    ret_family, ret_member = _wavelet_parser(wname.encode())
    if (np.any(ret_family == np.array([PYYAWT_FARRAS, PYYAWT_KINGSBURYQ, PYYAWT_NOT_DEFINED]))):
        raise Exception("Wrong wavelet name!")
    filterLength = _wfilters_length(wname.encode())
    Lo_D = np.zeros(filterLength,dtype=np.float64)
    Hi_D = np.zeros(filterLength,dtype=np.float64)
    Lo_R = np.zeros(filterLength,dtype=np.float64)
    Hi_R = np.zeros(filterLength,dtype=np.float64)
    if (ret_family == PYYAWT_DAUBECHIES):
        _dbwavf(wname.encode(),b'Lo_D',Lo_D)
        _dbwavf(wname.encode(),b'Hi_D',Hi_D)
        _dbwavf(wname.encode(),b'Lo_R',Lo_R)
        _dbwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_COIFLETS):
        _coifwavf(wname.encode(),b'Lo_D',Lo_D)
        _coifwavf(wname.encode(),b'Hi_D',Hi_D)
        _coifwavf(wname.encode(),b'Lo_R',Lo_R)
        _coifwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_SYMLETS):
        _symwavf(wname.encode(),b'Lo_D',Lo_D)
        _symwavf(wname.encode(),b'Hi_D',Hi_D)
        _symwavf(wname.encode(),b'Lo_R',Lo_R)
        _symwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_SPLINE_BIORTH):
        _biorwavf(wname.encode(),b'Lo_D',False,Lo_D)
        _biorwavf(wname.encode(),b'Hi_D',False,Hi_D)
        _biorwavf(wname.encode(),b'Lo_R',False,Lo_R)
        _biorwavf(wname.encode(),b'Hi_R',False,Hi_R)
    elif (ret_family == PYYAWT_BEYLKIN):
        _beylkinwavf(wname.encode(),b'Lo_D',Lo_D)
        _beylkinwavf(wname.encode(),b'Hi_D',Hi_D)
        _beylkinwavf(wname.encode(),b'Lo_R',Lo_R)
        _beylkinwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_VAIDYANATHAN):
        _vaidyanathanwavf(wname.encode(),b'Lo_D',Lo_D)
        _vaidyanathanwavf(wname.encode(),b'Hi_D',Hi_D)
        _vaidyanathanwavf(wname.encode(),b'Lo_R',Lo_R)
        _vaidyanathanwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_DMEY):
        _dmeywavf(wname.encode(),b'Lo_D',Lo_D)
        _dmeywavf(wname.encode(),b'Hi_D',Hi_D)
        _dmeywavf(wname.encode(),b'Lo_R',Lo_R)
        _dmeywavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_BATHLETS):
        _bathletswavf(wname.encode(),b'Lo_D',Lo_D)
        _bathletswavf(wname.encode(),b'Hi_D',Hi_D)
        _bathletswavf(wname.encode(),b'Lo_R',Lo_R)
        _bathletswavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_LEGENDRE):
        _legendrewavf(wname.encode(),b'Lo_D',Lo_D)
        _legendrewavf(wname.encode(),b'Hi_D',Hi_D)
        _legendrewavf(wname.encode(),b'Lo_R',Lo_R)
        _legendrewavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_SPLINE_RBIORTH):
        _rbiorwavf(wname.encode(),b'Lo_D',False,Lo_D)
        _rbiorwavf(wname.encode(),b'Hi_D',False,Hi_D)
        _rbiorwavf(wname.encode(),b'Lo_R',False,Lo_R)
        _rbiorwavf(wname.encode(),b'Hi_R',False,Hi_R)
    elif (ret_family == PYYAWT_HAAR):
        _haarwavf(wname.encode(),b'Lo_D',Lo_D)
        _haarwavf(wname.encode(),b'Hi_D',Hi_D)
        _haarwavf(wname.encode(),b'Lo_R',Lo_R)
        _haarwavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_FARRAS):
        _farraswavf(wname.encode(),b'Lo_D',Lo_D)
        _farraswavf(wname.encode(),b'Hi_D',Hi_D)
        _farraswavf(wname.encode(),b'Lo_R',Lo_R)
        _farraswavf(wname.encode(),b'Hi_R',Hi_R)
    elif (ret_family == PYYAWT_KINGSBURYQ):
        _kingsburyqwavf(wname.encode(),b'Lo_D',Lo_D)
        _kingsburyqwavf(wname.encode(),b'Hi_D',Hi_D)
        _kingsburyqwavf(wname.encode(),b'Lo_R',Lo_R)
        _kingsburyqwavf(wname.encode(),b'Hi_R',Hi_R)
    # _wfilters(wname.encode(),Lo_D,Hi_D,Lo_R,Hi_R)
    if (filterType is None):
        flow = 1
        return Lo_D, Hi_D, Lo_R, Hi_R
    elif (filterType == 'd'):
        flow = 2
        return Lo_D, Hi_D
    elif (filterType == 'r'):
        flow = 3
        return Lo_R, Hi_R
    elif (filterType == 'l'):
        flow = 4
        return Lo_D, Lo_R
    elif (filterType == 'h'):
        flow = 5
        return Hi_D, Hi_R
    else:
        raise Exception("Wrong input!")


def wmaxlev(signalLength,wname):
    """
    wmaxlev is the maximum decompostion level calculation utility.

    Parameters
    ----------
    signalLength: int
          signal length

    wname: str
          wavelet name,  wavelet name, haar( "haar"), daubechies ("db1" to "db20"), coiflets ("coif1" to "coif5"), symlets ("sym2" to "sym20"), legendre ("leg1" to "leg9"), bathlets("bath4.0" to "bath4.15" and "bath6.0" to "bath6.15"), dmey ("dmey"), beyklin ("beylkin"), vaidyanathan ("vaidyanathan"), biorthogonal B-spline wavelets ("bior1.1" to "bior6.8"), "rbior1.1" to "rbior6.8"

    Returns
    -------
    L: int
          decomposition level

    Examples
    --------
    L=wmaxlev(100,'db5')
    """

    if (np.size(signalLength) == 1):
        filterLength = _wfilters_length(wname.encode())
        stride, val = _wave_len_validate(signalLength,filterLength)
        if (val == 0):
            raise Exception("Unrecognized Input Pattern or parameter not valid for the algorithm! Please refer to help pages!\n")
        return stride
    elif (np.size(signalLength) == 2):
        filterLength = _wfilters_length(wname.encode())
        stride1, val1 = _wave_len_validate(signalLength[0],filterLength)
        if (val1 == 0):
            raise Exception("The wavelet you select is not appropriate for that row size of the matrix!\n")
        stride2, val2 = _wave_len_validate(signalLength[1],filterLength)
        if (val2 == 0):
            raise Exception("The wavelet you select is not appropriate for that column size of the matrix!\n")

        return np.min([stride1, stride2])
    else:
        return 0


def dwtmode(mode=None, nodisp=None):
    """
    dwtmode is to display or change extension mode.

    Parameters
    ----------
    mode: str
           'symh'('sym'), 'symw', 'asymh', 'asymw', 'zpd', 'zpd', 'per', 'ppd'.

    nodisp: str
          None or 'nodisp'

    Returns
    -------
    mode: str
          extension mode

    Examples
    --------
    dwtmode()
    dwtmode('status')
    mode=dwtmode('status','nodisp')
    dwtmode(mode)

    """

    mode_strings = np.array(['zpd','symh', 'symw', 'asymh', 'asymw', 'sp0', 'sp1', 'ppd', 'per'])
    modes = np.array([PYYAWT_ZPD, PYYAWT_SYMH, PYYAWT_SYMW, PYYAWT_ASYMH, PYYAWT_ASYMW, PYYAWT_SP0,
                          PYYAWT_SP1, PYYAWT_PPD, PYYAWT_PER])
    if (mode is None or mode == "status" and nodisp is None):
        dwtmode = _getdwtMode()
        if (dwtmode == PYYAWT_ZPD):
            print("**     DWT Extension Mode: Zero Padding   **\n")
        elif (dwtmode == PYYAWT_SYMH):
            print("**     DWT Extension Mode: Half Symmetrization  **\n")
        elif (dwtmode == PYYAWT_SYMW):
            print("**     DWT Extension Mode: Whole Symmetrization **\n")
        elif (dwtmode == PYYAWT_ASYMH):
            print("**     DWT Extension Mode: Half Asymmetrization   **\n")
        elif (dwtmode == PYYAWT_ASYMW):
            print("**     DWT Extension Mode: Whole Asymmetrization   **\n")
        elif (dwtmode == PYYAWT_SP0):
            print("**     DWT Extension Mode: order 0 smooth padding  **\n")
        elif (dwtmode == PYYAWT_SP1):
            print("**     DWT Extension Mode: order 1 smooth padding  **\n")
        elif (dwtmode == PYYAWT_PPD):
            print("**     DWT Extension Mode: Periodic Padding**\n")
        elif (dwtmode == PYYAWT_PER):
            print("**     DWT Extension Mode: Periodization    **\n")
    elif (np.any(mode == modes)):
        errCode = _dwtWrite(mode_strings[np.where(mode == modes)[0]].encode())
        if (errCode > 0):
            raise Exception("DWT Extension Mode error")
    elif (np.any(mode == mode_strings)):
        errCode = _dwtWrite(mode.encode())
        if (errCode > 0):
            raise Exception("DWT Extension Mode error")
    elif (mode == "status" and nodisp == "nodisp"):
        dwtmode = _getdwtMode()
        modestr = ""
        if (dwtmode == PYYAWT_ZPD):
            modestr = 'zpd'
        elif (dwtmode == PYYAWT_SYMH):
            modestr = 'symh'
        elif (dwtmode == PYYAWT_SYMW):
            modestr = 'symw'
        elif (dwtmode == PYYAWT_ASYMH):
            modestr = 'asymh'
        elif (dwtmode == PYYAWT_ASYMW):
            modestr = 'asymw'
        elif (dwtmode == PYYAWT_SP0):
            modestr = 'sp0'
        elif (dwtmode == PYYAWT_SP1):
            modestr = 'sp1'
        elif (dwtmode == PYYAWT_PPD):
            modestr = 'ppd'
        elif (dwtmode == PYYAWT_PER):
            modestr = 'per'
        return modestr
