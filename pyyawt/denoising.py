# -*- coding: utf-8 -*-

# Copyright (c) 2015 Holger Nahrstaedt
# See COPYING for license details.

"""
dwt1d function for wavelet
"""

from __future__ import division, print_function, absolute_import
import numpy as np
import sys as sys
from ._pyyawt import *
from .dwt1d import *

__all__ = ['wnoisest', 'wden', 'thselect', 'ValSUREThresh', 'dyadlength', 'wthresh', 'wnoise']


def wnoisest(C,L=None,S=None):
    """
    estimates of the detail coefficients' standard deviation for levels contained in the input vector S

    Parameters
    ----------
    C: array_like
         coefficent array
    L: array_like
         coefficent array
    S: array_like
         estimate noise for this decompostion levels
    Returns
    -------
    STDC: array_like
         STDC[k] is an estimate of the standard deviation of C[k]

    Examples
    --------
    [c,l] = wavedec(x,2,'db3')
    wnoisest(c,l,[0,1])
    """
    if (S is None and L is None):
        if (isinstance(C, list)):
            STDC = zeros(1,length(C))
            for k in np.arange(len(C)):
                STDC[k] = np.median(np.abs(C[k]))/0.6745
        else:
            STDC = np.median(np.abs(C))/0.6745
        return STDC

    maxLevel = np.size(L)-2
    if (S is None):
        S = np.arange(maxLevel)

    if(maxLevel < np.size(S)):
        print("C,L does not contain so much levels. reduce S!")
    STDC = np.zeros(np.size(S))
    for level in np.arange(np.size(S)):

        # STDC(level) = median(abs(C(sum(l(1:(maxLevel-level+1)))+1:sum(l(1:(maxLevel-level+2))))))/.6745;
        # STDC[level] = np.median(np.abs(C[maxLevel-S[level]]))/.6745
        STDC[level] = np.median(np.abs(detcoef(C,L,level + 1)))/.6745
        # STDC(level) = median(abs(C(sum(L(1:level))+(1:L(level+1)))))/.6745;

    return STDC


def wden(*args):
    """
    wden performs an automatic de-noising process of a one-dimensional signal using wavelets.
    Calling Sequence
    ----------------
    [XD,CXD,LXD] = wden(X,TPTR,SORH,SCAL,N,wname)
    [XD,CXD,LXD] = wden(C,L,TPTR,SORH,SCAL,N,wname)
    Parameters
    ----------
    x: array_like
          input vector
    C: array_like
         coefficent array
    L: array_like
         coefficent array
    TPTR: str threshold selection rule
         'rigrsure' uses the principle of Stein's Unbiased Risk.
         'heursure' is an heuristic variant of the first option.
         'sqtwolog' for universal threshold
         'minimaxi' for minimax thresholding
    SORH: str
         ('s' or 'h') soft or hard thresholding
    SCAL: str
         'one' for no rescaling
         'sln' for rescaling using a single estimation of level noise based on first-level coefficients
         'mln' for rescaling done using level-dependent estimation of level noise
    N: int
          N: decompostion level
    wname: str
          wavelet name
    Returns
    -------
    XD: array_like
         de-noised signal
    CXD: array_like
         de-noised coefficent array
    LXD: array_like
         de-noised length array

    Examples
    --------
    [xref,x] = wnoise(3,11,3)
    level = 4
    xd = wden(x,'heursure','s','one',level,'sym8')
    """
    if (len(args) == 6):
        x = args[0]
        TPTR = args[1]
        SORH = args[2]
        SCAL = args[3]
        N = args[4]
        wname = args[5]
        if (not isinstance(wname, str)):
            raise Exception("wname must be a string")
        C,L = wavedec(x,N,wname)
    elif (len(args) == 6):
        C = args[0]
        L = args[1]
        TPTR = args[2]
        SORH = args[3]
        SCAL = args[4]
        N = args[5]
        wname = args[6]
        if (not isinstance(wname, str)):
            raise Exception("wname must be a string")
    else:
        raise Exception("Wrong input!!")
    if (not isinstance(SCAL, str)):
        raise Exception("SCAL must be a string")
    if (not isinstance(SORH, str)):
        raise Exception("SORH must be a string")
    if (not isinstance(TPTR, str)):
        raise Exception("TPTR must be a string")

    if (SCAL == 'one'):
        sigma = np.ones(N)
    elif (SCAL == 'sln'):
        sigma = np.ones(N)*wnoisest(C,L,np.array([0]))
    elif (SCAL == 'mln'):
        sigma = wnoisest(C,L,np.arange(N))
    else:
        raise Exception("SCAL must be either ''one'',''sln'' or ''mln''")

    D = []
    for n in np.arange(N):
        D.append(detcoef(C,L,n + 1))
    CXD = C
    LXD = L
    i = 0
    for n in np.arange(N,0,-1)-1:
        i = i+1
        CXD[np.sum(L[:i])+np.arange(L[i])] = wthresh(D[n],SORH,thselect(D[n]/sigma[n],TPTR)*sigma[n])

    # for n in np.arange(N):
    #    CXD[N-n] = wthresh(C[N-n],SORH,thselect(C[N-n]/sigma[n],TPTR)*sigma[n])

    XD = waverec(CXD, LXD, wname)
    return XD,CXD,LXD


def thselect(X,TPTR):
    """
    Threshold selection for de-noising. The algorithm works only if the signal X has a white noise of N(0,1). Dealing with unscaled or nonwhite noise can be handled using rescaling of the threshold.

    Parameters
    ----------
    X: array
         input vector with scaled white noise (N(0,1))
    TPTR: str
         'rigrsure': adaptive threshold selection using principle of Stein's Unbiased Risk Estimate.
         'heursure': heuristic variant of the first option.
         'sqtwolog': threshold is sqrt(2*log(length(X))).
         'minimaxi': minimax thresholding.


    Returns
    -------
    THR: float
         threshold X-adapted value using selection rule defined by string TPTR

    Examples
    --------
    x = np.random.randn(1000)
    thr = thselect(x,'rigrsure')
    """
    if (TPTR == 'rigrsure'):
        THR = ValSUREThresh(X)
    elif (TPTR == 'heursure'):
        n,j = dyadlength(X)
        magic = np.sqrt(2*np.log(n))
        eta = (np.linalg.norm(X)**2 - n)/n
        crit = j**(1.5)/np.sqrt(n)
        if (eta < crit):
            THR = magic
        else:
            THR = np.min((ValSUREThresh(X), magic))
    elif (TPTR == 'sqtwolog'):
        n,j = dyadlength(X)
        THR = np.sqrt(2*np.log(n))
    elif (TPTR == 'minimaxi'):
        lamlist = [0, 0, 0, 0, 0, 1.27, 1.474, 1.669, 1.860, 2.048, 2.232, 2.414, 2.594, 2.773, 2.952, 3.131, 3.310, 3.49, 3.67, 3.85, 4.03, 4.21]
        n,j = dyadlength(X)
        if(j <= np.size(lamlist)):
            THR = lamlist[j-1]
        else:
            THR = 4.21 + (j-np.size(lamlist))*0.18
    return THR


def ValSUREThresh(X):
    """
    Adaptive Threshold Selection Using Principle of SURE

    Parameters
    ----------
    X: array
         Noisy Data with Std. Deviation = 1

    Returns
    -------
    tresh: float
         Value of Threshold

    """
    n = np.size(X)

    # a = mtlb_sort(abs(X)).^2
    a = np.sort(np.abs(X))**2

    c = np.linspace(n-1,0,n)
    s = np.cumsum(a)+c*a
    risk = (n - (2 * np.arange(n)) + s)/n
    # risk = (n-(2*(1:n))+(cumsum(a,'m')+c(:).*a))/n;
    ibest = np.argmin(risk)
    THR = np.sqrt(a[ibest])
    return THR


def dyadlength(x):
    """
    Find length and dyadic length of array

    Parameters
    ----------
    X: array
         array of length n = 2^J (hopefully)

    Returns
    -------
    n: int
         length(x)
    J: int
         least power of two greater than n
    """
    n = np.size(x)
    J = np.ceil(np.log(n)/np.log(2)).astype(np.int)
    return n,J


def wthresh(X,SORH,T):
    """
    doing either hard (if SORH = 'h') or soft (if SORH = 's') thresholding

    Parameters
    ----------
    X: array
         input data (vector or matrix)
    SORH: str
         's': soft thresholding
         'h' : hard thresholding
    T: float
          threshold value

    Returns
    -------
    Y: array_like
         output

    Examples
    --------
    y = np.linspace(-1,1,100)
    thr = 0.4
    ythard = wthresh(y,'h',thr)
    ytsoft = wthresh(y,'s',thr)
    """
    if ((SORH) != 'h' and (SORH) != 's'):
        print(' SORH must be either h or s')

    elif (SORH == 'h'):
        Y = X * (np.abs(X) > T)
        return Y
    elif (SORH == 's'):
        res = (np.abs(X) - T)
        res = (res + np.abs(res))/2.
        Y = np.sign(X)*res
        return Y


def wnoise(FUN,N,SQRT_SNR=1):
    """
    Noisy wavelet test data

    Parameters
    ----------
    FUN: str / int
        1 or 'blocks'
        2 or 'bumps'
        3 or 'heavy sine'
        4 or 'doppler'
        5 or 'quadchirp'
        6 or 'mishmash'
    N: int
         vector length of X = 2^N
    SQRT_SNR: float
          standard deviation of added noise

    Returns
    -------
    X: array_like
         test data
    XN: array_like
         noisy test data (rand(1,N,'normal') is added!)

    Examples
    --------
    [x,noisyx] = wnoise(4,10,7);
    """
    N = 2**N
    X = np.zeros(N)
    if (isinstance(FUN, str)):
        if (FUN == 'blocks'):
            FUN = 1
        elif(FUN == 'bumps'):
            FUN = 2
        elif(FUN == 'heavy sine'):
            FUN = 3
        elif(FUN == 'doppler'):
            FUN = 4
        elif(FUN == 'quadchirp'):
            FUN = 5
        elif(FUN == 'mishmash'):
            FUN = 6
    if (FUN == 1):
        tj = np.array([.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81])
        hj = np.array([4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2])
        for n in np.arange(N):
            t = (n-1)/(N-1)
            X[n] = np.sum(hj*((1+np.sign(t-tj))/2.))
    elif (FUN == 2):
        tj = np.array([.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81])
        hj = np.array([4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2])
        wj = np.array([.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005])
        for n in np.arange(N):
            t = (n-1)/(N-1)
            X[n] = np.sum(hj*((1+np.abs((t-tj)/wj)**4.)**(-1)))
    elif (FUN == 3):
        for n in np.arange(N):
            t = (n-1)/(N-1)
            X[n] = 4*np.sin(4*np.pi*t)-np.sign(t-0.3)-np.sign(.72-t)
    elif (FUN == 4):
        eps = 0.05
        for n in np.arange(N):
            t = (n-1)/(N-1)
            X[n] = np.sqrt(t*(1-t))*np.sin(2*np.pi*(1-eps)/(t+eps))
    elif (FUN == 5):
        t = np.arange(N)/N
        X = np.sin((np.pi/3.) * t * (N * t**2.))
    elif (FUN == 6):
        t = np.arange(N)/N
        X = np.sin((np.pi/3) * t * (N * t**2))
        X = X + np.sin(np.pi * (N * .6902) * t)
        X = X + np.sin(np.pi * t * (N * .125 * t))
    if (SQRT_SNR != 1):
        X = X/np.std(X)*SQRT_SNR
    XN = X+np.random.randn(N)

    return X,XN
