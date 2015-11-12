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


def qmf(x,even_odd=None):
    """
    quadrature mirror

    Calling Sequence
    ---------------
    Y=qmf(x,[EVEN_ODD])

    Parameters
    ----------
    x: double vector
    EVEN_ODD: even or odd integer

    Returns
    -------
    Y: quadrature mirror

    Description
    -----------
    qmf is a quadrature mirror utility function on time domain. If EVEN_ODD is an even integer, output would be reversed version of input with even index entries sign changed. Otherwise, odd index entries will be changed. Default is even.

    Examples
    --------
    a=np.random.rand(3)
    Y=qmf(a)
    """
    if (np.size(x.shape) == 1):
        x.shape = (x.shape[0], 1)
    (m1, n1) = x.shape
    output1 = np.zeros((m1,n1),dtype=np.float64)
    if (even_odd is None):
        _qmf_even(x, output1)
        return output1
    else:
        if ((even_odd % 2) == 0):
            _qmf_even(x, output1)
        else:
            _qmf_odd(x, output1)
        return output1


def dyaddown(x,*args):
    """
    dyadic downsampling
    Calling Sequence
    Y=dyaddown(x,[EVEN_ODD])
    Y=dyaddown(M,[EVEN_ODD],[type])
    Y=dyaddown(M,[type],[EVEN_ODD])
    Parameters
    x : double vector
    M : double matrix
    EVEN_ODD : even or odd integer
    type : downsampling manner, 'r' for row, 'c' for column, and 'm' for row and column simutaneously.
    Y : downsampling result
    Description
    dyaddown is an utility function for dyadic downsampling. if EVEN_ODD is even, even index entries of input will be kept. Otherwise, odd index entries will be kept. Default is even. Optional argumet type is especially for matrix input downsampling.
    Examples
    a=np.random.rand((1,100))
    Y=dyaddown(a)
    b=np.random.rand((25,25))
    Y=dyaddown(b,'r',0)
    """
    raise Exception("Not yet implemented!!")


def dyadup(*args):
    raise Exception("Not yet implemented!!")


def wkeep(x,size,index=None):
    """
    signal extraction
    Calling Sequence
    Y=wkeep(x,[L],[type])
    Y=wkeep(x,[L],[FIRST])
    Y=wkeep(M,[S],[indexVector])
    Parameters
    x : double vector
    M : double matrix
    L : length integer
    type: extraction manner, 'l' for left, 'r' for right, and 'c' for center
    FIRST: index integer from which extraction starts.
    S : size integer vector containing row size and column size wanted
    indexVector : row and column index integer vector from which extraction starts.
    Y : extraction result
    Description
    wkeep is an utility function for both vector and matrix extraction. For vector extraction, extractions will be aligned to the right, left or center based on optional argument type. So does matrix extraction.
    Examples
    a = np.linspace(1,8,8)
    X=np.dot(np.array([a]).T,np.array([a]))
    Y=wkeep(X,[4, 4])
    """
    raise Exception("Not yet implemented!!")


def wextend(dim,extMode,x,size,typeString=None):
    """
    signal extension
    Calling Sequence
    Y=wextend(onedim,extMode,x,L,[type])
    Y=wextend(twodim,extMode,M,sizeVector,[typeStringVector])
    Y=wextend(twodim,extMode,M,sizeVector,[typeString])
    Y=wextend(twodim,extMode,M,L)
    Y=wextend(row_col,extMode,M,L,[type])
    Parameters
    x : double vector
    M : double matrix
    L : length integer
    type : extraction manner, 'l' for left, 'r' for right, and 'b' for both left and right
    sizeVector : integer vector containing row and column size to extend
    typeString : string for extension, 'bb', 'll', 'rr', 'bl', 'lb', 'br', 'rb', 'lr', 'rl'.
    typeStringVector : string vector for extension, ['b' 'b'], ['l' 'l'], ['r' 'r'], ['b' 'l'], ['l' 'b'], ['b' 'r'], ['r' 'b'], ['r' 'l'], ['l' 'r'].
    extMode : extension method, 'symh'('sym'), 'symw', 'asymh', 'asymw', 'zpd', 'zpd', 'per', 'ppd'.
    row_col : adding row or adding column, 'ar' or 'addrow' for row, 'ac' or 'addcol' for column.
    onedim : one dimension indication, 1, '1', '1d' and '1D'
    twodim : two dimension indication, 2, '2', '2d' and '2D'
    Y : extension result
    Description
    wextend is an utility function for signal extension.
    Examples
    a=rand(1,100);
    Y=wextend(1,'symh',a,5,'b');
    b=rand(25,25);
    Y=wextend(2,'symh',b,[3,5],'lb');
    Y=wextend('ar','symh',b,3,'r');    
    """
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
