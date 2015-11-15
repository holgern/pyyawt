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
    n1_orig = 0
    if (np.size(a.shape) == 1):
        n1_orig = a.shape[0]
        a.shape = (1, n1_orig)
    m1,n1 = a.shape

    n2_orig = 0
    if (np.size(b.shape) == 1):
        n2_orig = b.shape[0]
        b.shape = (1, n2_orig)
    m2,n2 = b.shape
    m3 = 1
    n3 = m1 * n1 + m2 * n2 - 1
    Y = np.zeros((m3,n3),dtype=np.float64)
    _conv(a,b,Y)
    if (n1_orig > 0):
        Y = Y.flatten()
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
    output1 = np.zeros((m1,n1),dtype=np.float64,order="F")
    if (even_odd is None):
        _qmf_even(x.copy(order="F"), output1)
        return output1
    else:
        if ((even_odd % 2) == 0):
            _qmf_even(x.copy(order="F"), output1)
        else:
            _qmf_odd(x.copy(order="F"), output1)
        return output1.copy(order="C")


def dyaddown(x,*args):
    """
    dyadic downsampling

    Calling Sequence
    ----------------
    Y=dyaddown(x,[EVEN_ODD])
    Y=dyaddown(M,[EVEN_ODD],[type])
    Y=dyaddown(M,[type],[EVEN_ODD])

    Parameters
    ----------
    x : double vector
    M : double matrix
    EVEN_ODD : even or odd integer
    type : downsampling manner, 'r' for row, 'c' for column, and 'm' for row and column simutaneously.
    Y : downsampling result

    Description
    -----------
    dyaddown is an utility function for dyadic downsampling. if EVEN_ODD is even, even index entries of input will be kept. Otherwise, odd index entries will be kept. Default is even. Optional argumet type is especially for matrix input downsampling.

    Examples
    --------
    a=np.random.rand((1,100))
    Y=dyaddown(a)
    b=np.random.rand((25,25))
    Y=dyaddown(b,'r',0)
    """
    m_orig = 0
    n_orig = 0
    if (np.size(x.shape) == 2 and np.min(x.shape) == 1):
        (m_orig,n_orig) = x.shape
        x = x.flatten()
    if (len(args) == 0 and np.size(x.shape) == 1):
        m1 = 1
        n1 = x.shape[0]
        m2 = 1
        n2 = np.floor(n1 / 2).astype(int)
        output1 = np.zeros(n2*m2,dtype=np.float64)
        _dyaddown_1D_keep_even(x, output1)
        if (m_orig > 1):
            output1.shape = (n2, 1)
        elif (n_orig > 1):
            output1.shape = (1, n2)
        return output1
    elif (len(args) == 1 and np.size(x.shape) == 1):
        # isinstance(args[0], int)
        m1 = 1
        n1 = x.shape[0]
        if ((args[0] % 2) == 0):
            m3 = 1
            n3 = np.floor(n1 / 2).astype(int)
            output1 = np.zeros(n3*m3,dtype=np.float64)
            _dyaddown_1D_keep_even(x, output1)
        else:
            m3 = 1
            n3 = np.floor(n1 / 2).astype(int)
            if (n1 % 2 != 0):
                n3 += 1
            output1 = np.zeros(n3*m3,dtype=np.float64)
            _dyaddown_1D_keep_odd(x, output1)
        if (m_orig > 1):
            output1.shape = (n3, 1)
        elif (n_orig > 1):
            output1.shape = (1, n3)
        return output1
    elif (len(args) == 0 and np.size(x.shape) == 2):
        m1 = x.shape[0]
        n1 = x.shape[1]
        m2 = m1
        n2 = np.floor(n1 / 2).astype(int)
        output1 = np.zeros((m2,n2),dtype=np.float64,order="F")
        _dyaddown_2D_keep_even_col(x.copy(order="F"), output1)
        return output1.copy(order="C")
    elif (len(args) == 1 and np.size(x.shape) == 2 and isinstance(args[0], str)):
        m1 = x.shape[0]
        n1 = x.shape[1]
        if (args[0] == "r"):
            m3 = np.floor(m1 / 2).astype(int)
            n3 = n1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyaddown_2D_keep_even_row(x.copy(order="F"), output1)
        elif (args[0] == "c"):
            m3 = m1
            n3 = np.floor(n1 / 2).astype(int)
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyaddown_2D_keep_even_col(x.copy(order="F"), output1)
        elif (args[0] == "m"):
            m3 = np.floor(m1 / 2).astype(int)
            n3 = np.floor(n1 / 2).astype(int)
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyaddown_2D_keep_even(x.copy(order="F"), output1)
        else:
            raise Exception("Wrong input!!")
        return output1.copy(order="C")
    elif (len(args) == 1 and np.size(x.shape) == 2 and isinstance(args[0], int)):
        m1 = x.shape[0]
        n1 = x.shape[1]
        if ((args[0] % 2) == 0):
            m3 = m1
            n3 = np.floor(n1 / 2).astype(int)
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyaddown_2D_keep_even_col(x.copy(order="F"), output1)
        else:
            m3 = m1
            n3 = np.floor(n1 / 2).astype(int)
            if (n1 % 2 != 0):
                n3 += 1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyaddown_2D_keep_odd_col(x.copy(order="F"), output1)
        return output1.copy(order="C")
    elif (len(args) == 2 and np.size(x.shape) == 2):
        if (isinstance(args[0], int) and isinstance(args[1], str)):
            input_int = args[0]
            input_str = args[1]
        elif (isinstance(args[1], int) and isinstance(args[9], str)):
            input_int = args[1]
            input_str = args[0]
        else:
            raise Exception("Wrong input!!")
        m1 = x.shape[0]
        n1 = x.shape[1]
        if ((input_int % 2) == 0):
            if (input_str == "r"):
                m4 = np.floor(m1 / 2).astype(int)
                n4 = n1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_even_row(x.copy(order="F"), output1)
            elif (input_str == "c"):
                m4 = m1
                n4 = np.floor(n1 / 2).astype(int)
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_even_col(x.copy(order="F"), output1)
            elif (input_str == "m"):
                m4 = np.floor(m1 / 2).astype(int)
                n4 = np.floor(n1 / 2).astype(int)
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_even(x.copy(order="F"), output1)
            else:
                raise Exception("Wrong input!!")
        else:
            if (input_str == "r"):
                m4 = np.floor(m1 / 2).astype(int)
                n4 = n1
                if (m1 % 2 != 0):
                    m4 += 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_odd_row(x.copy(order="F"), output1)
            elif (input_str == "c"):
                m4 = m1
                n4 = np.floor(n1 / 2).astype(int)
                if (n1 % 2 != 0):
                    n4 += 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_odd_col(x.copy(order="F"), output1)
            elif (input_str == "m"):
                m4 = np.floor(m1 / 2).astype(int)
                n4 = np.floor(n1 / 2).astype(int)
                if (m1 % 2 != 0):
                    m4 += 1
                if (n1 % 2 != 0):
                    n4 += 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyaddown_2D_keep_odd(x.copy(order="F"), output1)
            else:
                raise Exception("Wrong input!!")
        return output1.copy(order="C")
    else:
        raise Exception("Wrong input!!")


def dyadup(x,*args):
    """
    dyadic upsampling

    Calling Sequence
    ----------------
    Y=dyadup(x,[EVEN_ODD])
    Y=dyadup(M,[EVEN_ODD],[type])
    Y=dyadup(M,[type],[EVEN_ODD])

    Parameters
    ----------
    x : double vector
    M : double matrix
    EVEN_ODD : even or odd integer
    type : upsampling manner, 'r' for row, 'c' for column, and 'm' for row and column simutaneously.
    Y : upsampling result

    Description
    -----------
    dyadup is an utility function for dyadic upsampling. if EVEN_ODD is even, zeors will be put between input entries and output length will be two times input length minus one. Otherwise, additional two zeros will be put at the head and tail of output so the output length will be two times input length plus one. Default is odd. Optional argumet type is especially for matrix input upsampling.

    Examples
    --------
    a=rand(1,100)
    Y=dyadup(a)
    b=rand(25,25)
    Y=dyadup(b,'r',0)
    """
    m_orig = 0
    n_orig = 0
    if (np.size(x.shape) == 2 and np.min(x.shape) == 1):
        (m_orig,n_orig) = x.shape
        x = x.flatten()
    if (len(args) == 0 and np.size(x.shape) == 1):
        m1 = 1
        n1 = x.shape[0]
        m2 = 1
        n2 = n1 * 2 + 1
        output1 = np.zeros(n2*m2,dtype=np.float64)
        _dyadup_1D_feed_even(x, output1)
        if (m_orig > 1):
            output1.shape = (n2, 1)
        elif (n_orig > 1):
            output1.shape = (1, n2)
        return output1
    elif (len(args) == 1 and np.size(x.shape) == 1):
        # isinstance(args[0], int)
        m1 = 1
        n1 = x.shape[0]
        if ((args[0] % 2) == 0):
            m3 = 1
            n3 = n1 * 2 - 1
            output1 = np.zeros(n3*m3,dtype=np.float64)
            _dyadup_1D_feed_odd(x, output1)
        else:
            m3 = 1
            n3 = n1 * 2 + 1
            output1 = np.zeros(n3*m3,dtype=np.float64)
            _dyadup_1D_feed_even(x, output1)
        if (m_orig > 1):
            output1.shape = (n3, 1)
        elif (n_orig > 1):
            output1.shape = (1, n3)
        return output1
    elif (len(args) == 0 and np.size(x.shape) == 2):
        m1 = x.shape[0]
        n1 = x.shape[1]
        m2 = m1
        n2 = n1 * 2 + 1
        output1 = np.zeros((m2,n2),dtype=np.float64,order="F")
        _dyadup_2D_feed_even_col(x.copy(order="F"), output1)
        return output1.copy(order="C")
    elif (len(args) == 1 and np.size(x.shape) == 2 and isinstance(args[0], str)):
        m1 = x.shape[0]
        n1 = x.shape[1]
        if (args[0] == "r"):
            m3 = m1 * 2 + 1
            n3 = n1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyadup_2D_feed_even_row(x.copy(order="F"), output1)
        elif (args[0] == "c"):
            m3 = m1
            n3 = n1 * 2 + 1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyadup_2D_feed_even_col(x.copy(order="F"), output1)
        elif (args[0] == "m"):
            m3 = m1 * 2 + 1
            n3 = n1 * 2 + 1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyadup_2D_feed_even(x.copy(order="F"), output1)
        else:
            raise Exception("Wrong input!!")
        return output1.copy(order="C")
    elif (len(args) == 1 and np.size(x.shape) == 2 and isinstance(args[0], int)):
        m1 = x.shape[0]
        n1 = x.shape[1]
        if ((args[0] % 2) == 0):
            m3 = m1
            n3 = n1 * 2 - 1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyadup_2D_feed_odd_col(x.copy(order="F"), output1)
        else:
            m3 = m1
            n3 = n1 * 2 + 1
            output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
            _dyadup_2D_feed_even_col(x.copy(order="F"), output1)
        return output1.copy(order="C")
    elif (len(args) == 2 and np.size(x.shape) == 2):
        if (isinstance(args[0], int) and isinstance(args[1], str)):
            input_int = args[0]
            input_str = args[1]
        elif (isinstance(args[1], int) and isinstance(args[9], str)):
            input_int = args[1]
            input_str = args[0]
        else:
            raise Exception("Wrong input!!")
        m1 = x.shape[0]
        n1 = x.shape[1]
        if ((input_int % 2) == 0):
            if (input_str == "r"):
                m4 = m1 * 2 - 1
                n4 = n1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_odd_row(x.copy(order="F"), output1)
            elif (input_str == "c"):
                m4 = m1
                n4 = n1 * 2 - 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_odd_col(x.copy(order="F"), output1)
            elif (input_str == "m"):
                m4 = m1 * 2 - 1
                n4 = n1 * 2 - 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_odd(x.copy(order="F"), output1)
            else:
                raise Exception("Wrong input!!")
        else:
            if (input_str == "r"):
                m4 = m1 * 2 + 1
                n4 = n1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_even_row(x.copy(order="F"), output1)
            elif (input_str == "c"):
                m4 = m1
                n4 = n1 * 2 + 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_even_col(x.copy(order="F"), output1)
            elif (input_str == "m"):
                m4 = m1 * 2 + 1
                n4 = n1 * 2 + 1
                output1 = np.zeros((m4,n4),dtype=np.float64,order="F")
                _dyadup_2D_feed_even(x.copy(order="F"), output1)
            else:
                raise Exception("Wrong input!!")
        return output1.copy(order="C")
    else:
        raise Exception("Wrong input!!")


def wkeep(x,*args):
    """
    signal extraction

    Calling Sequence
    ----------------
    Y=wkeep(x,L,[type])
    Y=wkeep(x,L,[FIRST])
    Y=wkeep(M,S,[indexVector])

    Parameters
    ----------
    x : double vector
    M : double matrix
    L : length integer
    type: extraction manner, 'l' for left, 'r' for right, and 'c' for center
    FIRST: index integer from which extraction starts.
    S : size integer vector containing row size and column size wanted
    indexVector : row and column index integer vector from which extraction starts.
    Y : extraction result

    Description
    -----------
    wkeep is an utility function for both vector and matrix extraction. For vector extraction, extractions will be aligned to the right, left or center based on optional argument type. So does matrix extraction.

    Examples
    --------
    a = np.linspace(1,8,8)
    X=np.dot(np.array([a]).T,np.array([a]))
    Y=wkeep(X,[4, 4])
    """
    m_orig = 0
    n_orig = 0
    if (np.size(x.shape) == 2 and np.min(x.shape) == 1):
        (m_orig,n_orig) = x.shape
        x = x.flatten()
    if (len(args) == 1 and np.size(x.shape) == 1 and isinstance(args[0], int)):
        m1 = 1
        n1 = x.shape[0]
        m3 = 1
        n3 = args[0]
        output1 = np.zeros(n3*m3,dtype=np.float64)
        _wkeep_1D_center(x,output1)
        if (m_orig > 1):
            output1.shape = (n3, 1)
        elif (n_orig > 1):
            output1.shape = (1, n3)
        return output1
    elif (len(args) == 2 and np.size(x.shape) == 1 and isinstance(args[0], int) and isinstance(args[1], str)):
        m1 = 1
        n1 = x.shape[0]
        m4 = 1
        n4 = args[0]
        output1 = np.zeros(n4*m4,dtype=np.float64)
        if (args[1] == 'l' or args[1] == 'L'):
            _wkeep_1D_left(x,output1)
        elif (args[1] == 'c' or args[1] == 'C'):
            _wkeep_1D_center(x,output1)
        elif (args[1] == 'r' or args[1] == 'R'):
            _wkeep_1D_right(x,output1)
        else:
            raise Exception("Wrong input!!")
        if (m_orig > 1):
            output1.shape = (n4, 1)
        elif (n_orig > 1):
            output1.shape = (1, n4)
        return output1
    elif (len(args) == 2 and np.size(x.shape) == 1 and isinstance(args[0], int) and isinstance(args[1], int)):
        m1 = 1
        n1 = x.shape[0]
        m4 = 1
        n4 = args[0]
        output1 = np.zeros(n4*m4,dtype=np.float64)
        _wkeep_1D_index(x,output1,args[1])
        if (m_orig > 1):
            output1.shape = (n4, 1)
        elif (n_orig > 1):
            output1.shape = (1, n4)
        return output1
    elif (len(args) == 1 and np.size(x.shape) == 2 and isinstance(args[0], list) and np.size(args[0]) == 2):
        m1 = x.shape[0]
        n1 = x.shape[1]
        m3 = args[0][0]
        n3 = args[0][1]
        output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
        _wkeep_2D_center(x.copy(order="F"), output1)
        return output1.copy(order="C")
    elif (len(args) == 2 and np.size(x.shape) == 2 and isinstance(args[0], list) and np.size(args[0]) == 2 and isinstance(args[1], list) and np.size(args[1]) == 2):
        m1 = x.shape[0]
        n1 = x.shape[1]
        m3 = args[0][0]
        n3 = args[0][1]
        output1 = np.zeros((m3,n3),dtype=np.float64,order="F")
        _wkeep_2D_index(x.copy(order="F"), output1, args[1][0], args[1][1])
        return output1.copy(order="C")
    else:
        raise Exception("Wrong Input!!")


def wextend(dim,extMethod,x,L,typeString=None):
    """
    signal extension

    Calling Sequence
    ----------------
    Y=wextend(onedim,extMode,x,L,[type])
    Y=wextend(twodim,extMode,M,sizeVector,[typeStringVector])
    Y=wextend(twodim,extMode,M,sizeVector,[typeString])
    Y=wextend(twodim,extMode,M,L)
    Y=wextend(row_col,extMode,M,L,[type])

    Parameters
    ----------
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
    -----------
    wextend is an utility function for signal extension.

    Examples
    --------
    a=rand(1,100);
    Y=wextend(1,'symh',a,5,'b');
    b=rand(25,25);
    Y=wextend(2,'symh',b,[3,5],'lb');
    Y=wextend('ar','symh',b,3,'r');
    """
    m_orig = 0
    n_orig = 0
    if (np.size(x.shape) == 2 and np.min(x.shape) == 1):
        (m_orig,n_orig) = x.shape
        x = x.flatten()
    if (np.size(x.shape) == 1 and typeString is not None):
        m3 = 1
        n3 = x.shape[0]
        if ((extMethod == 'per') and ((m3*n3) % 2 != 0)):
            if (typeString == 'l' or typeString == 'r'):
                m6 = 1
                n6 = n3 + L + 1
            elif (typeString == 'b'):
                m6 = 1
                n6 = n3 + 2*L + 1
        else:
            if (typeString == 'l' or typeString == 'r'):
                m6 = 1
                n6 = n3 + L
            elif (typeString == 'b'):
                m6 = 1
                n6 = n3 + 2*L
        output1 = np.zeros(n6*m6,dtype=np.float64)
        if (typeString == 'l'):
            _wextend_1D_left(x,output1,extMethod.encode())
        elif (typeString == 'r'):
            _wextend_1D_right(x,output1,extMethod.encode())
        else:
            _wextend_1D_center(x,output1,extMethod.encode())
        if (m_orig > 1):
            output1.shape = (n6, 1)
        elif (n_orig > 1):
            output1.shape = (1, n6)
        return output1
    elif (np.size(x.shape) == 1 and typeString is None):
        m3 = 1
        n3 = x.shape[0]
        if ((extMethod == 'per') and ((m3*n3) % 2 != 0)):
            m6 = 1
            n6 = n3 + 2*L + 1
        else:
            m6 = 1
            n6 = n3 + 2*L
        output1 = np.zeros(n6*m6,dtype=np.float64)
        _wextend_1D_center(x,output1,extMethod.encode())
        if (m_orig > 1):
            output1.shape = (n6, 1)
        elif (n_orig > 1):
            output1.shape = (1, n6)
        return output1
    elif (np.size(x.shape) == 2 and typeString is not None):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (m3 % 2 != 0)):
            if (typeString[0] == 'l' or typeString[0] == 'r'):
                m6 = m3 + L + 1
            elif (typeString[0] == 'r'):
                m6 = m3 + 2*L + 1
        else:
            if (typeString[0] == 'l' or typeString[0] == 'r'):
                m6 = m3 + L
            elif (typeString[0] == 'r'):
                m6 = m3 + 2*L
        if ((extMethod == 'per') and (n3 % 2 != 0)):
            if (typeString[1] == 'l' or typeString[1] == 'r'):
                n6 = n3 + L + 1
            elif (typeString[1] == 'r'):
                n6 = n3 + 2*L + 1
        else:
            if (typeString[1] == 'l' or typeString[1] == 'r'):
                n6 = n3 + L
            elif (typeString[1] == 'r'):
                n6 = n3 + 2*L
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        _wextend_2D(x.copy(order="F"), output1, extMethod.encode(), typeString[0].encode(), typeString[1].encode())
        return output1.copy(order="C")
    elif (np.size(x.shape) == 2 and typeString is None):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (m3 % 2 != 0)):
            m6 = m3 + 2*L + 1
        else:
            m6 = m3 + 2*L
        if ((extMethod == 'per') and (n3 % 2 != 0)):
            n6 = n3 + 2*L + 1
        else:
            n6 = n3 + 2*L
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        _wextend_2D(x.copy(order="F"), output1, extMethod.encode(), b'b', b'b')
        return output1.copy(order="C")
    elif (np.size(x.shape) == 2 and typeString is None and dim == "ar"):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (m3 % 2 != 0)):
            m6 = m3 + 2*L + 1
        else:
            m6 = m3 + 2*L
        if ((extMethod == 'per') and (n3 % 2 != 0)):
            n6 = n3 + 2*L + 1
        else:
            n6 = n3 + 2*L
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        wextend_2D_row(x.copy(order="F"), output1, extMethod.encode(), b'b')
        return output1.copy(order="C")
    elif (np.size(x.shape) == 2 and typeString is None and dim == "ac"):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (m3 % 2 != 0)):
            m6 = m3 + 2*L + 1
        else:
            m6 = m3 + 2*L
        if ((extMethod == 'per') and (n3 % 2 != 0)):
            n6 = n3 + 2*L + 1
        else:
            n6 = n3 + 2*L
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        wextend_2D_col(x.copy(order="F"), output1, extMethod.encode(), b'b')
        return output1.copy(order="C")
    elif (np.size(x.shape) == 2 and typeString is not None and dim == "ar"):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (m3 % 2 != 0)):
            if (typeString == 'l' or typeString == 'r'):
                m6 = m3 + L + 1
            elif (typeString == 'r'):
                m6 = m3 + 2*L + 1
        else:
            if (typeString == 'l' or typeString == 'r'):
                m6 = m3 + L
            elif (typeString == 'r'):
                m6 = m3 + 2*L
        n6 = n3
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        _wextend_2D_row(x.copy(order="F"), output1, extMethod.encode(), typeString.encode())
        return output1.copy(order="C")
    elif (np.size(x.shape) == 2 and typeString is not None and dim == "ac"):
        m3 = x.shape[0]
        n3 = x.shape[1]
        if ((extMethod == 'per') and (n3 % 2 != 0)):
            if (typeString == 'l' or typeString == 'r'):
                n6 = n3 + L + 1
            elif (typeString == 'r'):
                n6 = n3 + 2*L + 1
        else:
            if (typeString == 'l' or typeString == 'r'):
                n6 = n3 + L
            elif (typeString == 'r'):
                n6 = n3 + 2*L
        m6 = m3
        output1 = np.zeros((m6,n6),dtype=np.float64,order="F")
        _wextend_2D_col(x.copy(order="F"), output1, extMethod.encode(), typeString.encode())
        return output1.copy(order="C")
    else:
        raise Exception("Wrong input!!")


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
