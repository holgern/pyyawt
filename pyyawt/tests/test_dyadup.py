# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDyadup(unittest.TestCase):
    def test_vector_odd(self):
        a = np.random.randn(51)
        ind = np.arange(51)
        b = pyyawt.dyadup(a)
        b1 = pyyawt.dyadup(a,0)
        b2 = pyyawt.dyadup(a,1)
        c = np.zeros(103)
        c[2*ind+1] = a[ind]
        c1 = np.zeros(101)
        c1[2*ind] = a[ind]

        np.testing.assert_almost_equal(b, b2)
        np.testing.assert_almost_equal(b1, c1)
        np.testing.assert_almost_equal(b2, c)

        a = np.random.randn(1,51)
        ind = np.arange(51)
        b = pyyawt.dyadup(a)
        b1 = pyyawt.dyadup(a,0)
        b2 = pyyawt.dyadup(a,1)
        b3 = pyyawt.dyadup(a.T)
        b4 = pyyawt.dyadup(a.T,0)
        b5 = pyyawt.dyadup(a.T,1)
        c = np.zeros((1,103))
        c[:,2*ind+1] = a[:,ind]
        c1 = np.zeros((1,101))
        c1[:,2*ind] = a[:,ind]

        np.testing.assert_almost_equal(b, b2)
        np.testing.assert_almost_equal(b1, c1)
        np.testing.assert_almost_equal(b2, c)
        np.testing.assert_almost_equal(b3, b2.T)
        np.testing.assert_almost_equal(b4, c1.T)
        np.testing.assert_almost_equal(b5, c.T)

    def test_vector_even(self):
        a = np.random.randn(50)
        ind = np.arange(50)
        b = pyyawt.dyadup(a)
        b1 = pyyawt.dyadup(a,0)
        b2 = pyyawt.dyadup(a,1)
        c = np.zeros(101)
        c[2*ind+1] = a[ind]
        c1 = np.zeros(99)
        c1[2*ind] = a[ind]

        np.testing.assert_almost_equal(b, b2)
        np.testing.assert_almost_equal(b1, c1)
        np.testing.assert_almost_equal(b2, c)

        a = np.random.randn(1,50)
        ind = np.arange(50)
        b = pyyawt.dyadup(a)
        b1 = pyyawt.dyadup(a,0)
        b2 = pyyawt.dyadup(a,1)
        b3 = pyyawt.dyadup(a.T)
        b4 = pyyawt.dyadup(a.T,0)
        b5 = pyyawt.dyadup(a.T,1)
        c = np.zeros((1,101))
        c[:,2*ind+1] = a[:,ind]
        c1 = np.zeros((1,99))
        c1[:,2*ind] = a[:,ind]

        np.testing.assert_almost_equal(b, b2)
        np.testing.assert_almost_equal(b1, c1)
        np.testing.assert_almost_equal(b2, c)
        np.testing.assert_almost_equal(b3, b2.T)
        np.testing.assert_almost_equal(b4, c1.T)
        np.testing.assert_almost_equal(b5, c.T)

    def test_matrix(self):
        a = np.random.randn(50,51)
        b = pyyawt.dyadup(a)
        b1 = pyyawt.dyadup(a,0)
        b2 = pyyawt.dyadup(a,1)
        b3 = pyyawt.dyadup(a,'r')
        b4 = pyyawt.dyadup(a,'c')
        b5 = pyyawt.dyadup(a,'m')
        b6 = pyyawt.dyadup(a,0,'r')
        b7 = pyyawt.dyadup(a,0,'c')
        b8 = pyyawt.dyadup(a,0,'m')
        b9 = pyyawt.dyadup(a,1,'r')
        b10 = pyyawt.dyadup(a,1,'c')
        b11 = pyyawt.dyadup(a,1,'m')
        ind = np.arange(50)
        ind1 = np.arange(51)
        c = np.zeros((99,51))
        c[2*ind,:] = a[ind,:]
        c1 = np.zeros((101,51))
        c1[2*ind+1,:] = a[ind,:]

        c2 = np.zeros((50,101))
        c2[:,2*ind1] = a[:,ind1]
        c3 = np.zeros((50,103))
        c3[:,2*ind1+1] = a[:,ind1]

        c4 = np.zeros((99,101))
        c4[2*ind,:][:,2*ind1] = a[ind,:][:,ind1]
        c5 = np.zeros((101,103))
        c5[2*ind+1,:][:,2*ind1+1] = a[ind,:][:,ind1]
        np.testing.assert_almost_equal(b, b10)
        np.testing.assert_almost_equal(b, c3)
        np.testing.assert_almost_equal(b1, c2)
        np.testing.assert_almost_equal(b2, b)
        np.testing.assert_almost_equal(b3, c1)
        np.testing.assert_almost_equal(b4, b2)
        # np.testing.assert_almost_equal(b5, c5)
        np.testing.assert_almost_equal(b6, c)
        np.testing.assert_almost_equal(b7, c2)
        # np.testing.assert_almost_equal(b8, c4)
        np.testing.assert_almost_equal(b9, c1)
        np.testing.assert_almost_equal(b10, c3)
        np.testing.assert_almost_equal(b11, b5)
