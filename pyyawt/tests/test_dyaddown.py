# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDyaddown(unittest.TestCase):
    def test_vector_odd(self):
        a = np.random.randn(51)
        ind = np.arange(25)
        ind1 = np.arange(26)
        b = pyyawt.dyaddown(a)
        b1 = pyyawt.dyaddown(a,0)
        b2 = pyyawt.dyaddown(a,1)

        c = a[2*ind+1]
        c1 = a[2*ind1]

        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)

        a = np.random.randn(1,51)
        ind = np.arange(25)
        ind1 = np.arange(26)
        b = pyyawt.dyaddown(a)
        b1 = pyyawt.dyaddown(a,0)
        b2 = pyyawt.dyaddown(a,1)

        b3 = pyyawt.dyaddown(a.T)
        b4 = pyyawt.dyaddown(a.T,0)
        b5 = pyyawt.dyaddown(a.T,1)
        c = a[:,2*ind+1]
        c1 = a[:,2*ind1]

        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, b4)
        np.testing.assert_almost_equal(b3, b.T)
        np.testing.assert_almost_equal(b5, b2.T)

    def test_vector_even(self):
        a = np.random.randn(50)
        ind = np.arange(25)
        b = pyyawt.dyaddown(a)
        b1 = pyyawt.dyaddown(a,0)
        b2 = pyyawt.dyaddown(a,1)

        c = a[2*ind+1]
        c1 = a[2*ind]

        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)

        a = np.random.randn(1,50)
        ind = np.arange(25)
        b = pyyawt.dyaddown(a)
        b1 = pyyawt.dyaddown(a,0)
        b2 = pyyawt.dyaddown(a,1)

        b3 = pyyawt.dyaddown(a.T)
        b4 = pyyawt.dyaddown(a.T,0)
        b5 = pyyawt.dyaddown(a.T,1)
        c = a[:,2*ind+1]
        c1 = a[:,2*ind]

        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, b4)
        np.testing.assert_almost_equal(b3, b.T)
        np.testing.assert_almost_equal(b5, b2.T)

    def test_matrix(self):
        a = np.random.randn(50,51)
        b = pyyawt.dyaddown(a)
        b1 = pyyawt.dyaddown(a,0)
        b2 = pyyawt.dyaddown(a,1)
        b3 = pyyawt.dyaddown(a, 'r')
        b4 = pyyawt.dyaddown(a, 'c')
        b5 = pyyawt.dyaddown(a, 'm')
        b6 = pyyawt.dyaddown(a, 0, 'r')
        b7 = pyyawt.dyaddown(a, 0, 'c')
        b8 = pyyawt.dyaddown(a, 0, 'm')
        b9 = pyyawt.dyaddown(a, 1, 'r')
        b10 = pyyawt.dyaddown(a, 1, 'c')
        b11 = pyyawt.dyaddown(a, 1, 'm')
        ind = np.arange(25)
        ind1 = np.arange(26)
        c = a[2*ind+1,:]
        c1 = a[2*ind,:]
        c2 = a[:, 2*ind+1]
        c3 = a[:,2*ind1]
        c4 = a[2*ind+1,:][:,2*ind+1]
        c5 = a[2*ind,:][:,2*ind+1]
        c7 = a[2*ind,:][:,2*ind1]

        np.testing.assert_almost_equal(b, c2)
        np.testing.assert_almost_equal(b1, b)
        np.testing.assert_almost_equal(b2, c3)
        np.testing.assert_almost_equal(b3, c)
        np.testing.assert_almost_equal(b4, c2)
        np.testing.assert_almost_equal(b5, c4)
        np.testing.assert_almost_equal(b6, c)
        np.testing.assert_almost_equal(b7, c2)
        np.testing.assert_almost_equal(b8, c4)
        np.testing.assert_almost_equal(b9, c1)
        np.testing.assert_almost_equal(b10, c3)
        np.testing.assert_almost_equal(b11, c7)
