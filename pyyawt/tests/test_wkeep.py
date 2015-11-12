# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestWkeep(unittest.TestCase):
    def test_vector_even(self):
        a = np.random.randn(50)
        b = pyyawt.wkeep(a,27)
        b1 = pyyawt.wkeep(a,27,'c')
        b2 = pyyawt.wkeep(a,27,'l')
        b3 = pyyawt.wkeep(a,27,'r')
        b4 = pyyawt.wkeep(a,27)
        b5 = pyyawt.wkeep(a,27,'c')
        b6 = pyyawt.wkeep(a,27,'l')
        b7 = pyyawt.wkeep(a,27,'r')
        b8 = pyyawt.wkeep(a,27,2)
        b9 = pyyawt.wkeep(a,27,2)
        c = a[11:38]
        c1 = a[0:27]
        c2 = a[-27:]
        c3 = a[1:28]
        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b1, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)
        np.testing.assert_almost_equal(b4, b1)
        np.testing.assert_almost_equal(b5, b1)
        np.testing.assert_almost_equal(b6, b2)
        np.testing.assert_almost_equal(b7, b3)
        np.testing.assert_almost_equal(b8, c3)
        np.testing.assert_almost_equal(b9, c3)

        a = np.random.randn(1,50)
        b = pyyawt.wkeep(a,27)
        b1 = pyyawt.wkeep(a,27,'c')
        b2 = pyyawt.wkeep(a,27,'l')
        b3 = pyyawt.wkeep(a,27,'r')
        b4 = pyyawt.wkeep(a.T,27)
        b5 = pyyawt.wkeep(a.T,27,'c')
        b6 = pyyawt.wkeep(a.T,27,'l')
        b7 = pyyawt.wkeep(a.T,27,'r')
        b8 = pyyawt.wkeep(a,27,2)
        b9 = pyyawt.wkeep(a.T,27,2)
        c = a[:,11:38]
        c1 = a[:,0:27]
        c2 = a[:,-27:]
        c3 = a[:,1:28]
        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b1, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)
        np.testing.assert_almost_equal(b4, b1.T)
        np.testing.assert_almost_equal(b5, b1.T)
        np.testing.assert_almost_equal(b6, b2.T)
        np.testing.assert_almost_equal(b7, b3.T)
        np.testing.assert_almost_equal(b8, c3)
        np.testing.assert_almost_equal(b9, c3.T)

    def test_vector_odd(self):
        a = np.random.randn(51)
        b = pyyawt.wkeep(a,26)
        b1 = pyyawt.wkeep(a,26,'c')
        b2 = pyyawt.wkeep(a,26,'l')
        b3 = pyyawt.wkeep(a,26,'r')
        b4 = pyyawt.wkeep(a,26,2)
        b5 = pyyawt.wkeep(a,26)
        b6 = pyyawt.wkeep(a,26,'c')
        b7 = pyyawt.wkeep(a,26,'l')
        b8 = pyyawt.wkeep(a,26,'r')
        b9 = pyyawt.wkeep(a,26,2)
        c = a[12:38]
        c1 = a[0:26]
        c2 = a[-26:]
        c3 = a[1:27]
        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)
        np.testing.assert_almost_equal(b4, c3)
        np.testing.assert_almost_equal(b5, b)
        np.testing.assert_almost_equal(b6, b)
        np.testing.assert_almost_equal(b7, b2)
        np.testing.assert_almost_equal(b8, b3)
        np.testing.assert_almost_equal(b9, c3)

        a = np.random.randn(1,51)
        b = pyyawt.wkeep(a,26)
        b1 = pyyawt.wkeep(a,26,'c')
        b2 = pyyawt.wkeep(a,26,'l')
        b3 = pyyawt.wkeep(a,26,'r')
        b4 = pyyawt.wkeep(a,26,2)
        b5 = pyyawt.wkeep(a.T,26)
        b6 = pyyawt.wkeep(a.T,26,'c')
        b7 = pyyawt.wkeep(a.T,26,'l')
        b8 = pyyawt.wkeep(a.T,26,'r')
        b9 = pyyawt.wkeep(a.T,26,2)
        c = a[:,12:38]
        c1 = a[:,0:26]
        c2 = a[:,-26:]
        c3 = a[:,1:27]
        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)
        np.testing.assert_almost_equal(b4, c3)
        np.testing.assert_almost_equal(b5, b.T)
        np.testing.assert_almost_equal(b6, b.T)
        np.testing.assert_almost_equal(b7, b2.T)
        np.testing.assert_almost_equal(b8, b3.T)
        np.testing.assert_almost_equal(b9, c3.T)

    def test_matrix(self):
        a = np.random.randn(50,51)
        b = pyyawt.wkeep(a,[27,26])
        b1 = pyyawt.wkeep(a,[27,26],[2,3])
        c = a[11:38,:][:,12:38]
        c1 = a[1:28,:][:,2:28]
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b1, c1)
