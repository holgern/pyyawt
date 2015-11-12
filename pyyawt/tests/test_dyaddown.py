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
