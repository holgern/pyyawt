# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestConv(unittest.TestCase):
    def test1(self):
        a = np.random.random((51,1))
        b = np.random.randn(100,1)
        c = pyyawt.conv(a,b)
        c1 = pyyawt.conv(a.T,b).T
        c2 = pyyawt.conv(a,b.T)
        c3 = pyyawt.conv(a.T,b.T).T
        expected = np.convolve(a.flatten(),b.flatten())
        np.testing.assert_almost_equal(c.flatten(), expected)
        np.testing.assert_almost_equal(c1.flatten(), expected)
        np.testing.assert_almost_equal(c2.flatten(), expected)
        np.testing.assert_almost_equal(c3.flatten(), expected)

if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
