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
        a = np.random.random((1,5))
        b = pyyawt.qmf(a)
        b1 = pyyawt.qmf(a,0)
        b2 = pyyawt.qmf(a,1)
        b3 = pyyawt.qmf(a.T)
        b4 = pyyawt.qmf(a.T,0)
        b5 = pyyawt.qmf(a.T,1)
        c = np.array([[0, 0, 0, 0, 1],[0, 0, 0, -1, 0],[0, 0, 1, 0, 0],[0, -1, 0, 0, 0],[1, 0, 0, 0, 0]])
        c1 = np.array([[0, 0, 0, 0, -1],[0, 0, 0, 1, 0],[0, 0, -1, 0, 0],[0, 1, 0, 0, 0],[-1,0, 0, 0, 0]])
        d = np.dot(a,c)
        d1 = np.dot(a,c1)
        np.testing.assert_almost_equal(b, b1)
        np.testing.assert_almost_equal(b, d)
        np.testing.assert_almost_equal(b2, d1)
        np.testing.assert_almost_equal(b3, b.T)
        np.testing.assert_almost_equal(b4, b.T)
        np.testing.assert_almost_equal(b5, d1.T)

if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
