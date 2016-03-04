# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestCoifwavf(unittest.TestCase):
    def setUp(self):
        self.accuracy_decimal = 15

    def test_sum(self):
        for N in np.arange(1,18):
            w = pyyawt.coifwavf("coif" + str(N))
            np.testing.assert_almost_equal(np.sum(w)-np.sqrt(2),0, decimal=self.accuracy_decimal)

    def test_sumEven(self):
        for N in np.arange(1,18):
            w = pyyawt.coifwavf("coif" + str(N))
            np.testing.assert_almost_equal(np.sum(w[::2]), 1./np.sqrt(2), decimal=self.accuracy_decimal)

    def test_sumOdd(self):
        for N in np.arange(1,18):
            w = pyyawt.coifwavf("coif" + str(N))
            np.testing.assert_almost_equal(np.sum(w[1::2]),1./np.sqrt(2), decimal=self.accuracy_decimal)


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
