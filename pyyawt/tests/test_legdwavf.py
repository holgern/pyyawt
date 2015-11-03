# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestLegdwavf(unittest.TestCase):
    def test_sum(self):
        for N in np.arange(1,10):
            w = pyyawt.legdwavf("legd" + str(N))
            np.testing.assert_almost_equal(np.sum(w)-np.sqrt(2),0,decimal=5)

    def test_sumEven(self):
        for N in np.arange(1,10):
            w = pyyawt.legdwavf("legd" + str(N))
            np.testing.assert_almost_equal(np.sum(w[::2]), 1./np.sqrt(2),decimal=5)

    def test_sumOdd(self):
        for N in np.arange(1,10):
            w = pyyawt.legdwavf("legd" + str(N))
            np.testing.assert_almost_equal(np.sum(w[1::2]),1./np.sqrt(2),decimal=5)


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
