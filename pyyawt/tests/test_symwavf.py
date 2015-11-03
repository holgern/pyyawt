# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestSymwavf(unittest.TestCase):
    def test_sum(self):
        for N in np.arange(2,20):
            w = pyyawt.symwavf("sym" + str(N))
            np.testing.assert_almost_equal(np.sum(w)-np.sqrt(2),0)

    def test_sumEven(self):
        for N in np.arange(2,20):
            w = pyyawt.symwavf("sym" + str(N))
            np.testing.assert_almost_equal(np.sum(w[::2]), 1./np.sqrt(2))

    def test_sumOdd(self):
        for N in np.arange(2,20):
            w = pyyawt.symwavf("sym" + str(N))
            np.testing.assert_almost_equal(np.sum(w[1::2]),1./np.sqrt(2))

    def test_ZeroOne(self):
        for N in np.arange(2,20):
            m = 0
            w = pyyawt.symwavf("sym" + str(N))
            np.testing.assert_almost_equal(np.sum(w[2*m:(2*N+2*m)]*w[0:2*N]),1)
            for m in np.arange(1, N-1):
                np.testing.assert_almost_equal(np.sum(w[2*m:(2*N)]*w[0:2*N-2*m]),0)


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
