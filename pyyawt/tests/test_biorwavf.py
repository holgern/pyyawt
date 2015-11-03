# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestBiorwavf(unittest.TestCase):
    def setUp(self):
        self.wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']

    def test_sum(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.biorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF),np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF),np.sqrt(2))

    def test_sumEven(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.biorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF[::2]), 1./np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF[::2]), 1./np.sqrt(2))

    def test_sumOdd(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.biorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF[1::2]), 1./np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF[1::2]), 1./np.sqrt(2))


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
