# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestRBiorwavf(unittest.TestCase):
    def setUp(self):
        self.wnames = ['rbior1.1', 'rbior1.3', 'rbior1.5', 'rbior2.2', 'rbior2.4', 'rbior2.6',
                  'rbior2.8', 'rbior3.1', 'rbior3.3', 'rbior3.5', 'rbior3.7', 'rbior3.9',
                  'rbior4.4', 'rbior5.5', 'rbior6.8']

    def test_sum(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.rbiorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF),np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF),np.sqrt(2))

    def test_sumEven(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.rbiorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF[::2]), 1./np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF[::2]), 1./np.sqrt(2))

    def test_sumOdd(self):
        for N in np.arange(len(self.wnames)):
            RF,DF = pyyawt.rbiorwavf(self.wnames[N])
            np.testing.assert_almost_equal(np.sum(RF[1::2]), 1./np.sqrt(2))
            np.testing.assert_almost_equal(np.sum(DF[1::2]), 1./np.sqrt(2))


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
