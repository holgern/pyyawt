# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestWmaxlev(unittest.TestCase):
    def test1(self):
        l1 = pyyawt.wmaxlev(50,'db4')
        l2 = pyyawt.wmaxlev(65,'db4')
        l3 = pyyawt.wmaxlev([50, 65],'db4')
        l4 = pyyawt.wmaxlev([65, 50],'db4')
        l5 = pyyawt.wmaxlev([50, 65],'db4')
        l6 = pyyawt.wmaxlev([65, 50],'db4')
        l = np.min([l1,l2])

        np.testing.assert_almost_equal(l3, l)
        np.testing.assert_almost_equal(l4, l)
        np.testing.assert_almost_equal(l5, l)
        np.testing.assert_almost_equal(l6, l)

if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
