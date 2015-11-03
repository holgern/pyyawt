# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDWT(unittest.TestCase):
    def test_wavelet_parser(self):
        family = -1
        member = -1
        wname = "db2"
        family,member = pyyawt._wavelet_parser(wname.encode())
        np.testing.assert_equal(family, 1)
        np.testing.assert_equal(member, 2)

        wname = "coif2"
        family,member = pyyawt._wavelet_parser(wname.encode())
        np.testing.assert_equal(family, 2)

        wname = "sym7"
        family,member = pyyawt._wavelet_parser(wname.encode())
        np.testing.assert_equal(family, 3)

    def test_dbwavf(self):
        F = pyyawt.dbwavf("db2")
        F_req = np.array([0.4829629, 0.8365163, 0.2241439,-0.1294095])

        np.testing.assert_almost_equal(F,F_req)

    def test_coifwavf(self):
        F = pyyawt.coifwavf("coif1")
        F_req = np.array([- 0.0727326, 0.3378977, 0.8525720, 0.3848648, -0.0727326, -0.0156557])

        np.testing.assert_almost_equal(F,F_req)


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
