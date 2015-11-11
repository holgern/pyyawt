# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDwtMode(unittest.TestCase):
    def test1(self):

        a = 'asymh'
        pyyawt.dwtmode(a)
        ST = pyyawt.dwtmode('status','nodisp')
        b = (ST == a)

        x = np.random.randn(50)
        [cA,cD] = pyyawt.dwt(x,'db2')
        [caa,cdd] = pyyawt.dwt(x,'db2','mode',a)
        # x0 = pyyawt.idwt(cA,cD,'db2',len(x))

        np.testing.assert_almost_equal(cA, caa)
        np.testing.assert_almost_equal(cD, cdd)
        # np.testing.assert_almost_equal(x, x0)

        np.testing.assert_equal(b,True)

        a = 'sp1'
        pyyawt.dwtmode(a)
        ST = pyyawt.dwtmode('status','nodisp')
        b = (ST == a)
        x = np.random.randn(50)
        [cA,cD] = pyyawt.dwt(x,'db2')
        [caa,cdd] = pyyawt.dwt(x,'db2','mode',a)
        # x0=pyyawt.idwt(cA,cD,'db2',length(x))

        np.testing.assert_almost_equal(cA, caa)
        np.testing.assert_almost_equal(cD, cdd)
        # np.testing.assert_almost_equal(x, x0)

        np.testing.assert_equal(b,True)

        pyyawt.dwtmode("symh")


if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
