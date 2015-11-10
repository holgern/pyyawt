# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestWfilters(unittest.TestCase):
    def test_haar(self):
        Lo_D,Hi_D,Lo_R,Hi_R = pyyawt.wfilters('haar')
        w = pyyawt.dbwavf('db1')
        lo_d,hi_d,lo_r,hi_r = pyyawt.orthfilt(w)
        np.testing.assert_almost_equal(lo_d, Lo_D)
        np.testing.assert_almost_equal(lo_r, Lo_R)
        np.testing.assert_almost_equal(hi_d, Hi_D)
        np.testing.assert_almost_equal(hi_r, Hi_R)

        Lo_D,Hi_D = pyyawt.wfilters('haar','d')
        Lo_R,Hi_R = pyyawt.wfilters('haar','r')
        np.testing.assert_almost_equal(lo_d, Lo_D)
        np.testing.assert_almost_equal(lo_r, Lo_R)
        np.testing.assert_almost_equal(hi_d, Hi_D)
        np.testing.assert_almost_equal(hi_r, Hi_R)

        Lo_D,Lo_R = pyyawt.wfilters('haar','l')
        Hi_D,Hi_R = pyyawt.wfilters('haar','h')
        np.testing.assert_almost_equal(lo_d, Lo_D)
        np.testing.assert_almost_equal(lo_r, Lo_R)
        np.testing.assert_almost_equal(hi_d, Hi_D)
        np.testing.assert_almost_equal(hi_r, Hi_R)

    def test_dbwavf(self):
        for N in np.arange(1,37):
            Lo_D,Hi_D,Lo_R,Hi_R = pyyawt.wfilters("db" + str(N))
            w = pyyawt.dbwavf("db" + str(N))
            lo_d,hi_d,lo_r,hi_r = pyyawt.orthfilt(w)
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

            Lo_D,Hi_D = pyyawt.wfilters("db" + str(N),'d')
            Lo_R,Hi_R = pyyawt.wfilters("db" + str(N),'r')
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

            Lo_D,Lo_R = pyyawt.wfilters("db" + str(N),'l')
            Hi_D,Hi_R = pyyawt.wfilters("db" + str(N),'h')
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

    def test_biorwavf(self):
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                       'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                       'bior4.4', 'bior5.5', 'bior6.8']
        for N in np.arange(len(wnames)):
            Lo_D,Hi_D,Lo_R,Hi_R = pyyawt.wfilters(wnames[N])
            rf,df = pyyawt.biorwavf(wnames[N])
            lo_d,hi_d,lo_r,hi_r = pyyawt.biorfilt(df,rf)
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

            Lo_D,Hi_D = pyyawt.wfilters(wnames[N],'d')
            Lo_R,Hi_R = pyyawt.wfilters(wnames[N],'r')
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

            Lo_D,Lo_R = pyyawt.wfilters(wnames[N],'l')
            Hi_D,Hi_R = pyyawt.wfilters(wnames[N],'h')
            np.testing.assert_almost_equal(lo_d, Lo_D)
            np.testing.assert_almost_equal(lo_r, Lo_R)
            np.testing.assert_almost_equal(hi_d, Hi_D)
            np.testing.assert_almost_equal(hi_r, Hi_R)

if __name__ == '__main__':
    # run_module_suite(argv=sys.argv)
    unittest.main()
