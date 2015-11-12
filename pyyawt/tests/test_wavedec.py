# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestWavedec(unittest.TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        test_data_file = os.path.join(data_dir, 'Data.npz')
        testFile = np.load(test_data_file)
        self.x1 = testFile["x1"]
        self.x2 = testFile["x2"]
        self.s1 = testFile["s1"]

    def test_haar(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        [cA1,cD1] = pyyawt.dwt(s1,'haar')
        [cA2,cD2] = pyyawt.dwt(cA1,'haar')
        [cA3,cD3] = pyyawt.dwt(cA2,'haar')
        c0 = np.concatenate([cA3, cD3, cD2, cD1])
        l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
        c, l = pyyawt.wavedec(s1,3,'haar')
        np.testing.assert_almost_equal(c, c0)
        np.testing.assert_almost_equal(l, l0)

    def test_db(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,37):
            wname = "db" + str(N)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,wname)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

    def test_coif(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,18):
            wname = "coif" + str(N)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,wname)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

    def test_sym(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(2,21):
            wname = "sym" + str(N)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,wname)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

    def test_bior(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']
        for N in np.arange(len(wnames)):
            wname = wnames[N]
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,wname)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

    def test_type2(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        Lo_D = np.random.randn(20)
        Hi_D = np.random.randn(20)
        [cA1,cD1] = pyyawt.dwt(s1,Lo_D,Hi_D)
        [cA2,cD2] = pyyawt.dwt(cA1,Lo_D,Hi_D)
        [cA3,cD3] = pyyawt.dwt(cA2,Lo_D,Hi_D)
        c0 = np.concatenate([cA3, cD3, cD2, cD1])
        l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
        c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
        np.testing.assert_almost_equal(c, c0)
        np.testing.assert_almost_equal(l, l0)

    def test_bior2(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        ST = pyyawt.dwtmode("status","nodisp")
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']
        for N in np.arange(len(wnames)):
            wname = wnames[N]
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("asymh")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','asymh')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','asymh')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','asymh')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("symw")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','symw')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','symw')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','symw')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("asymw")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','asymw')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','asymw')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','asymw')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("zpd")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','zpd')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','zpd')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','zpd')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("sp0")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','sp0')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','sp0')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','sp0')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("sp1")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','sp1')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','sp1')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','sp1')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("ppd")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','ppd')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','ppd')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','ppd')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode("per")
            [cA1,cD1] = pyyawt.dwt(s1,wname,'mode','per')
            [cA2,cD2] = pyyawt.dwt(cA1,wname,'mode','per')
            [cA3,cD3] = pyyawt.dwt(cA2,wname,'mode','per')
            Lo_D,Hi_D = pyyawt.wfilters(wname,'d')
            c0 = np.concatenate([cA3, cD3, cD2, cD1])
            l0 = [np.size(cA3), np.size(cD3), np.size(cD2), np.size(cD1), np.size(s1)]
            c, l = pyyawt.wavedec(s1,3,Lo_D,Hi_D)
            np.testing.assert_almost_equal(c, c0)
            np.testing.assert_almost_equal(l, l0)

            pyyawt.dwtmode(ST)
