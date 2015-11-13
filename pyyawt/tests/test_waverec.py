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
        [c,l] = pyyawt.wavedec(s1,3,'haar')
        [cA1,cD1] = pyyawt.dwt(s1,'haar')
        [cA2,cD2] = pyyawt.dwt(cA1,'haar')
        [cA3,cD3] = pyyawt.dwt(cA2,'haar')
        ca2 = pyyawt.idwt(cA3,cD3,'haar',np.size(cA2))
        ca1 = pyyawt.idwt(ca2,cD2,'haar',np.size(cA1))
        a0 = pyyawt.idwt(ca1,cD1,'haar',np.size(s1))
        s0 = pyyawt.waverec(c,l,'haar')
        np.testing.assert_almost_equal(a0, s0)

    def test_db(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,37):
            wname = "db" + str(N)
            [c,l] = pyyawt.wavedec(s1,3,wname)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            ca2 = pyyawt.idwt(cA3,cD3,wname,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,wname,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,wname,np.size(s1))
            s0 = pyyawt.waverec(c,l,wname)
            np.testing.assert_almost_equal(a0, s0)

    def test_coif(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,18):
            wname = "coif" + str(N)
            [c,l] = pyyawt.wavedec(s1,3,wname)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            ca2 = pyyawt.idwt(cA3,cD3,wname,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,wname,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,wname,np.size(s1))
            s0 = pyyawt.waverec(c,l,wname)
            np.testing.assert_almost_equal(a0, s0)

    def test_sym(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(2,21):
            wname = "sym" + str(N)
            [c,l] = pyyawt.wavedec(s1,3,wname)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            ca2 = pyyawt.idwt(cA3,cD3,wname,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,wname,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,wname,np.size(s1))
            s0 = pyyawt.waverec(c,l,wname)
            np.testing.assert_almost_equal(a0, s0)

    def test_bior(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']
        for N in np.arange(len(wnames)):
            wname = wnames[N]
            [c,l] = pyyawt.wavedec(s1,3,wname)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            ca2 = pyyawt.idwt(cA3,cD3,wname,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,wname,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,wname,np.size(s1))
            s0 = pyyawt.waverec(c,l,wname)
            np.testing.assert_almost_equal(a0, s0)

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
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            [c,l] = pyyawt.wavedec(s1,3,wname)
            s0 = pyyawt.waverec(c,l,wname)
            [cA1,cD1] = pyyawt.dwt(s1,wname)
            [cA2,cD2] = pyyawt.dwt(cA1,wname)
            [cA3,cD3] = pyyawt.dwt(cA2,wname)
            ca2 = pyyawt.idwt(cA3,cD3,wname,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,wname,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,wname,np.size(s1))
            np.testing.assert_almost_equal(a0, s0)

            Lo_D = np.random.randn(np.size(Lo_D))
            Hi_D = np.random.randn(np.size(Lo_D))
            ca2 = pyyawt.idwt(cA3,cD3,Lo_R,Hi_R,np.size(cA2))
            ca1 = pyyawt.idwt(ca2,cD2,Lo_R,Hi_R,np.size(cA1))
            a0 = pyyawt.idwt(ca1,cD1,Lo_R,Hi_R,np.size(s1))
            s0 = pyyawt.waverec(c,l,Lo_R,Hi_R)
            np.testing.assert_almost_equal(a0, s0)

    def test_bior3(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        ST = pyyawt.dwtmode("status","nodisp")
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']
        ST = pyyawt.dwtmode("status","nodisp")
        for N in np.arange(len(wnames)):
            wname = wnames[N]
            [c,l] = pyyawt.wavedec(s1,3,wname)
            pyyawt.dwtmode("symh")
            a0 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("symw")
            a1 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("asymh")
            a2 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("asymw")
            a3 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("zpd")
            a4 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("sp0")
            a5 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("sp1")
            a6 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("ppd")
            a7 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode("per")
            a8 = pyyawt.waverec(c,l,wname)
            pyyawt.dwtmode(ST)
            np.testing.assert_almost_equal(a0, a1)
            np.testing.assert_almost_equal(a0, a2)
            np.testing.assert_almost_equal(a0, a3)
            np.testing.assert_almost_equal(a0, a4)
            np.testing.assert_almost_equal(a0, a5)
            np.testing.assert_almost_equal(a0, a6)
            np.testing.assert_almost_equal(a0, a7)
            np.testing.assert_almost_equal(a0, a8)
