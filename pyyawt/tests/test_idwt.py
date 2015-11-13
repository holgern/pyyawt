# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDwt(unittest.TestCase):
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
        wname = 'haar'
        [cA,cD] = pyyawt.dwt(x1,wname)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,wname)
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(x1,wname)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,wname)
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(x1,wname)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,wname)
        np.testing.assert_almost_equal(r, x0)

    def test_db(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,11):
            wname = "db" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

    def test_coif(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,6):
            wname = "coif" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

    def test_sym(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(4,11):
            wname = "sym" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

    def test_bior(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wnames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                  'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                  'bior4.4', 'bior5.5', 'bior6.8']
        for N in np.arange(1,13):
            wname = wnames[N]
            [cA,cD] = pyyawt.dwt(x1,wname)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
            d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
            x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
            r = pyyawt.idwt(cA,cD,wname)
            np.testing.assert_almost_equal(r, x0)

    def test_type2(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wname = "bior3.9"
        [cA,cD] = pyyawt.dwt(x1,wname)
        Lo_R = np.random.randn(20)
        Hi_R = np.random.randn(20)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,Lo_R.flatten(),Hi_R.flatten())
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(x2,wname)
        Lo_R = np.random.randn(20)
        Hi_R = np.random.randn(20)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,Lo_R.flatten(),Hi_R.flatten())
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(s1,wname)
        Lo_R = np.random.randn(20)
        Hi_R = np.random.randn(20)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        r = pyyawt.idwt(cA,cD,Lo_R.flatten(),Hi_R.flatten())
        np.testing.assert_almost_equal(r, x0)

    def test_type3(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wname = "sym8"
        [cA,cD] = pyyawt.dwt(x1,wname)
        r = pyyawt.idwt(cA,cD,wname,50)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,50)
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(x2,wname)
        r = pyyawt.idwt(cA,cD,wname,50)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,50)
        np.testing.assert_almost_equal(r, x0)

        [cA,cD] = pyyawt.dwt(s1,wname)
        r = pyyawt.idwt(cA,cD,wname,50)
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,50)
        np.testing.assert_almost_equal(r, x0)

        Lo_R = np.random.randn(50)
        Hi_R = np.random.randn(50)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,50)
        r = pyyawt.idwt(cA,cD,Lo_R.flatten(),Hi_R.flatten(),50)
        np.testing.assert_almost_equal(r, x0)

    def test_type4(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wname = "db7"
        dwtModes = ['symw', 'asymh', 'asymw', 'zpd', 'sp0', 'sp1',
                    'ppd','per']
        [cA,cD] = pyyawt.dwt(x1,wname)
        [Lo_R,Hi_R] = pyyawt.wfilters(wname,'r')
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA)-np.size(Lo_R)+2)
        for N in np.arange(len(dwtModes)-1):
            r = pyyawt.idwt(cA,cD,wname,'mode',dwtModes[N])
            np.testing.assert_almost_equal(r, x0)
        r = pyyawt.idwt(cA,cD,wname,'mode','per')
        x0 = pyyawt.wkeep(a0+d0,2*np.size(cA))
        np.testing.assert_almost_equal(r, x0)

    def test_type5(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wname = "db7"
        dwtModes = ['symw', 'asymh', 'asymw', 'zpd', 'sp0', 'sp1',
                    'ppd','per']
        [cA,cD] = pyyawt.dwt(x1,wname)
        [Lo_R,Hi_R] = pyyawt.wfilters(wname,'r')
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R)
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R)
        x0 = pyyawt.wkeep(a0+d0,50)
        for N in np.arange(len(dwtModes)-1):
            r = pyyawt.idwt(cA,cD,wname,50,'mode',dwtModes[N])
            np.testing.assert_almost_equal(r, x0)
        r = pyyawt.idwt(cA,cD,wname,50,'mode','per')
        x0 = pyyawt.wkeep(a0+d0,50)
        np.testing.assert_almost_equal(r, x0)

    def test_type6(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        wname = "db7"
        dwtModes = ['symw', 'asymh', 'asymw', 'zpd', 'sp0', 'sp1',
                    'ppd','per']
        [cA,cD] = pyyawt.dwt(x1,wname)
        Lo_R = np.random.rand(14)
        Hi_R = np.random.rand(14)
        a0 = pyyawt.conv(pyyawt.dyadup(cA),Lo_R.flatten())
        d0 = pyyawt.conv(pyyawt.dyadup(cD),Hi_R.flatten())
        x0 = pyyawt.wkeep(a0+d0,50)
        for N in np.arange(len(dwtModes)-1):
            r = pyyawt.idwt(cA,cD,Lo_R,Hi_R,50,'mode',dwtModes[N])
            np.testing.assert_almost_equal(r, x0)
        r = pyyawt.idwt(cA,cD,Lo_R,Hi_R,50,'mode','per')
        x0 = pyyawt.wkeep(a0+d0,50)
        np.testing.assert_almost_equal(r, x0)
