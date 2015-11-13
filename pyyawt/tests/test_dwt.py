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
        [cA,cD] = pyyawt.dwt(x1,'haar')
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters('haar')
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)

        [cA,cD] = pyyawt.dwt(x2,'haar')
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)

        [cA,cD] = pyyawt.dwt(s1,'haar')
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)

    def test_db(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,31):
            wname = "db" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(x2,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(s1,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

    def test_coif(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(1,11):
            wname = "coif" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(x2,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(s1,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

    def test_sym(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        for N in np.arange(4,11):
            wname = "sym" + str(N)
            [cA,cD] = pyyawt.dwt(x1,wname)
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(x2,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(s1,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

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
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(x2,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(s1,wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

    def test_dwtMode(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        # ST = pyyawt.dwtmode("status","nodisp")
        dwtModes = ['symw', 'asymh', 'asymw', 'zpd', 'sp0', 'sp1',
                  'ppd','per']
        wname = "db10"
        for N in np.arange(len(dwtModes)-1):
            # pyyawt.dwtmode(dwtModes[N])
            [cA,cD] = pyyawt.dwt(x1,wname,"mode",dwtModes[N])
            [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(x2,wname,"mode",dwtModes[N])
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

            [cA,cD] = pyyawt.dwt(s1,wname,"mode",dwtModes[N])
            caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1)+np.size(Lo_D)-1)))
            cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1)+np.size(Lo_D)-1)))
            np.testing.assert_almost_equal(caa.flatten(), cA)
            np.testing.assert_almost_equal(cdd.flatten(), cD)

        N = len(dwtModes) - 1
        [cA,cD] = pyyawt.dwt(x1,wname,"mode",dwtModes[N])
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters(wname)
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x1,np.size(Lo_D),'b'),Lo_D),(np.size(x1))))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x1,np.size(Lo_D),'b'),Hi_D),(np.size(x1))))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)

        [cA,cD] = pyyawt.dwt(x2,wname,"mode",dwtModes[N])
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x2,np.size(Lo_D),'b'),Lo_D),(np.size(x2))))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],x2,np.size(Lo_D),'b'),Hi_D),(np.size(x2))))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)

        [cA,cD] = pyyawt.dwt(s1,wname,"mode",dwtModes[N])
        caa = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],s1,np.size(Lo_D),'b'),Lo_D),(np.size(s1))))
        cdd = pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,dwtModes[N],s1,np.size(Lo_D),'b'),Hi_D),(np.size(s1))))
        np.testing.assert_almost_equal(caa.flatten(), cA)
        np.testing.assert_almost_equal(cdd.flatten(), cD)
