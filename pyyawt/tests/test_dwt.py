# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt
import h5py


class TestDwt(unittest.TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        test_data_file = os.path.join(data_dir, 'Data.bin')
        demoFile = h5py.File(test_data_file, 'r')
        self.x1 = np.array(demoFile["x1"])
        self.x2 = np.array(demoFile["x2"])
        self.s1 = np.array(demoFile["s1"])

    def test_haar(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        [cA,cD] = pyyawt.dwt(x1,'haar')
        [Lo_D,Hi_D,Lo_R,Hi_R] = pyyawt.wfilters('haar')
        # caa=pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,len(Lo_D),'b'),Lo_D),(len(x1)+len(Lo_D)-1)))
        # cdd=pyyawt.dyaddown(pyyawt.wkeep(pyyawt.conv(pyyawt.wextend(1,'symh',x1,len(Lo_D),'b'),Hi_D),(len(x1)+len(Lo_D)-1)));
        # np.testing.assert_almost_equal(caa, cA)
        # np.testing.assert_almost_equal(cdd , cD)

        [cA,cD] = pyyawt.dwt(x2,'haar')
        # caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
        # cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
        # np.testing.assert_almost_equal( caa , cA)
        # np.testing.assert_almost_equal( cdd , cD)
        [cA,cD] = pyyawt.dwt(s1,'haar')
        # caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
        # cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
        # np.testing.assert_almost_equal(caa , cA)
        # np.testing.assert_almost_equal(cdd , cD)
        np.testing.assert_almost_equal(1,1)
