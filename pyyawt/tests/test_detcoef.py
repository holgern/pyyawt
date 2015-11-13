# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestDetcoef(unittest.TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        test_data_file = os.path.join(data_dir, 'Data.npz')
        testFile = np.load(test_data_file)
        self.x1 = testFile["x1"]
        self.x2 = testFile["x2"]
        self.s1 = testFile["s1"]

    def test1(self):
        s1 = self.s1
        x1 = self.x1
        x2 = self.x2
        level = 10
        [c,l] = pyyawt.wavedec(s1,level,'sym10')
        cA = []
        cD = []
        [ca,cd] = pyyawt.dwt(s1,'sym10')
        cA.append(ca)
        cD.append(cd)
        for i in np.arange(1,level):
            [ca,cd] = pyyawt.dwt(cA[i-1],'sym10')
            cA.append(ca)
            cD.append(cd)
        cddetMax = pyyawt.detcoef(c,l)
        cdet = []
        for i in np.arange(level):
            cdet.append(pyyawt.detcoef(c,l,i+1))
        np.testing.assert_almost_equal(cddetMax, cD[level-1])
        for i in np.arange(level):
            np.testing.assert_almost_equal(cdet[i], cD[i])
