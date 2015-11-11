# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt
import hdf5


class TestDwt(unittest.TestCase):
    def setUp(self):
        demoFile = h5py.File("Data.bin", 'r')
        self.x1 = np.array(demoFile["x1"])
        self.x2 = np.array(demoFile["x2"])
        self.s1 = np.array(demoFile["s1"])

    def test_dwt(self):
        np.testing.assert_almost_equal(1,1)
