# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


class TestWextend(unittest.TestCase):
    def test_vector(self):
        a = np.random.randn(1,51)
        b = pyyawt.wextend(1,'symh',a,7)
        b1 = pyyawt.wextend(1,'symh',a,7,'b')
        b2 = pyyawt.wextend(1,'symh',a,7,'l')
        b3 = pyyawt.wextend(1,'symh',a,7,'r')
        c = np.hstack((a[:,6::-1], a, a[:,-1:-8:-1]))
        c1 = np.hstack((a[:,6::-1], a))
        c2 = np.hstack((a, a[:,-1:-8:-1]))
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b1, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)

        b = pyyawt.wextend(1,'sym',a,7)
        b1 = pyyawt.wextend(1,'sym',a,7,'b')
        b2 = pyyawt.wextend(1,'sym',a,7,'l')
        b3 = pyyawt.wextend(1,'sym',a,7,'r')
        c = np.hstack((a[:,6::-1], a, a[:,-1:-8:-1]))
        c1 = np.hstack((a[:,6::-1], a))
        c2 = np.hstack((a, a[:,-1:-8:-1]))
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b1, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)

        b = pyyawt.wextend(1,'symw',a,7)
        b1 = pyyawt.wextend(1,'symw',a,7,'b')
        b2 = pyyawt.wextend(1,'symw',a,7,'l')
        b3 = pyyawt.wextend(1,'symw',a,7,'r')
        c = np.hstack((a[:,7:0:-1], a, a[:,-2:-9:-1]))
        c1 = np.hstack((a[:,7:0:-1], a))
        c2 = np.hstack((a, a[:,-2:-9:-1]))
        np.testing.assert_almost_equal(b, c)
        np.testing.assert_almost_equal(b1, c)
        np.testing.assert_almost_equal(b2, c1)
        np.testing.assert_almost_equal(b3, c2)
