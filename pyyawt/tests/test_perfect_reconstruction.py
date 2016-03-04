# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt


def check_reconstruction(pmode, wavelet, dtype, epsilon=None):
    data_size = [500, 1000, 2000, 10000, 50000, 100000]
    np.random.seed(12345)
    # TODO: smoke testing - more failures for different seeds
    if epsilon is None:
        if dtype == np.float32:
            epsilon = 3e-7
        else:
            # FIXME: limit was 5e-11, but gave failures.  Investigate
            epsilon = 5e-11

    for N in data_size:
        data = np.asarray(np.random.random(N), dtype)

        # compute dwt coefficients
        pa, pd = pyyawt.dwt(data, wavelet, 'mode', pmode)

        # compute reconstruction
        rec = pyyawt.idwt(pa, pd, wavelet,'mode', pmode)

        if len(data) % 2:
            rec = rec[:len(data)]

        rms_rec = np.sqrt(np.mean((data-rec)**2))
        msg = ('[RMS_REC > EPSILON] for Mode: %s, Wavelet: %s, '
               'Length: %d, rms=%.3g' % (pmode, wavelet, len(data), rms_rec))
        np.testing.assert_equal(rms_rec < epsilon,True, err_msg=msg)


class TestPerfectReconstruction(unittest.TestCase):
    def setUp(self):
        self.dwtModes = ['symw', 'asymh', 'asymw', 'zpd', 'sp0', 'sp1','ppd']
        self.biornames = ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6',
                       'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',
                       'bior4.4', 'bior5.5', 'bior6.8']
        self.rbiornames = ['rbior1.1', 'rbior1.3', 'rbior1.5', 'rbior2.2', 'rbior2.4', 'rbior2.6',
                       'rbior2.8', 'rbior3.1', 'rbior3.3', 'rbior3.5', 'rbior3.7', 'rbior3.9',
                       'rbior4.4', 'rbior5.5', 'rbior6.8']

    def test_db(self):
        for Nmodes in np.arange(len(self.dwtModes)-1):
            for N in np.arange(1,37):
                check_reconstruction(self.dwtModes[Nmodes], "db" + str(N), np.float64)

    def test_coif(self):
        for Nmodes in np.arange(len(self.dwtModes)-1):
            for N in np.arange(1,18):
                check_reconstruction(self.dwtModes[Nmodes], "coif" + str(N), np.float64)

    def test_bior(self):
        for Nmodes in np.arange(len(self.dwtModes)-1):
            for N in np.arange(len(self.biornames)):
                check_reconstruction(self.dwtModes[Nmodes], self.biornames[N], np.float64)

    def test_rbior(self):
        for Nmodes in np.arange(len(self.dwtModes)-1):
            for N in np.arange(len(self.rbiornames)):
                check_reconstruction(self.dwtModes[Nmodes], self.rbiornames[N], np.float64)

    def test_sym(self):
        for Nmodes in np.arange(len(self.dwtModes)-1):
            for N in np.arange(2,20):
                check_reconstruction(self.dwtModes[Nmodes], "sym" + str(N), np.float64)
