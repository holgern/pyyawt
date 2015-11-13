# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt
from __future__ import division, print_function, absolute_import

import os
import numpy as np

import sys
import unittest
import pyyawt
import matplotlib.pyplot as pyplot

snr = 3

# generate test signal
[xref,x] = pyyawt.wnoise(3,11,snr)

level = 5
xd,cxd,lxd = pyyawt.wden(x,'heursure','s','one',level,'sym8')

fig = pyplot.figure()
ax1 = fig.add_subplot(611)
pyplot.plot(xref)
pyplot.title('Original signal')
ax2 = fig.add_subplot(612)
pyplot.plot(x)
pyplot.title('Noisy signal - Signal to noise ratio = '+str(snr))
ax1 = fig.add_subplot(613)
pyplot.plot(xd)
pyplot.title('De-noised signal - heuristic SURE')

xd,cxd,lxd = pyyawt.wden(x,'heursure','s','one',level,'sym8')

ax1 = fig.add_subplot(614)
pyplot.plot(xd)
pyplot.title('De-noised signal - SURE')

xd,cxd,lxd = pyyawt.wden(x,'sqtwolog','s','sln',level,'sym8')

ax1 = fig.add_subplot(615)
pyplot.plot(xd)
pyplot.title('De-noised signal - Fixed form threshold')

xd,cxd,lxd = pyyawt.wden(x,'minimaxi','s','sln',level,'sym8')

ax1 = fig.add_subplot(616)
pyplot.plot(xd)
pyplot.title('De-noised signal - Minimax')

pyplot.show()
