# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt

__doc__ = """Pyrex wrapper for low-level C edflib implementation."""
__all__ = ["PYYAWT_HAAR", "PYYAWT_DAUBECHIES", "PYYAWT_COIFLETS", "PYYAWT_SYMLETS",
        "PYYAWT_SPLINE_BIORTH", "PYYAWT_BEYLKIN", "PYYAWT_VAIDYANATHAN", "PYYAWT_DMEY",
        "PYYAWT_BATHLETS", "PYYAWT_LEGENDRE", "PYYAWT_SPLINE_RBIORTH", "PYYAWT_FARRAS",
        "PYYAWT_KINGSBURYQ", "PYYAWT_NOT_DEFINED", "PYYAWT_ORTH", "PYYAWT_BIORTH",
        "PYYAWT_ZPD", "PYYAWT_SYMH", "PYYAWT_SYMW", "PYYAWT_ASYMH", "PYYAWT_ASYMW",
        "PYYAWT_SP0", "PYYAWT_SP1", "PYYAWT_PPD", "PYYAWT_PER",
        "_wavelet_parser", "_dbwavf","_coifwavf","_symwavf", "_legdwavf", "_biorwavf",  "_rbiorwavf",
        "_wfilters_length", "_conv", "_orthfilt", "_biorfilt","_haarwavf",
        "_beylkinwavf", "_vaidyanathanwavf", "_dmeywavf", "_bathletswavf", "_legendrewavf",
        "_farraswavf", "_kingsburyqwavf", "_wave_len_validate", "_getdwtMode", "_dwtWrite",
        "_dwt_neo","_qmf_odd","_qmf_even", "_wkeep_1D_center", "_wkeep_1D_left", "_wkeep_1D_right",
        "_wkeep_1D_index", "_wkeep_2D_center", "_wkeep_2D_index", "_dyaddown_1D_keep_odd",
        "_dyaddown_1D_keep_even", "_dyaddown_2D_keep_odd_row", "_dyaddown_2D_keep_odd_col", "_dyaddown_2D_keep_even_row",
        "_dyaddown_2D_keep_even_col", "_dyaddown_2D_keep_odd", "_dyaddown_2D_keep_even",
        "_dyadup_1D_feed_odd", "_dyadup_1D_feed_even", "_dyadup_2D_feed_odd_row", "_dyadup_2D_feed_odd_col", "_dyadup_2D_feed_even_row",
        "_dyadup_2D_feed_even_col", "_dyadup_2D_feed_odd", "_dyadup_2D_feed_even",
        "_wextend_1D_center", "_wextend_1D_left", "_wextend_1D_right", "_wextend_2D",
        "_wextend_2D_col", "_wextend_2D_row", "_waverec", "_wavedec", "_wave_dec_len_cal",
        "_idwt_neo", "_idwt_detail_neo", "_idwt_approx_neo", "_detcoef"]


from c_pyyawt cimport *
cimport cpython
import numpy as np
cimport numpy as np
from cpython.version cimport PY_MAJOR_VERSION
include "pyyawt.pxi"
# from libc.stdlib cimport malloc, free

# constants are redeclared here so we can access them from Python
PYYAWT_HAAR = HAAR
PYYAWT_DAUBECHIES = DAUBECHIES
PYYAWT_COIFLETS = COIFLETS
PYYAWT_SYMLETS = SYMLETS
PYYAWT_SPLINE_BIORTH = SPLINE_BIORTH
PYYAWT_BEYLKIN = BEYLKIN
PYYAWT_VAIDYANATHAN = VAIDYANATHAN
PYYAWT_DMEY = DMEY
PYYAWT_BATHLETS = BATHLETS
PYYAWT_LEGENDRE = LEGENDRE
PYYAWT_SPLINE_RBIORTH = SPLINE_RBIORTH
PYYAWT_FARRAS = FARRAS
PYYAWT_KINGSBURYQ = KINGSBURYQ
PYYAWT_NOT_DEFINED = NOT_DEFINED
PYYAWT_ORTH = ORTH
PYYAWT_BIORTH = BIORTH

PYYAWT_ZPD = 0
PYYAWT_SYMH = 1
PYYAWT_SYMW = 2
PYYAWT_ASYMH = 3
PYYAWT_ASYMW = 4
PYYAWT_SP0 = 5
PYYAWT_SP1 = 6
PYYAWT_PPD = 7
PYYAWT_PER = 8


def _wavelet_parser(char *wname):
        cdef int ret_family = -1
        cdef int ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        return ret_family, ret_member

def _wave_len_validate(int sigInLen, int waveLength):
        cdef int lev = -1
        cdef int val = -1
        wave_len_validate(sigInLen, waveLength, &lev, &val)
        return lev, val

def _dbwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                daubechies_synthesis_initialize(ret_member,&thisptr)
        else:
                daubechies_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()

def _coifwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                coiflets_synthesis_initialize(ret_member,&thisptr)
        else:
                coiflets_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)

        filter_clear()
  
def _symwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                symlets_synthesis_initialize(ret_member,&thisptr)
        else:
                symlets_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)

        filter_clear()

def _legdwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                legendre_synthesis_initialize(ret_member,&thisptr)
        else:
                legendre_analysis_initialize(ret_member,&thisptr)
        m2 = 1
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)

        filter_clear()

def _biorwavf(char *wname, char *filtertype, flip, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                sp_bior_synthesis_initialize(ret_member,&thisptr)
        else:
                sp_bior_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (flip):
                if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                        wrev (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
                else:
                        wrev (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
                else:
                        verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()

def _rbiorwavf(char *wname, char *filtertype, flip, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                sp_rbior_synthesis_initialize(ret_member,&thisptr)
        else:
                sp_rbior_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (flip):
                if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                        wrev (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
                else:
                        wrev (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
                else:
                        verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()

def _wfilters_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (ret_family == DAUBECHIES):
                daubechies_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == COIFLETS):
                coiflets_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == SYMLETS):
                symlets_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == SPLINE_BIORTH):
                sp_bior_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == BEYLKIN):
                beylkin_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == VAIDYANATHAN):
                vaidyanathan_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == DMEY):
                dmey_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == LEGENDRE):
                legendre_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == SPLINE_RBIORTH):
                sp_rbior_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == HAAR):
                haar_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == BATHLETS):
                bathlets_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == FARRAS):
                farras_synthesis_initialize(ret_member,&thisptr)
        elif (ret_family == KINGSBURYQ):
                kingsburyq_synthesis_initialize(ret_member,&thisptr)
        filter_clear()
        return thisptr.length

        
def _conv(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] input2, np.ndarray[np.float64_t, ndim=2] output1):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = input2.shape[0]
        n2 = input2.shape[1]
        cdef int m3, n3
        m3 = output1.shape[0]
        n3 = output1.shape[1]

        conv (<double*>input1.data, m1*n1, <double*>output1.data, m3*n3, <double*>input2.data, m2*n2)

def _orthfilt(np.ndarray[np.float64_t, ndim=1] input1, np.ndarray[np.float64_t, ndim=1] output1, np.ndarray[np.float64_t, ndim=1] output2, np.ndarray[np.float64_t, ndim=1] output3, np.ndarray[np.float64_t, ndim=1] output4):
        cdef int m1, n1
        m1 = 1
        n1 = input1.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output1.shape[0]
        cdef int m3, n3
        m3 = 1
        n3 = output2.shape[0]
        cdef int m4, n4
        m4 = 1
        n4 = output3.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output4.shape[0]

        orth_filt_group (<double*>input1.data, m1*n1, <double*>output3.data, <double*>output1.data, <double*>output4.data, <double*>output2.data)

def _biorfilt(np.ndarray[np.float64_t, ndim=1] input1, np.ndarray[np.float64_t, ndim=1] input2, np.ndarray[np.float64_t, ndim=1] output1, np.ndarray[np.float64_t, ndim=1] output2, np.ndarray[np.float64_t, ndim=1] output3, np.ndarray[np.float64_t, ndim=1] output4):
        cdef int m1, n1
        m1 = 1
        n1 = input1.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = input2.shape[0]
        cdef int m3, n3
        m3 = 1
        n3 = output1.shape[0]
        cdef int m4, n4
        m4 = 1
        n4 = output2.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output3.shape[0]
        cdef int m6, n6
        m6 = 1
        n6 = output4.shape[0]


        bior_filt_group (<double*>input1.data, m1*n1, <double*>input2.data, m2*n2, <double*>output1.data, m3*n3, <double*>output2.data, m4*n4, <double*>output3.data, m5*n5, <double*>output4.data, m6*n6)


def _haarwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                haar_synthesis_initialize(ret_member,&thisptr)
        else:
                haar_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()
        

def _beylkinwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                beylkin_synthesis_initialize(ret_member,&thisptr)
        else:
                beylkin_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()


def _vaidyanathanwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                vaidyanathan_synthesis_initialize(ret_member,&thisptr)
        else:
                vaidyanathan_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()



def _dmeywavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                dmey_synthesis_initialize(ret_member,&thisptr)
        else:
                dmey_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()


def _bathletswavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                bathlets_synthesis_initialize(ret_member,&thisptr)
        else:
                bathlets_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()
        

def _legendrewavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                legendre_synthesis_initialize(ret_member,&thisptr)
        else:
                legendre_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()


def _farraswavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                farras_synthesis_initialize(ret_member,&thisptr)
        else:
                farras_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()


def _kingsburyqwavf(char *wname, char *filtertype, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        if (filtertype == b'Lo_R' or filtertype == b'Hi_R'):
                kingsburyq_synthesis_initialize(ret_member,&thisptr)
        else:
                kingsburyq_analysis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        if (filtertype == b'Lo_D' or filtertype == b'Lo_R'):
                verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
        else:
                verbatim_copy (thisptr.pHiPass, m2*n2, <double*>sigbuf.data, m2*n2)
        filter_clear()

def _getdwtMode():
        if (getdwtMode() == ZPD):
                return PYYAWT_ZPD
        elif (getdwtMode() == SYMH):
                return PYYAWT_SYMH
        elif (getdwtMode() == SYMW):
                return PYYAWT_SYMW
        elif (getdwtMode() == ASYMH):
                return PYYAWT_ASYMH
        elif (getdwtMode() == ASYMW):
                return PYYAWT_ASYMW
        elif (getdwtMode() == SP0):
                return PYYAWT_SP0
        elif (getdwtMode() == SP1):
                return PYYAWT_SP1
        elif (getdwtMode() == PPD):
                return PYYAWT_PPD
        elif (getdwtMode() == PER):
                return PYYAWT_PER

def _dwtWrite(char* mode):
        cdef int errCode
        dwt_write(mode, &errCode)
        #cdef int dwtMode
        #if (status == PYYAWT_ZPD):
                #dwtMode = ZPD
        #elif (status == PYYAWT_SYMH):
                #dwtMode = SYMH
        #elif (status == PYYAWT_SYMW):
                #dwtMode = SYMW
        #elif (status == PYYAWT_ASYMH):
                #dwtMode = ASYMH
        #elif (status == PYYAWT_ASYMW):
                #dwtMode = ASYMW
        #elif (status == PYYAWT_SP0):
                #dwtMode = SP0
        #elif (status == PYYAWT_SP1):
                #dwtMode = SP1
        #elif (status == PYYAWT_PPD):
                #dwtMode = PPD
        #elif (status == PYYAWT_PER):
                #dwtMode = PER
        #setdwtMode(dwtMode)
        return errCode
        
def _dwt_neo(np.ndarray[np.float64_t, ndim=1] input1, np.ndarray[np.float64_t, ndim=1] input2, np.ndarray[np.float64_t, ndim=1] input3, np.ndarray[np.float64_t, ndim=1] output1, np.ndarray[np.float64_t, ndim=1] output2):
        cdef int m1, n1
        m1 = 1
        n1 = input1.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = input2.shape[0]
        cdef int m3, n3
        m3 = 1
        n3 = input3.shape[0]
        cdef int m4, n4
        m4 = 1
        n4 = output1.shape[0]

        dwt_neo (<double*>input1.data, m1*n1, <double*>input2.data,  <double*>input3.data,m3*n3,  <double*>output1.data, <double*>output2.data, m4*n4, getdwtMode());
        filter_clear()
        
def _qmf_odd(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        qmf_odd (<double*>input.data, m1*n1, <double*>output.data, m2*n2)

def _qmf_even(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        qmf_even (<double*>input.data, m1*n1, <double*>output.data, m2*n2)

def _wkeep_1D_center(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wkeep_1D_center (<double*>input.data, m1*n1, <double*>output.data, m2*n2)

def _wkeep_1D_left(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wkeep_1D_left (<double*>input.data, m1*n1, <double*>output.data, m2*n2)
        
def _wkeep_1D_right(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wkeep_1D_right (<double*>input.data, m1*n1, <double*>output.data, m2*n2)
        
def _wkeep_1D_index(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output, first):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wkeep_1D_index (<double*>input.data, m1*n1, <double*>output.data, m2*n2, first)

def _wkeep_2D_center(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        wkeep_2D_center (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _wkeep_2D_index(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output, rowFirst, colFirst):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        wkeep_2D_index (<double*>input.data, m1, n1, <double*>output.data, m2, n2, rowFirst, colFirst)

def _dyaddown_1D_keep_odd(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        dyaddown_1D_keep_odd (<double*>input.data, m1*n1, <double*>output.data, m2*n2);
        
def _dyaddown_1D_keep_even(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        dyaddown_1D_keep_even (<double*>input.data, m1*n1, <double*>output.data, m2*n2);
        
def _dyaddown_2D_keep_odd_row(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_odd_row (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyaddown_2D_keep_odd_col(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_odd_col (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyaddown_2D_keep_even_row(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_even_row (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyaddown_2D_keep_even(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_even (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyaddown_2D_keep_even_col(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_even_col (<double*>input1.data, m1, n1, <double*>output.data, m2, n2)


def _dyaddown_2D_keep_odd(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyaddown_2D_keep_odd (<double*>input1.data, m1, n1, <double*>output.data, m2, n2)


def _dyadup_1D_feed_odd(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        dyadup_1D_feed_odd (<double*>input.data, m1*n1, <double*>output.data, m2*n2);
        
def _dyadup_1D_feed_even(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        dyadup_1D_feed_even (<double*>input.data, m1*n1, <double*>output.data, m2*n2);
        
def _dyadup_2D_feed_odd_row(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_odd_row (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyadup_2D_feed_odd_col(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_odd_col (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyadup_2D_feed_even_row(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_even_row (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyadup_2D_feed_even(np.ndarray[np.float64_t, ndim=2] input, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input.shape[0]
        n1 = input.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_even (<double*>input.data, m1, n1, <double*>output.data, m2, n2)

def _dyadup_2D_feed_even_col(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_even_col (<double*>input1.data, m1, n1, <double*>output.data, m2, n2)


def _dyadup_2D_feed_odd(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]

        dyadup_2D_feed_odd (<double*>input1.data, m1, n1, <double*>output.data, m2, n2)



def _wextend_1D_center(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output, char* extModeStr):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wextend_1D_center(<double*>input.data, m1*n1, <double*>output.data, m2*n2, char_to_extend_method (extModeStr));


def _wextend_1D_left(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output, char* extModeStr):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]

        wextend_1D_left(<double*>input.data, m1*n1, <double*>output.data, m2*n2, char_to_extend_method (extModeStr));


def _wextend_1D_right(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output, char* extModeStr):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]
        
        wextend_1D_right(<double*>input.data, m1*n1, <double*>output.data, m2*n2, char_to_extend_method (extModeStr));


def _wextend_2D(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output, char* extModeStr, char *rowOpt, char *colOpt):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]
        wextend_2D (<double*>input1.data, m1, n1, <double*>output.data, m2, n2, char_to_extend_method (extModeStr), rowOpt, colOpt)


def _wextend_2D_col(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output, char* extModeStr, char *Opt):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]
        wextend_2D_col (<double*>input1.data, m1, n1, <double*>output.data, m2, n2, char_to_extend_method (extModeStr), Opt)


def _wextend_2D_row(np.ndarray[np.float64_t, ndim=2] input1, np.ndarray[np.float64_t, ndim=2] output, char* extModeStr, char *Opt):
        cdef int m1, n1
        m1 = input1.shape[0]
        n1 = input1.shape[1]
        cdef int m2, n2
        m2 = output.shape[0]
        n2 = output.shape[1]
        wextend_2D_row (<double*>input1.data, m1, n1, <double*>output.data, m2, n2, char_to_extend_method (extModeStr), Opt)


def _waverec(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output,np.ndarray[np.float64_t, ndim=1] lowRe, np.ndarray[np.float64_t, ndim=1] hiRe, np.ndarray[np.int32_t, ndim=1] waveDecLengthArray, stride):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output.shape[0]
        cdef int filterLength
        filterLength = lowRe.shape[0]
        cdef int waveDecLength
        waveDecLength = waveDecLengthArray.shape[0]
        
        waverec (<double*>input.data, m1*n1, <double*>output.data, m2*n2, <double*>lowRe.data, <double*>hiRe.data, filterLength, <int*>waveDecLengthArray.data, waveDecLength, stride, getdwtMode())
        filter_clear()


def _wavedec(np.ndarray[np.float64_t, ndim=1] input, np.ndarray[np.float64_t, ndim=1] output1, np.ndarray[np.float64_t, ndim=1] lowRe, np.ndarray[np.float64_t, ndim=1] hiRe, np.ndarray[np.int32_t, ndim=1] output2, stride):
        cdef int m1, n1
        m1 = 1
        n1 = input.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = output1.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output2.shape[0]
        cdef int filterLength
        filterLength = lowRe.shape[0]
        
        wavedec (<double*>input.data, m1*n1, <double*>output1.data, m2*n2, <double*>lowRe.data, <double*>hiRe.data, filterLength, <int*>output2.data, m5*n5, stride, getdwtMode())
        filter_clear()
        
def _wave_dec_len_cal(filterLen, sigLength, stride, np.ndarray[np.int32_t, ndim=1] waveDecLengthArray):
        cdef int waveDecLength
        waveDecLength = waveDecLengthArray.shape[0]
        wave_dec_len_cal(filterLen, sigLength, stride, <int*>waveDecLengthArray.data)
        #count = 0
        #waveDecLengthArray[stride + 1] = sigLength
        #if (getdwtMode()!=PER):
        #        calLen = sigLength
        #        for count in np.arange(stride):
        #                calLen += (filterLen - 1)
        #                waveDecLengthArray[stride-count]=int(np.floor(calLen/2))
        #                calLen = waveDecLengthArray[stride - count]
        #        waveDecLengthArray[0] = waveDecLengthArray[1]
        #else:
        #        for count in np.arange(count,0,-1):
        #                waveDecLengthArray[count] = int(np.ceil(((waveDecLengthArray[count+1]))/2.0))
        #        waveDecLengthArray[0] = waveDecLengthArray[1]

def _idwt_neo(np.ndarray[np.float64_t, ndim=1] cA, np.ndarray[np.float64_t, ndim=1] cD, np.ndarray[np.float64_t, ndim=1] lowRe, np.ndarray[np.float64_t, ndim=1] hiRe, np.ndarray[np.float64_t, ndim=1] output1):
        cdef int m1, n1
        m1 = 1
        n1 = cA.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = cD.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output1.shape[0]
        cdef int filterLength
        filterLength = lowRe.shape[0]
        
        idwt_neo (<double*>cA.data, <double*>cD.data, m1*n1, <double*>lowRe.data, <double*>hiRe.data, filterLength, <double*>output1.data, m5*n5)
        filter_clear()

def _idwt_detail_neo(np.ndarray[np.float64_t, ndim=1] cD, np.ndarray[np.float64_t, ndim=1] hiRe, np.ndarray[np.float64_t, ndim=1] output1):
        cdef int m2, n2
        m2 = 1
        n2 = cD.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output1.shape[0]
        cdef int filterLength
        filterLength = hiRe.shape[0]
        
        idwt_detail_neo (<double*>cD.data, m2*n2, <double*>hiRe.data, filterLength, <double*>output1.data, m5*n5)
        filter_clear()

def _idwt_approx_neo(np.ndarray[np.float64_t, ndim=1] cA, np.ndarray[np.float64_t, ndim=1] lowRe,  np.ndarray[np.float64_t, ndim=1] output1):
        cdef int m1, n1
        m1 = 1
        n1 = cA.shape[0]
        cdef int m5, n5
        m5 = 1
        n5 = output1.shape[0]
        cdef int filterLength
        filterLength = lowRe.shape[0]
        
        idwt_approx_neo (<double*>cA.data, m1*n1, <double*>lowRe.data,  filterLength, <double*>output1.data, m5*n5)
        filter_clear()

def _detcoef(np.ndarray[np.float64_t, ndim=1] C, np.ndarray[np.int32_t, ndim=1] L,  np.ndarray[np.float64_t, ndim=1] output1, stride, level):
        cdef int m1, n1
        m1 = 1
        n1 = C.shape[0]
        cdef int m2, n2
        m2 = 1
        n2 = L.shape[0]
        cdef int m5, n5
        m4 = 1
        n4 = output1.shape[0]
        
        detcoef(<double*>C.data, m1*n1, <int*>L.data, m2*n2, <double*>output1.data, m4*n4, stride, level)
