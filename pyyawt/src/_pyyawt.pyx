# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt

__doc__ = """Pyrex wrapper for low-level C edflib implementation."""
__all__ = ["PYYAWT_HAAR", "PYYAWT_DAUBECHIES", "PYYAWT_COIFLETS", "PYYAWT_SYMLETS",
        "PYYAWT_SPLINE_BIORTH", "PYYAWT_BEYLKIN", "PYYAWT_VAIDYANATHAN", "PYYAWT_DMEY",
        "PYYAWT_BATHLETS", "PYYAWT_LEGENDRE", "PYYAWT_SPLINE_RBIORTH", "PYYAWT_FARRAS",
        "PYYAWT_KINGSBURYQ", "PYYAWT_NOT_DEFINED", "PYYAWT_ORTH", "PYYAWT_BIORTH",
        "_wavelet_parser", "_dbwavf_length", "_dbwavf","_coifwavf",
        "_coifwavf_length", "_symwavf", "_symwavf_length", "_legdwavf",
        "_legdwavf_length", "_biorwavf", "_biorwavf_length", "_rbiorwavf", "_rbiorwavf_length"]


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


def _wavelet_parser(char *wname):
        cdef int ret_family = -1
        cdef int ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        return ret_family, ret_member

def _dbwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        daubechies_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length

def _dbwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        daubechies_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)

def _coifwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        coiflets_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
  
def _coifwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        coiflets_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length

def _symwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        symlets_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
  
def _symwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        symlets_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length

def _legdwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        legendre_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)

def _legdwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        legendre_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length

def _biorwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf1, np.ndarray[np.float64_t, ndim=1] sigbuf2):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        sp_bior_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf1.data, m2*n2)
        
        sp_bior_analysis_initialize(ret_member,&thisptr)
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf2.data, m2*n2)

def _biorwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        sp_bior_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length

def _rbiorwavf(char *wname, np.ndarray[np.float64_t, ndim=1] sigbuf1, np.ndarray[np.float64_t, ndim=1] sigbuf2):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr

        sp_rbior_synthesis_initialize(ret_member,&thisptr)
        m2 = 1;
        n2 = thisptr.length
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf1.data, m2*n2)
        
        sp_rbior_analysis_initialize(ret_member,&thisptr)
        verbatim_copy (thisptr.pLowPass, m2*n2, <double*>sigbuf2.data, m2*n2)

def _rbiorwavf_length(char *wname):
        cdef int ret_family
        cdef int ret_member
        ret_family = -1
        ret_member = -1
        wavelet_parser(wname, &ret_family, &ret_member)
        cdef swt_wavelet thisptr
        sp_rbior_synthesis_initialize(ret_member,&thisptr)
        return thisptr.length
