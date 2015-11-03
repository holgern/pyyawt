# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt

__doc__ = """Pyrex wrapper for low-level C edflib implementation."""
__all__ = ["PYYAWT_HAAR", "PYYAWT_DAUBECHIES", "PYYAWT_COIFLETS", "PYYAWT_SYMLETS",
        "PYYAWT_SPLINE_BIORTH", "PYYAWT_BEYLKIN", "PYYAWT_VAIDYANATHAN", "PYYAWT_DMEY",
        "PYYAWT_BATHLETS", "PYYAWT_LEGENDRE", "PYYAWT_SPLINE_RBIORTH", "PYYAWT_FARRAS",
        "PYYAWT_KINGSBURYQ", "PYYAWT_NOT_DEFINED", "PYYAWT_ORTH", "PYYAWT_BIORTH",
        "_wavelet_parser", "_dbwavf_length", "_dbwavf","_coifwavf",
        "_coifwavf_length", "_symwavf", "_symwavf_length", "_legdwavf",
        "_legdwavf_length"]


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


#cdef class CyDaubechies:
        #cdef swt_wavelet *_thisptr
        #cdef int family
        #cdef int member
        #def __cinit__(self, char *wname):
                #self.family = -1
                #self.member = -1
                #wavelet_parser(wname, &self.family, &self.member)
                #self._thisptr = <swt_wavelet *> malloc (sizeof (swt_wavelet))
                #daubechies_synthesis_initialize(self.member,self._thisptr)
                #if self._thisptr == NULL:
                        #msg = "Insufficient memory."
                        #raise MemoryError(msg)
                        
        #def __dealloc__(self):
                #if self._thisptr is not NULL:
                        #free (self._thisptr)
                        #self._thisptr = NULL

        #def getLength(self):
                #if self._thisptr is not NULL:
                        #return self._thisptr.length
            
        #def getLowPass(self, np.ndarray[np.float64_t, ndim=1] sigbuf):
                #if self._thisptr is not NULL:
                        #m2 = 1;
                        #n2 = self.getLength()
                        #verbatim_copy (self._thisptr.pLowPass, m2*n2, <double*>sigbuf.data, m2*n2)
