# -*- coding: utf-8 -*-
# Copyright (c) 2015 Holger Nahrstaedt

# get constants
include "pyyawt.pxi"

cdef extern from "swtlib.h":
      cdef enum extend_method:
            ZPD, SYMH, SYMW, ASYMH, ASYMW,
            SP0, SP1, PPD, PER
           
            
            
      ctypedef struct  extension_identity:
            char extMethodName[6]
            extend_method extMethod
            
      # ctypedef struct  str_error_notification:
      #      int   errorNumber;
      #      char  message[150];
            
      ctypedef void(*WScaleFunc)(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      
      ctypedef struct  cwt_identity:
            char wname[20];
            int     realOrComplex;
            int     family;
            int     phipsi;
            double  lb;
            double  ub;
            double cpsi;
            WScaleFunc scalef;
            
      ctypedef struct  cwt_family:
            char wname[20];
            char     realOrComplex[20];
            char     family[20];
            
      ctypedef struct  swt_wavelet:
            int     length;
            double  *pLowPass;
            double  *pHiPass;
      
      
      ctypedef void(*Func)(int member, swt_wavelet *pWaveStruct);
            
      ctypedef struct  wavelet_identity:
            char  wname[20];
            int   rOrB;
            int   family;
            int   member;
            Func  analysis;
            Func  synthesis;
            
      ctypedef struct   wavelet_family:
            char  wname[20];
            char   rOrB[20];
            char   family[20];
            


            
      void sinus(double *x, int sigInLength, double *psi, int sigOutLength, double ys);

      void poisson(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void mexihat(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      
      void morlet(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      
      void DOGauss(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      
      void Gauss(double *x, int sigInLength, double *psi, int sigOutLength, int n, double ys);
      void Gaus1(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus2(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus3(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus4(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus5(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus6(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus7(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      void Gaus8(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
      
      void cgauss(double *x, int sigInLength, int p, double *psir, double *psii, int sigOutLength, double ys);
      void cgau1(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau1_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau2(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau2_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau3(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau3_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau4(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau4_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      
      void cgau5(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau5_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau6(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau6_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau7(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau7_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      void cgau8(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cgau8_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
      
      
      void cmorlet(double *x, int sigInLength, double fb, double fc, double *psir, double *psii,  int sigOutLength, double ys);
      void cmorlet_packet(double *x, int sigInLength,  double *f, int sigOutLength, double ys);
      
      void shanwavf(double *x, int sigInLength,  double fb, double fc, double *psir, double *psii,   int sigOutLength, double ys);
      void shanwavf_packet(double *x, int sigInLength,  double *f, int sigOutLength, double ys);
      void fbspwavf(double *x, int sigInLength,int m,  double fb, double fc, double *psir, double *psii,   int sigOutLength, double ys);
      void fbspwavf_packet(double *x, int sigInLength,  double *f, int sigOutLength, double ys);
      
      void cauchy(double *x, int sigInLength, double fb, double fc, double *psir, double *psii,  int sigOutLength, double ys);
      void cauchy_neo(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
      void cauchy_packet(double *x, int sigInLength,  double *f, int sigOutLength, double ys);
      
      
      void meyeraux(double x, double *y);
      void meyer_phi(double *x, int sigInLength,    double lb, double ub, double *phir, double *phii,    int sigOutLength, double ys);
       
       
      void haar_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void haar_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void daubechies_analysis_initialize (int memeber, swt_wavelet *pWaveStruct);
      void daubechies_synthesis_initialize (int memeber, swt_wavelet *pWaveStruct);
      void symlets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void symlets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void coiflets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void coiflets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void sp_bior_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void sp_bior_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void sp_rbior_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void sp_rbior_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void beylkin_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void beylkin_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void vaidyanathan_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void vaidyanathan_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void dmey_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void dmey_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void bathlets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void bathlets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void legendre_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void legendre_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void farras_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void farras_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
      void kingsburyq_analysis_initialize (int member, swt_wavelet *pWaveStruct);
      void kingsburyq_synthesis_initialize (int member, swt_wavelet *pWaveStruct);

      void  cowavedec (double *sigIn, int sigInLength, double *sigOutR,  double *sigOutI, int sigOutLength,        double *lowDTree1S1, double *hiDTree1S1,     double *lowDTree2S1, double *hiDTree2S1,   double *lowDTree1S2, double *hiDTree1S2,     double *lowDTree2S2, double *hiDTree2S2,  int filterLen, int *waveDecLengthArray,           int lengthArrayLengh, int stride, extend_method extMethod);

      void  cowaverec (double *sigInR, double *sigInI, int sigInLength,    double *sigOut, int sigOutLength,  double *lowRTree1S1, double *hiRTree1S1,    double *lowRTree2S1, double *hiRTree2S1,     double *lowRTree1S2, double *hiRTree1S2,         double *lowRTree2S2, double *hiRTree2S2,          int filterLen, int *waveDecLengthArray,          int lengthArraylength, int stride,          extend_method extMethod);

      void  cowavedec2 (double *matrixIn, int matrixInRow, int matrixInCol,            double *lowDTree1S1, double *hiDTree1S1,            double *lowDTree1S2, double *hiDTree1S2,            int filterLen, int *pLen, double *coef,            int sigOutLength, int stride, extend_method extMethod);
      void  cowavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,            double *lowDTree1S1R, double *hiDTree1S1R,            double *lowDTree1S1C, double *hiDTree1S1C,            double *lowDTree1S2R, double *hiDTree1S2R,            double *lowDTree1S2C, double *hiDTree1S2C,            int filterLen, int *pLen, double *coef,            int sigOutLength, int stride, extend_method extMethod);

      void  cowaverec2 (double *coef, int sigInLength,            double *lowRTree1S1, double *hiRTree1S1,            double *lowRTree1S2, double *hiRTree1S2,   int filterLen, double *matrixOut, int matrixOutRow,  int matrixOutCol, int *pLen, int stride,   extend_method extMethod);

      void cowaverec2a (double *coef, int sigInLength,            double *lowRTree1S1R, double *hiRTree1S1R,            double *lowRTree1S1C, double *hiRTree1S1C,            double *lowRTree1S2R, double *hiRTree1S2R,            double *lowRTree1S2C, double *hiRTree1S2C,            int filterLen, double *matrixOut, int matrixOutRow,            int matrixOutCol, int *pLen, int stride,            extend_method extMethod);

      void  copmd (double *matrixInR, double *matrixInI, int sigInLength,            int InRow, int InCol, double *matrixOutR, double *matrixOutI);

      void  copmr (double *matrixInR, double *matrixInI, int sigInLength,            int InRow, int InCol, double *matrixOutR, double *matrixOutI);

      void filter_clear ();
      void orth_filt_group (double *filterIn, int sigInLength,            double *filterLowRec,            double *filterLowDec,            double *filterHiRec,            double *filterHiDec);             
      void bior_filt_group (double *f1, int sigInLength1,            double *f2, int sigInLength2,            double *lowDecom, int sigOutLength1,            double *hiDecom, int sigOutLength2,            double *lowRecon, int sigOutLength3,            double *hiRecon, int sigOutLength4);
      void wavelet_parser (char *wname, int *family, int *member);
      void wavelet_fun_parser (char *wname, int *ii);
      void wave_len_validate (int sigInLen, int waveLength, int *lev, int *val);

      extend_method getdwtMode();
      void setdwtMode(extend_method mode);
      void dwt_write (char *mode, int *errCode);
      extend_method char_to_extend_method(char *mode);
      void dwt_parse(char **strr);
      void dwt (double *sigIn, int sigInLength, double *lowDe,            double *hiDe, int filterLen, double *approx,            double *detail, int sigOutLength,            extend_method extMethod);             
      void dwt_neo (double *sigIn, int sigInLength, double *lowDe,            double *hiDe, int filterLen, double *approx,            double *detail, int sigOutLength,            extend_method extMethod);
      void dwt_nex (double *sigIn, int sigInLength, double *lowDe,            double *hiDe, int filterLen, double *approx,            double *detail, int sigOutLength);
      void dwt_no_extension (double *sigIn, int sigInLength, double *lowDe,            double *hiDe, int filterLen, double *approx,            double *detail, int sigOutLength);
      void dwt_conv (double *sigIn, int sigInLength, double *lowDe,            double *hiDe, int filterLen, double *approx,            double *detail, int sigOutLength);
      void idwt_complete (double *approx, double *detail,            int sigInLength, double *lowRe,           double *hiRe, int filterLen,            double *sigOut, int sigOutLength);
      void idwt_neo (double *approx, double *detail,            int sigInLength, double *lowRe,            double *hiRe, int filterLen,            double *sigOut, int sigOutLength);
      void idwt_complete_ex (double *approx, double *detail,            int sigInLength, double *lowRe,            double *hiRe, int filterLen,           double *sigOut, int sigOutLength,            extend_method extMethod);
      void idwt_approx (double *approx, int sigInLength,            double *lowRe, int filterLen,            double *sigOut, int sigOutLength);
      void idwt_approx_ex (double *approx, int sigInLength,            double *lowRe, int filterLen,            double *sigOut, int sigOutLength,            extend_method extMethod);
      void idwt_approx_neo (double *approx, int sigInLength,            double *lowRe, int filterLen,            double *sigOut, int sigOutLength);
      void idwt_detail (double *detail, int sigInLength,            double *hiRe, int filterLen,            double *sigOut, int sigOutLength);
      void idwt_detail_ex (double *detail, int sigInLength,            double *hiRe, int filterLen,            double *sigOut, int sigOutLength,            extend_method extMethod);
      void idwt_detail_neo (double *detail, int sigInLength,            double *hiRe, int filterLen,            double *sigOut, int sigOutLength);
      void wave_dec_len_cal (int filterLen, int sigLength,            int stride, int *waveDecLengthArray);
      void wavedec (double *sigIn, int sigInLength, double *sigOut,            int sigOutLength, double *lowDe, double *hiDe,            int filterLen, int *waveDecLengthArray,            int lengthArrayLengh, int stride,            extend_method extMethod);
      void waverec (double *sigIn, int sigInLength, double *sigOut,            int sigOutLength, double *lowRe, double *hiRe,            int filterLen, int *waveDecLengthArray,            int lengthArraylength, int stride,            extend_method extMethod);
      void wenergy (double *coef, int coefLen, int *lenArray,            int arrayLen, double *aE, int aELen,            double *dE, int dELen);
      void detcoef (double *sigIn, int sigInLength,            int *waveDecLengthArray, int arrayLen,            double *sigOut, int sigOutLength,            int stride, int level);
      void appcoef (double *sigIn, int sigInLength, double *sigOut,            int sigOutLength, double *lowRe, double *hiRe,            int filterLen, int *waveDecLengthArray,            int lengthArraylength, int stride, int level,            extend_method extMethod);
      void wrcoef (double *sigIn, int sigInLength, double *lowRe,            double *hiRe, int filterLen,            int *waveDecLengthArray, int arrayLen,            double *sigOut, int sigOutLength,            char *coefType, int stride, int level,            extend_method extMethod);
      void upcoef_len_cal (int sigInLength, int filterLen,            int stride, int *sigOutLength,            int *sigOutLengthDefault);
      void upwlev (double *coefArray, int coefLen,            int *waveDecLengthArray,	int arrayLen,            double *lowRe, double *hiRe, int filterLen,            double *newCoefArray, int newCoefLen,            int *newLenArray, int newArrayLen,            double *approx, int approxLen, int stride,            extend_method extMethod);
      void upcoef (double *sigIn, int sigInLength, double *lowRe,            double *hiRe, int filterLen, double *sigOut,            int sigOutLength, int defaultLength,            char *coefType, int step);

      void dwt2D (double *matrixIn, int matrixInRow,            int matrixInCol, double *matrixOutApprox,            double *matrixOutColDetail,            double *matrixOutRowDetail,            double *matrixOutDetail, int matrixOutRow,            int matrixOutCol, double *lowDe, double *hiDe,            int filterLen, extend_method extMethod);
      void            dwt2D_neo_a (double *matrixIn, int matrixInRow, int matrixInCol,            double *matrixOutApprox, double *matrixOutColDetail,            double *matrixOutRowDetail, double *matrixOutDetail,            int matrixOutRow, int matrixOutCol, double *lowDeR,            double *hiDeR, double *lowDeC, double *hiDeC,            int filterLen, extend_method extMethod);
      void dwt2D_neo (double *matrixIn, int matrixInRow,            int matrixInCol, double *matrixOutApprox,            double *matrixOutColDetail,            double *matrixOutRowDetail,            double *matrixOutDetail, int matrixOutRow,            int matrixOutCol, double *lowDe, double *hiDe,            int filterLen, extend_method extMethod);             
      void idwt2D (double *matrixInApprox,            double *matrixInColDetail,            double *matrixInRowDetail,            double *matrixInDetail,            int matrixInRow, int matrixInCol, double *lowRe,            double *hiRe, int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol,            extend_method extMethod);
      void idwt2D_neo (double *matrixInApprox, double *matrixInColDetail,            double *matrixInRowDetail, double *matrixInDetail,            int matrixInRow, int matrixInCol, double *lowRe,            double *hiRe, int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol);
      void            idwt2D_neo_a (double *matrixInApprox, double *matrixInColDetail,            double *matrixInRowDetail, double *matrixInDetail,            int matrixInRow, int matrixInCol, double *lowReR,            double *hiReR, double *lowReC, double *hiReC,            int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol);

      void wave_mem_cal (int *pLen, int stride, int *total);
      void matrix_wavedec_len_cal (int matrixInRow, int matrixInCol,            int stride, int filterLen,            int *pLen);
      void matrix_locate (int stride, int *pLen, int *pH,            int *pV, int *pD);
      void wavedec2 (double *matrixIn, int matrixInRow,            int matrixInCol, double *lowDe, double *hiDe,            int filterLen, int *pLen, double *coef,            int sigOutLength, int stride,            extend_method extMethod);
      void            wavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,            double *lowDeR, double *hiDeR, double *lowDeC,            double *hiDeC, int filterLen, int *pLen,            double *coef, int sigOutLength, int stride,            extend_method extMethod);
      void waverec2 (double *coef, int sigInLength, double *lowRe,            double *hiRe, int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol, int *pLen,            int stride, extend_method extMethod);
      void            waverec2a (double *coef, int sigInLength, double *lowReR,            double *hiReR, double *lowReC, double *hiReC,            int filterLen, double *matrixOut, int matrixOutRow,            int matrixOutCol, int *pLen, int stride,            extend_method extMethod);

      void wenergy_2output (double *coef, int sigInLength,            int *pLen, double *ae, double *de,            int deLength, int stride);
      void wenergy_4output (double *coef, int sigInLength,            int *pLen, double *ae, double *he,            double *ve, double *de, int deLength,            int stride);
      void detcoef2 (double *coef, int sigInLength, double *coefOut,            int sigOutLength, int *pLen, int stride,            int level, char *coefType);
      void appcoef2 (double *coef, int sigInLength, double *lowRe,            double *hiRe, int filterLen, double *coefOut,            int matrixOutRow, int matrixOutCol, int *pLen,           int stride, int level, extend_method extMethod);
      void wrcoef2 (double *coef, int sigInLength, double *lowRe,            double *hiRe, int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol, int *pLen,            int stride, int level, char *type,            extend_method extMethod);
      void upwlev2 (double *coef, int sigInLength, double *lowRe,            double *hiRe,            int filterLen, int *pLen, int matrixRow, int matrixCol,            double *approx, int approxLen, double *newCoef,            int newCoefLen, int *newLenMatrix, int lenMatrixRow,            int lenMatrixCol, int stride, extend_method extMethod);
      void upcoef2 (double *matrixIn, int matrixInRow,            int matrixInCol, double *lowRe, double *hiRe,            int filterLen, double *matrixOut,            int matrixOutRow, int matrixOutCol,            int matrixOutDefaultRow,            int matrixOutDefaultCol,            int step, char *type);

      void dwt3d_tran(double *mat3DIn, int row1, int col1, int sli1,            double *mat3DOut, int row2, int col2, int sli2);

      void dwt3d_line_forward(double *mat3DIn, int row1, int col1, int sli1,            double *mat3DOutApp, double *mat3DOutDet,            int row2, int col2, int sli2,            double *loDe, double *hiDe, int filterLen,            extend_method extMethod);

      void dwt3d_tran_z(double *mat3DIn, int row1, int col1, int sli1,            double *mat3DOut, int row2, int col2, int sli2);
      void dwt3d_tran_z_inv(double *mat3DIn, int row1, int col1, int sli1,            double *mat3DOut, int row2, int col2, int sli2);

      void dwt3d_combine(double *mat1, double *mat2, double *mat3,            double *mat4, double *mat5, double *mat6,            double *mat7, double *mat8, int rowIn,            int colIn, int sliIn, double *matOut,            int rowOut, int colOut, int sliOut);
      void dwt3d_line_reverse(double *mat3DInApp, double *mat3DInDet,            int row1, int col1, int sli1,            double *mat3DOut, int row2, int col2,            int sli2,            double *loDe, double *hiDe, int filterLen);


      void dwt3d_split(double *matIn, int rowIn, int colIn, int sliIn,            double *mat1, double *mat2, double *mat3,            double *mat4, double *mat5, double *mat6,            double *mat7, double *mat8, int rowOut,            int colOut, int sliOut);
      void dwt3(double *mat3DIn, int row, int col, int sli,            double *mat3DOut, int row2, int col2, int sli2,            int r, int c, int s, double *Lo1, double *Hi1,            double *Lo2, double *Hi2, double *Lo3, double *Hi3,            int fLen1, int fLen2, int fLen3, extend_method extMethod);

      void cwt_fun_parser(char *wname, int *ind);
      void cwt_len_cal (int sigInLength, int scale, int *sigOutLength, double *delta);
      void full_range_scalef (char *wname, double *f, int sigOutLength);
      void scale_real (double *f, int sigInLength, double delta, double *fout, int sigOutLength);
      void cwt_conv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
      void cwt_iconv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
      void cwt_conv_complex (double *sigIn, int sigInLength, double *fr, double *fi, int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);
      void cwt_conv_complex_complex (double *a, double *b, int sigInLength,double *c, double *d,  int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);
      void swt_conv(double *sigIn, int sigInLength,             double *approx, int approxLength,             double *detail, int detailLength,             double *filterLow, double *filterHi,             int filterLength);
      void  swt_out1 (double *sigIn, int sigInLength,             double *sigOutMatrix, int rowLength,             int colLength, double *filterLow,             double *filterHi, int filterLength, int step);

      void swt_out2 (double *sigIn, int sigInLength,             double *approxMatrix, double *detailMatrix,             int rowLength, int colLength, double *filterLow,             double *filterHi, int filterLength, int step);

      void iswt_conv (double *approx, double *detail, int sigInLength,             double *sigOut, int sigOutLength, double *filterLow,             double *filterHi, int filterLength);              
      void iswt_conv_step (double *approx, double *detail, int sigInLength,
      double *sigOut, int sigOutLength, double *filterLow,             double *filterHi, int filterLength, int level);

      void iswt_input1 (double *matrixIn, int rowLength, int colLength,             double *sigOut, int sigOutLength, double *filterLow,             double *filterHi, int filterLength);

      void iswt_input2 (double *matrixApproxIn, double *matrixDetailIn,             int rowLength, int colLength,             double *sigOut, int sigOutLength, double *filterLow,             double *filterHi, int filterLength);
      void swt2_output4(double *matrixIn, int matrixInRow, int matrixInCol,             double *matrixOutApprox, double *matrixOutColDetail,             double *matrixOutRowDetail, double *matrixOutDetail,            int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);
      void swt2_output4_step(double *matrixIn, int matrixInRow, int matrixInCol,             double *matrixOutApprox, double *matrixOutColDetail,             double *matrixOutRowDetail, double *matrixOutDetail,             int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);
      void swt2_output1_step(double *matrixIn, int matrixInRow,
      int matrixInCol,  double *matrixOut,             int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);              
      void iswt2(double *matrixInApprox, double *matrixInColDetail,             double *matrixInRowDetail, double *matrixInDetail,             int matrixInRow, int matrixInCol,             double *matrixOut, int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);
      void iswt2_input4_step(double *matrixInApprox, double *matrixInColDetail,             double *matrixInRowDetail, double *matrixInDetail,             int matrixInRow, int matrixInCol,             double *matrixOut, int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);
      void iswt2_input1_step(double *matrixIn,  int matrixInRow, int matrixInCol,             double *matrixOut, int matrixOutRow, int matrixOutCol,             double *filterLow, double *filterHi,             int filterLength, int step);
             
      void matrix_tran (double *matrixIn, int matrixInRow,                                 int matrixInCol, double *matrixOut,                                 int matrixOutRow, int matrixOutCol);             
      void wrev (const double *sigIn, int sigInLength,                          double *sigOut, int sigOutLength);
      void qmf_even (const double *sigIn, int sigInLength,                              double *sigOut, int sigOutLength);
      void qmf_odd (double *sigIn, int sigInLength,                             double *sigOut, int sigOutLength);
      void qmf_wrev (const double *sigIn, int sigInLength,                              double *sigOut, int sigOutLength);
      void verbatim_copy (const double *sigIn, int sigInLength,                                   double *sigOut, int sigOutLength);
      void dyaddown_1D_keep_odd (double *sigIn, int sigInLength,                                       double *sigOut, int sigOutLength);
      void dyaddown_1D_keep_even (double *sigIn, int sigInLength,                                        double *sigOut, int sigOutLength);
      void dyaddown_2D_keep_odd_row (double *matrixIn,                                              int matrixInRow,                                              int matrixInCol,                                              double *matrixOut,                                              int matrixOutRow,                                              int matrixOutCol);
      void dyaddown_2D_keep_odd_col (double *matrixIn,                                              int matrixInRow,                                              int matrixInCol,                                              double *matrixOut,                                              int matrixOutRow,                                              int matrixOutCol);
      void dyaddown_2D_keep_even_row (double *matrixIn,                                               int matrixInRow,                                               int matrixInCol,                                               double *matrixOut,                                               int matrixOutRow,                                               int matrixOutCol);
      void dyaddown_2D_keep_even_col (double *matrixIn,                                               int matrixInRow,                                               int matrixInCol,                                               double *matrixOut,                                               int matrixOutRow,                                               int matrixOutCol);
      void dyaddown_2D_keep_odd (double *matrixIn,                                          int matrixInRow,                                          int matrixInCol,                                          double *matrixOut,                                          int matrixOutRow,                                          int matrixOutCol);
      void dyaddown_2D_keep_even (double *matrixIn,                                           int matrixInRow,                                           int matrixInCol,                                           double *matrixOut,                                           int matrixOutRow,                                           int matrixOutCol);
      void dyadup_1D_feed_odd (double *sigIn, int sigInLength,                                        double *sigOut, int sigOutLength);              
      void dyadup_1D_feed_even (double *sigIn, int sigInLength,                                         double *sigOut, int sigOutLength);
      void dyadup_2D_feed_odd_row (double *matrixIn,                                            int matrixInRow,                                            int matrixInCol,                                            double *matrixOut,                                            int matrixOutRow,                                            int matrixOutCol);
      void dyadup_2D_feed_odd_col (double *matrixIn,                                            int matrixInRow,                                            int matrixInCol,                                            double *matrixOut,                                            int matrixOutRow,                                            int matrixOutCol);
      void dyadup_2D_feed_even_row (double *matrixIn,                                             int matrixInRow,                                             int matrixInCol,                                             double *matrixOut,                                             int matrixOutRow,                                             int matrixOutCol);
      void dyadup_2D_feed_even_col (double *matrixIn,                                             int matrixInRow,                                             int matrixInCol,                                             double *matrixOut,                                             int matrixOutRow,                                             int matrixOutCol);
      void dyadup_2D_feed_odd (double *matrixIn,                                        int matrixInRow,                                        int matrixInCol,                                        double *matrixOut,                                        int matrixOutRow,                                        int matrixOutCol);
      void dyadup_2D_feed_even (double *matrixIn,                                         int matrixInRow,                                         int matrixInCol,                                         double *matrixOut,                                         int matrixOutRow,                                         int matrixOutCol);
      void extend_method_parse (char *mode, extend_method *extMethod);              
      void wextend_1D_center (double *sigIn, int sigInLength,                                       double *sigOut, int sigOutLength,                                       extend_method method);
      void wextend_1D_left (double *sigIn, int sigInLength,                                       double *sigOut, int sigOutLength,                                       extend_method method);
      void wextend_1D_right (double *sigIn, int sigInLength,                                       double *sigOut, int sigOutLength,                                       extend_method method);
      void wextend_2D (double *matrixIn, int matrixInRow,                                int matrixInCol, double *matrixOut,                                int matrixOutRow, int matrixOutCol,                                extend_method extMethod, char *rowOpt,                                char *colOpt);
      void wextend_2D_row (double *matrixIn, int matrixInRow,                                    int matrixInCol, double *matrixOut,                                    int matrixOutRow, int matrixOutCol,                                    extend_method extMethod, char *Opt);
      void wextend_2D_col (double *matrixIn, int matrixInRow,                                    int matrixInCol, double *matrixOut,                                    int matrixOutRow, int matrixOutCol,                                    extend_method extMethod, char *Opt);
      void wkeep_1D_center (double *sigIn, int sigInLength,                                     double *sigOut, int sigOutLength);              
      void wkeep_1D_left (double *sigIn, int sigInLength,                                   double *sigOut, int sigOutLength);
      void wkeep_1D_right (double *sigIn, int sigInLength,                                    double *sigOut, int sigOutLength);
      void wkeep_1D_index (double *sigIn, int sigInLength,                                   double *sigOut, int sigOutLength,                                    int first);
      void wkeep_2D_center (double *matrixIn, int matrixInRow,                                     int matrixInCol, double *matrixOut,                                     int matrixOutRow, int matrixOutCol);
      void wkeep_2D_index (double *matrixIn, int matrixInRow,                                    int matrixInCol, double *matrixOut,                                    int matrixOutRow, int matrixOutCol,                                    int rowFirst, int colFirst);
      void conv (double *sigIn, int sigInLength,                          double *sigOut, int sigOutLength,                          double *fiter, int filterLength);
      void i_conv (double *sigIn, int sigInLength,                          double *sigOut, int sigOutLength,                          double *fiter, int filterLength);               
              
      void swt_exp2(int lev, int *outputV);
      void linspace(double lb, double ub, int n, double *sigOut, int sigOutLength);
      void ocumsum (double *sigIn, int sigInLength);
      void swt_max(double *sigIn, int sigInLength, double *sigMax);
      void swt_min(double *sigIn, int sigInLength, double *sigMin);
      void wcodemat_abs(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv);
      void swt_max_abs(double *sigIn, int sigInLength, double *sigMax);
      void swt_min_abs(double *sigIn, int sigInLength, double *sigMin);
      double swt_abs(double sigIn);
      void wcodemat(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv);
      void wcodematd(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, double minv, double maxv);
      void wcodemat_matrix (double *matrixIn, int matrixInRow, int matrixInCol,                                   double *matrixOut, int matrixOutRow, int matrixOutCol,                                                   int minv, int maxv, int abso);
      void wcodemat_matrix_col (double *matrixIn, int matrixInRow, int matrixInCol,                                       double *matrixOut, int matrixOutRow, int matrixOutCol,                                                           int minv, int maxv, int abso);
      void wcodemat_matrix_row (double *matrixIn, int matrixInRow, int matrixInCol,                                       double *matrixOut, int matrixOutRow, int matrixOutCol,                                                           int minv, int maxv, int abso);
