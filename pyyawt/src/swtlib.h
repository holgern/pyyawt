/*
 * -------------------------------------------------------------------------
 * swtlib.h --yet another python wavelet toolbox header file
 * PYYAWT - yet another python wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
 * Copyright (C) 2010-2015  Holger Nahrstaedt
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------------
 */
#ifndef _SWTLIB_H_
#define _SWTLIB_H_



  #include <stdio.h>
  #include <stdlib.h>
  #include <stdarg.h>
  #include <math.h>
  /*********************************************
  * Macro
  ********************************************/

  #define SUCCESS      0
  #define DIM_ERR_ONE  1
  #define DIM_ERR_VEC  2
  #define DIM_ERR_MAT  3


  #define POSITIVE_INTEGER_ONLY                               1
  #define LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION          2
  #define SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION            3
  #define OPT_CHAR_NOT_VALID                                  4
  #define EXTENSION_OPT_NOT_VALID                             5
  #define WAVELET_NAME_NOT_VALID                              6
  #define DECOMPOSITION_LEVEL_NOT_VALID                       7
  #define MULTI_DECOM_LEVEL_LESS_THAN_TWO                     8
  #define WRONG_LHS                                           9
  #define UNKNOWN_INPUT_ERR                                   20

  #define PI   3.1415926535897931159980

  /*********************************************
  * Macros CWT
  ********************************************/

  #define REAL    0
  #define COMPLEX 1

  #define PHI_ONLY     0
  #define PSI_ONLY     1
  #define PHI_PSI_BOTH 2

  #define SINUS           0
  #define POISSON         1
  #define MEXICAN_HAT     2
  #define MORLET          3
  #define DOGAUSS         4
  #define CMORLET         5
  #define SHANNON         6
  #define FBSP            7
  #define CAUCHY          8
  #define GAUSS           9
  #define CGAUSS          10

    /*********************************************
    * Macros DWT
    ********************************************/

    #define HAAR           0
    #define DAUBECHIES     1
    #define COIFLETS       2
    #define SYMLETS        3
    #define SPLINE_BIORTH  4
    #define BEYLKIN        5
    #define VAIDYANATHAN   6
    #define DMEY           7
    #define BATHLETS       8
    #define LEGENDRE       9
    #define SPLINE_RBIORTH 10
    #define FARRAS         11
    #define KINGSBURYQ     12
    #define NOT_DEFINED    99

    #define ORTH       0
    #define BIORTH     1

#ifdef __cplusplus
extern "C" {
#endif

  /*********************************************
  * Extension Type
  ********************************************/

  typedef enum {
    ZPD, SYMH, SYMW, ASYMH, ASYMW,
    SP0, SP1, PPD, PER} extend_method;

    /*********************************************
    * Structure Declarations
    ********************************************/
    // #ifndef __USE_DEPRECATED_STACK_FUNCTIONS__
    // typedef struct sciintmat {
    //   int m,n;
    //   int it ;
    //   int l;
    //   void *D;
    // } SciIntMat ;
    // #endif



    typedef struct {
      char extMethodName[6];
      extend_method extMethod;
    } extension_identity;


    typedef struct {
      int   errorNumber;
      char  message[150];
    } str_error_notification;




    /*********************************************
    * Structures CWT
    ********************************************/

    typedef void(*WScaleFunc)(double *x, int sigInLength, double *psi, int sigOutLength, double ys);


    typedef struct {
      char wname[20];
      int     realOrComplex;
      int     family;
      int     phipsi;
      double  lb;
      double  ub;
      double cpsi;
      WScaleFunc scalef;
    } cwt_identity;

    typedef struct {
      char wname[20];
      char     realOrComplex[20];
      char     family[20];
    } cwt_family;





    /*********************************************
    * Wavelet Structure Declarations
    ********************************************/

    typedef struct {
      int     length;
      double  *pLowPass;
      double  *pHiPass;
    } swt_wavelet;

    typedef void(*Func)(int member, swt_wavelet *pWaveStruct);

    typedef struct {
      char  wname[20];
      int   rOrB;
      int   family;
      int   member;
      Func  analysis;
      Func  synthesis;
    } wavelet_identity;

    typedef struct {
      char  wname[20];
      char   rOrB[20];
      char   family[20];
    } wavelet_family;


    /*********************************************
    * Global Variable Declaration
    ********************************************/

    //  double LowDecomFilCoef[80];
    //  double LowReconFilCoef[80];
    //  double HiDecomFilCoef[80];
    //  double HiReconFilCoef[80];

    /*********************************************
    * swt Variable Declaration
    ********************************************/


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

  static  cwt_identity ci[] = {
      {"sinus", REAL, SINUS, PSI_ONLY, -0.5, 0.5, 1, sinus},
      {"poisson", REAL, POISSON, PSI_ONLY,-10, 10, 1, poisson},
      {"mexh", REAL, MEXICAN_HAT,PSI_ONLY, -5, 5, 1.0, mexihat},
      {"morl",REAL,  MORLET, PSI_ONLY, -4, 4, 1.0, morlet},
      {"DOG", REAL, DOGAUSS, PSI_ONLY, -5, 5, 0.6455109, DOGauss},
      {"gaus1", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus1},
      {"gaus2", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus2},
      {"gaus3", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus3},
      {"gaus4", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus4},
      {"gaus5", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus5},
      {"gaus6", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus6},
      {"gaus7", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus7},
      {"gaus8", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus8},
      {"cmor",COMPLEX, CMORLET, PSI_ONLY, -8, 8, 1, cmorlet_packet},
      {"shan",  COMPLEX, SHANNON, PSI_ONLY, -20, 20, 1, shanwavf_packet},
      {"fbsp", COMPLEX, FBSP, PSI_ONLY, -20, 20, 1, fbspwavf_packet},
      {"cauchy", COMPLEX, CAUCHY, PSI_ONLY, -5, 5, 1, cauchy_packet},
      {"cgau1", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau1_packet},
      {"cgau2", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau2_packet},
      {"cgau3", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau3_packet},
      {"cgau4", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau4_packet},
      {"cgau5", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau5_packet},
      {"cgau6", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau6_packet},
      {"cgau7", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau7_packet},
      {"cgau8", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau8_packet}
    };
    static int cwtIdentityNum = sizeof(ci)/sizeof(cwt_identity);

    static cwt_family cif[] = {
      {"sinus", "REAL", "SINUS"},
      {"poisson", "REAL", "POISSON"},
      {"mexh", "REAL", "MEXICAN_HAT"},
      {"morl","REAL",  "MORLET"},
      {"DOG", "REAL", "DOGAUSS"},
      {"cmor","COMPLEX", "CMORLET"},
      {"shan",  "COMPLEX", "SHANNON"},
      {"fbsp", "COMPLEX", "FBSP"},
      {"cauchy", "COMPLEX", "CAUCHY"},
      {"gaus", "REAL", "GAUSS"},
      {"cgau", "COMPLEX", "CGAUSS"}
    };
    static int cwtFamilyNum = sizeof(cif)/sizeof(cwt_family);

    #define ISODD(x)        ((x/2.0)== ((int)(x/2)) ? 0 : 1)

    /*********************************************
    * Function Prototype
    ********************************************/

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

    static wavelet_identity wi[] = {
      {"haar",ORTH, HAAR, 0, haar_analysis_initialize , haar_synthesis_initialize},
      {"db1", ORTH, DAUBECHIES, 1, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db2", ORTH, DAUBECHIES, 2, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db3", ORTH, DAUBECHIES, 3, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db4", ORTH, DAUBECHIES, 4, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db5", ORTH, DAUBECHIES, 5, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db6", ORTH, DAUBECHIES, 6, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db7", ORTH, DAUBECHIES, 7, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db8", ORTH, DAUBECHIES, 8, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db9", ORTH, DAUBECHIES, 9, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db10", ORTH, DAUBECHIES, 10, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db11", ORTH, DAUBECHIES, 11, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db12", ORTH, DAUBECHIES, 12, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db13", ORTH, DAUBECHIES, 13, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db14", ORTH, DAUBECHIES, 14, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db15", ORTH, DAUBECHIES, 15, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db16", ORTH, DAUBECHIES, 16, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db17", ORTH, DAUBECHIES, 17, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db18", ORTH, DAUBECHIES, 18, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db19", ORTH, DAUBECHIES, 19, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db20", ORTH, DAUBECHIES, 20, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db21", ORTH, DAUBECHIES, 21, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db22", ORTH, DAUBECHIES, 22, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db23", ORTH, DAUBECHIES, 23, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db24", ORTH, DAUBECHIES, 24, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db25", ORTH, DAUBECHIES, 25, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db26", ORTH, DAUBECHIES, 26, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db27", ORTH, DAUBECHIES, 27, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db28", ORTH, DAUBECHIES, 28, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db29", ORTH, DAUBECHIES, 29, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db30", ORTH, DAUBECHIES, 30, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db31", ORTH, DAUBECHIES, 31, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db32", ORTH, DAUBECHIES, 32, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db33", ORTH, DAUBECHIES, 33, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db34", ORTH, DAUBECHIES, 34, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db35", ORTH, DAUBECHIES, 35, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db36", ORTH, DAUBECHIES, 36, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db37", ORTH, DAUBECHIES, 37, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"db38", ORTH, DAUBECHIES, 38, daubechies_analysis_initialize, daubechies_synthesis_initialize},
      {"coif1", ORTH, COIFLETS, 1, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif2", ORTH, COIFLETS, 2, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif3", ORTH, COIFLETS, 3, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif4", ORTH, COIFLETS, 4, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif5", ORTH, COIFLETS, 5, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif6", ORTH, COIFLETS, 6, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif7", ORTH, COIFLETS, 7, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif8", ORTH, COIFLETS, 8, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif9", ORTH, COIFLETS, 9, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif10", ORTH, COIFLETS, 10, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif11", ORTH, COIFLETS, 11, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif12", ORTH, COIFLETS, 12, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif13", ORTH, COIFLETS, 13, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif14", ORTH, COIFLETS, 14, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif15", ORTH, COIFLETS, 15, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif16", ORTH, COIFLETS, 16, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"coif17", ORTH, COIFLETS, 17, coiflets_analysis_initialize, coiflets_synthesis_initialize},
      {"sym2", ORTH, SYMLETS, 2, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym3", ORTH, SYMLETS, 3, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym4", ORTH, SYMLETS, 4, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym5", ORTH, SYMLETS, 5, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym6", ORTH, SYMLETS, 6, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym7", ORTH, SYMLETS, 7, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym8", ORTH, SYMLETS, 8, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym9", ORTH, SYMLETS, 9, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym10", ORTH, SYMLETS, 10, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym11", ORTH, SYMLETS, 11, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym12", ORTH, SYMLETS, 12, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym13", ORTH, SYMLETS, 13, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym14", ORTH, SYMLETS, 14, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym15", ORTH, SYMLETS, 15, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym16", ORTH, SYMLETS, 16, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym17", ORTH, SYMLETS, 17, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym18", ORTH, SYMLETS, 18, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym19", ORTH, SYMLETS, 19, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"sym20", ORTH, SYMLETS, 20, symlets_analysis_initialize, symlets_synthesis_initialize},
      {"bior1.1", BIORTH, SPLINE_BIORTH, 11, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior1.3", BIORTH,SPLINE_BIORTH, 13, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior1.5", BIORTH,SPLINE_BIORTH, 15, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior2.2", BIORTH,SPLINE_BIORTH, 22, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior2.4", BIORTH,SPLINE_BIORTH, 24, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior2.6", BIORTH,SPLINE_BIORTH, 26, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior2.8", BIORTH,SPLINE_BIORTH, 28, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior3.1", BIORTH,SPLINE_BIORTH, 31, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior3.3", BIORTH,SPLINE_BIORTH, 33, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior3.5", BIORTH,SPLINE_BIORTH, 35, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior3.7", BIORTH,SPLINE_BIORTH, 37, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior3.9", BIORTH,SPLINE_BIORTH, 39, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior4.4", BIORTH,SPLINE_BIORTH, 44, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior5.5", BIORTH,SPLINE_BIORTH, 55, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"bior6.8", BIORTH,SPLINE_BIORTH, 68, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
      {"rbior1.1", BIORTH,SPLINE_RBIORTH, 11, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior1.3", BIORTH,SPLINE_RBIORTH, 13, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior1.5", BIORTH,SPLINE_RBIORTH, 15, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior2.2", BIORTH,SPLINE_RBIORTH, 22, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior2.4", BIORTH,SPLINE_RBIORTH, 24, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior2.6", BIORTH,SPLINE_RBIORTH, 26, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior2.8", BIORTH,SPLINE_RBIORTH, 28, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior3.1", BIORTH,SPLINE_RBIORTH, 31, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior3.3", BIORTH,SPLINE_RBIORTH, 33, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior3.5", BIORTH,SPLINE_RBIORTH, 35, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior3.7", BIORTH,SPLINE_RBIORTH, 37, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior3.9", BIORTH,SPLINE_RBIORTH, 39, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior4.4", BIORTH,SPLINE_RBIORTH, 44, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior5.5", BIORTH,SPLINE_RBIORTH, 55, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"rbior6.8", BIORTH,SPLINE_RBIORTH, 68, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
      {"beylkin", ORTH, BEYLKIN, 0, beylkin_analysis_initialize, beylkin_synthesis_initialize},
      {"vaidyanathan", ORTH, VAIDYANATHAN, 0, vaidyanathan_analysis_initialize, vaidyanathan_synthesis_initialize},
      {"dmey", ORTH, DMEY, 0, dmey_analysis_initialize, dmey_synthesis_initialize},
      {"bath4.0", ORTH, BATHLETS, 40, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.1", ORTH, BATHLETS, 41, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.2", ORTH, BATHLETS, 42, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.3", ORTH, BATHLETS, 43, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.4", ORTH, BATHLETS, 44, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.5", ORTH, BATHLETS, 45, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.6", ORTH, BATHLETS, 46, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.7", ORTH, BATHLETS, 47, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.8", ORTH, BATHLETS, 48, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.9", ORTH, BATHLETS, 49, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.10", ORTH, BATHLETS, 410, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.11", ORTH, BATHLETS, 411, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.12", ORTH, BATHLETS, 412, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.13", ORTH, BATHLETS, 413, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.14", ORTH, BATHLETS, 414, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath4.15", ORTH, BATHLETS, 415, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.0", ORTH, BATHLETS, 60, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.1", ORTH, BATHLETS, 61, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.2", ORTH, BATHLETS, 62, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.3", ORTH, BATHLETS, 63, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.4", ORTH, BATHLETS, 64, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.5", ORTH, BATHLETS, 65, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.6", ORTH, BATHLETS, 66, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.7", ORTH, BATHLETS, 67, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.8", ORTH, BATHLETS, 68, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.9", ORTH, BATHLETS, 69, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.10", ORTH, BATHLETS, 610, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.11", ORTH, BATHLETS, 611, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.12", ORTH, BATHLETS, 612, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.13", ORTH, BATHLETS, 613, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.14", ORTH, BATHLETS, 614, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"bath6.15", ORTH, BATHLETS, 615, bathlets_analysis_initialize, bathlets_synthesis_initialize},
      {"legd1", ORTH, LEGENDRE, 1, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd2", ORTH, LEGENDRE, 2, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd3", ORTH, LEGENDRE, 3, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd4", ORTH, LEGENDRE, 4, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd5", ORTH, LEGENDRE, 5, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd6", ORTH, LEGENDRE, 6, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd7", ORTH, LEGENDRE, 7, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd8", ORTH, LEGENDRE, 8, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"legd9", ORTH, LEGENDRE, 9, legendre_analysis_initialize, legendre_synthesis_initialize},
      {"fa1", ORTH, FARRAS, 1, farras_analysis_initialize, farras_synthesis_initialize},
      {"fa2", ORTH, FARRAS, 2, farras_analysis_initialize, farras_synthesis_initialize},
      {"ksq1", ORTH, KINGSBURYQ, 1, kingsburyq_analysis_initialize, kingsburyq_synthesis_initialize},
      {"ksq2", ORTH, KINGSBURYQ, 2, kingsburyq_analysis_initialize, kingsburyq_synthesis_initialize}
    };

    static int waveletIdentityNum = sizeof(wi)/sizeof(wavelet_identity);

    static wavelet_family wif[] = {
      {"haar","ORTH", "HAAR"},
      {"db", "ORTH", "DAUBECHIES"},
      {"coif", "ORTH", "COIFLETS"},
      {"sym", "ORTH", "SYMLETS"},
      {"bior", "BIORTH", "SPLINE_BIORTH"},
      {"beylkin", "ORTH", "BEYLKIN"},
      {"vaidyanathan", "ORTH", "VAIDYANATHAN"},
      {"dmey", "ORTH", "DMEY"},
      {"bath", "ORTH", "BATHLETS"},
      {"legd", "ORTH", "LEGENDRE"},
      {"rbior", "BIORTH","SPLINE_RBIORTH"},
      {"fa", "ORTH", "FARRAS"},
      {"ksq", "ORTH", "KINGSBURYQ"}
    };
    static int waveletFamilyIdentityNum = sizeof(wif)/sizeof(wavelet_family);

    static extend_method dwtMode = SYMH;




    // wavelet_identity wi[];

    static double LowDecomFilCoef[80] = {0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0};

      static double LowReconFilCoef[80] = {0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0};

        static double HiDecomFilCoef[80] = {0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0};

          static double HiReconFilCoef[80] = {0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0};

            static extension_identity ei[] = {
              {"zpd", ZPD}, {"symh", SYMH}, {"symw", SYMW},
              {"asymh", ASYMH}, {"asymw", ASYMW}, {"sp0", SP0},
              {"sp1", SP1}, {"ppd", PPD}, {"per", PER},
              {"spd", SP1}, {"sym", SYMH}, {"asym", ASYMH},
              {"ZPD", ZPD}, {"SYMH", SYMH}, {"SYMW", SYMW},
              {"ASYMH", ASYMH}, {"ASYMW", ASYMW}, {"SP0", SP0},
              {"SP1", SP1}, {"PPD", PPD}, {"PER", PER},
              {"SPD", SP1},  {"SYM", SYMH}, {"ASYM", ASYMH}
            };

            static int extensionIdentityNum = sizeof(ei)/sizeof(extension_identity);


            void
            cowavedec (double *sigIn, int sigInLength, double *sigOutR,
            double *sigOutI, int sigOutLength,
            double *lowDTree1S1, double *hiDTree1S1,
            double *lowDTree2S1, double *hiDTree2S1,
            double *lowDTree1S2, double *hiDTree1S2,
            double *lowDTree2S2, double *hiDTree2S2,
            int filterLen, int *waveDecLengthArray,
            int lengthArrayLengh, int stride, extend_method extMethod);


            void
            cowaverec (double *sigInR, double *sigInI, int sigInLength,
            double *sigOut, int sigOutLength,
            double *lowRTree1S1, double *hiRTree1S1,
            double *lowRTree2S1, double *hiRTree2S1,
            double *lowRTree1S2, double *hiRTree1S2,
            double *lowRTree2S2, double *hiRTree2S2,
            int filterLen, int *waveDecLengthArray,
            int lengthArraylength, int stride,
            extend_method extMethod);

             void
            cowavedec2 (double *matrixIn, int matrixInRow, int matrixInCol,
            double *lowDTree1S1, double *hiDTree1S1,
            double *lowDTree1S2, double *hiDTree1S2,
            int filterLen, int *pLen, double *coef,
            int sigOutLength, int stride, extend_method extMethod);

             void
            cowavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
            double *lowDTree1S1R, double *hiDTree1S1R,
            double *lowDTree1S1C, double *hiDTree1S1C,
            double *lowDTree1S2R, double *hiDTree1S2R,
            double *lowDTree1S2C, double *hiDTree1S2C,
            int filterLen, int *pLen, double *coef,
            int sigOutLength, int stride, extend_method extMethod);

             void
            cowaverec2 (double *coef, int sigInLength,
            double *lowRTree1S1, double *hiRTree1S1,
            double *lowRTree1S2, double *hiRTree1S2,
            int filterLen, double *matrixOut, int matrixOutRow,
            int matrixOutCol, int *pLen, int stride,
            extend_method extMethod);

             void
            cowaverec2a (double *coef, int sigInLength,
            double *lowRTree1S1R, double *hiRTree1S1R,
            double *lowRTree1S1C, double *hiRTree1S1C,
            double *lowRTree1S2R, double *hiRTree1S2R,
            double *lowRTree1S2C, double *hiRTree1S2C,
            int filterLen, double *matrixOut, int matrixOutRow,
            int matrixOutCol, int *pLen, int stride,
            extend_method extMethod);

             void
            copmd (double *matrixInR, double *matrixInI, int sigInLength,
            int InRow, int InCol, double *matrixOutR, double *matrixOutI);

             void
            copmr (double *matrixInR, double *matrixInI, int sigInLength,
            int InRow, int InCol, double *matrixOutR, double *matrixOutI);

            /*------------------------------------------*/
            /* Wavelet Family Function                  */
            /* -----------------------------------------*/
             void filter_clear ();
             void orth_filt_group (double *filterIn, int sigInLength,
            double *filterLowRec,
            double *filterLowDec,
            double *filterHiRec,
            double *filterHiDec);
             void bior_filt_group (double *f1, int sigInLength1,
            double *f2, int sigInLength2,
            double *lowDecom, int sigOutLength1,
            double *hiDecom, int sigOutLength2,
            double *lowRecon, int sigOutLength3,
            double *hiRecon, int sigOutLength4);
             void wavelet_parser (char *wname, int *family, int *member);
             void wavelet_fun_parser (char *wname, int *ii);
             void wave_len_validate (int sigInLen, int waveLength, int *lev, int *val);

             extend_method getdwtMode();
	     void setdwtMode(extend_method mode);
             void dwt_write (char *mode, int *errCode);
	     extend_method char_to_extend_method(char *mode);
             void dwt_parse(char **strr);
             void dwt (double *sigIn, int sigInLength, double *lowDe,
            double *hiDe, int filterLen, double *approx,
            double *detail, int sigOutLength,
            extend_method extMethod);
             void dwt_neo (double *sigIn, int sigInLength, double *lowDe,
            double *hiDe, int filterLen, double *approx,
            double *detail, int sigOutLength,
            extend_method extMethod);
             void dwt_nex (double *sigIn, int sigInLength, double *lowDe,
            double *hiDe, int filterLen, double *approx,
            double *detail, int sigOutLength);
             void dwt_no_extension (double *sigIn, int sigInLength, double *lowDe,
            double *hiDe, int filterLen, double *approx,
            double *detail, int sigOutLength);
             void dwt_conv (double *sigIn, int sigInLength, double *lowDe,
            double *hiDe, int filterLen, double *approx,
            double *detail, int sigOutLength);
             void idwt_complete (double *approx, double *detail,
            int sigInLength, double *lowRe,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength);
             void idwt_neo (double *approx, double *detail,
            int sigInLength, double *lowRe,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength);
             void idwt_complete_ex (double *approx, double *detail,
            int sigInLength, double *lowRe,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength,
            extend_method extMethod);
             void idwt_approx (double *approx, int sigInLength,
            double *lowRe, int filterLen,
            double *sigOut, int sigOutLength);
             void idwt_approx_ex (double *approx, int sigInLength,
            double *lowRe, int filterLen,
            double *sigOut, int sigOutLength,
            extend_method extMethod);
             void idwt_approx_neo (double *approx, int sigInLength,
            double *lowRe, int filterLen,
            double *sigOut, int sigOutLength);
             void idwt_detail (double *detail, int sigInLength,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength);
             void idwt_detail_ex (double *detail, int sigInLength,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength,
            extend_method extMethod);
             void idwt_detail_neo (double *detail, int sigInLength,
            double *hiRe, int filterLen,
            double *sigOut, int sigOutLength);
             void wave_dec_len_cal (int filterLen, int sigLength,
            int stride, int *waveDecLengthArray);
             void wavedec (double *sigIn, int sigInLength, double *sigOut,
            int sigOutLength, double *lowDe, double *hiDe,
            int filterLen, int *waveDecLengthArray,
            int lengthArrayLengh, int stride,
            extend_method extMethod);
             void waverec (double *sigIn, int sigInLength, double *sigOut,
            int sigOutLength, double *lowRe, double *hiRe,
            int filterLen, int *waveDecLengthArray,
            int lengthArraylength, int stride,
            extend_method extMethod);
             void wenergy (double *coef, int coefLen, int *lenArray,
            int arrayLen, double *aE, int aELen,
            double *dE, int dELen);
             void detcoef (double *sigIn, int sigInLength,
            int *waveDecLengthArray, int arrayLen,
            double *sigOut, int sigOutLength,
            int stride, int level);
             void appcoef (double *sigIn, int sigInLength, double *sigOut,
            int sigOutLength, double *lowRe, double *hiRe,
            int filterLen, int *waveDecLengthArray,
            int lengthArraylength, int stride, int level,
            extend_method extMethod);
             void wrcoef (double *sigIn, int sigInLength, double *lowRe,
            double *hiRe, int filterLen,
            int *waveDecLengthArray, int arrayLen,
            double *sigOut, int sigOutLength,
            char *coefType, int stride, int level,
            extend_method extMethod);
             void upcoef_len_cal (int sigInLength, int filterLen,
            int stride, int *sigOutLength,
            int *sigOutLengthDefault);
             void upwlev (double *coefArray, int coefLen,
            int *waveDecLengthArray,	int arrayLen,
            double *lowRe, double *hiRe, int filterLen,
            double *newCoefArray, int newCoefLen,
            int *newLenArray, int newArrayLen,
            double *approx, int approxLen, int stride,
            extend_method extMethod);
             void upcoef (double *sigIn, int sigInLength, double *lowRe,
            double *hiRe, int filterLen, double *sigOut,
            int sigOutLength, int defaultLength,
            char *coefType, int step);

             void dwt2D (double *matrixIn, int matrixInRow,
            int matrixInCol, double *matrixOutApprox,
            double *matrixOutColDetail,
            double *matrixOutRowDetail,
            double *matrixOutDetail, int matrixOutRow,
            int matrixOutCol, double *lowDe, double *hiDe,
            int filterLen, extend_method extMethod);
             void
            dwt2D_neo_a (double *matrixIn, int matrixInRow, int matrixInCol,
            double *matrixOutApprox, double *matrixOutColDetail,
            double *matrixOutRowDetail, double *matrixOutDetail,
            int matrixOutRow, int matrixOutCol, double *lowDeR,
            double *hiDeR, double *lowDeC, double *hiDeC,
            int filterLen, extend_method extMethod);
             void dwt2D_neo (double *matrixIn, int matrixInRow,
            int matrixInCol, double *matrixOutApprox,
            double *matrixOutColDetail,
            double *matrixOutRowDetail,
            double *matrixOutDetail, int matrixOutRow,
            int matrixOutCol, double *lowDe, double *hiDe,
            int filterLen, extend_method extMethod);
             void idwt2D (double *matrixInApprox,
            double *matrixInColDetail,
            double *matrixInRowDetail,
            double *matrixInDetail,
            int matrixInRow, int matrixInCol, double *lowRe,
            double *hiRe, int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol,
            extend_method extMethod);
             void idwt2D_neo (double *matrixInApprox, double *matrixInColDetail,
            double *matrixInRowDetail, double *matrixInDetail,
            int matrixInRow, int matrixInCol, double *lowRe,
            double *hiRe, int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol);
             void
            idwt2D_neo_a (double *matrixInApprox, double *matrixInColDetail,
            double *matrixInRowDetail, double *matrixInDetail,
            int matrixInRow, int matrixInCol, double *lowReR,
            double *hiReR, double *lowReC, double *hiReC,
            int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol);

             void wave_mem_cal (int *pLen, int stride, int *total);
             void matrix_wavedec_len_cal (int matrixInRow, int matrixInCol,
            int stride, int filterLen,
            int *pLen);
             void matrix_locate (int stride, int *pLen, int *pH,
            int *pV, int *pD);
             void wavedec2 (double *matrixIn, int matrixInRow,
            int matrixInCol, double *lowDe, double *hiDe,
            int filterLen, int *pLen, double *coef,
            int sigOutLength, int stride,
            extend_method extMethod);
             void
            wavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
            double *lowDeR, double *hiDeR, double *lowDeC,
            double *hiDeC, int filterLen, int *pLen,
            double *coef, int sigOutLength, int stride,
            extend_method extMethod);
             void waverec2 (double *coef, int sigInLength, double *lowRe,
            double *hiRe, int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol, int *pLen,
            int stride, extend_method extMethod);
             void
            waverec2a (double *coef, int sigInLength, double *lowReR,
            double *hiReR, double *lowReC, double *hiReC,
            int filterLen, double *matrixOut, int matrixOutRow,
            int matrixOutCol, int *pLen, int stride,
            extend_method extMethod);

             void wenergy_2output (double *coef, int sigInLength,
            int *pLen, double *ae, double *de,
            int deLength, int stride);
             void wenergy_4output (double *coef, int sigInLength,
            int *pLen, double *ae, double *he,
            double *ve, double *de, int deLength,
            int stride);
             void detcoef2 (double *coef, int sigInLength, double *coefOut,
            int sigOutLength, int *pLen, int stride,
            int level, char *coefType);
             void appcoef2 (double *coef, int sigInLength, double *lowRe,
            double *hiRe, int filterLen, double *coefOut,
            int matrixOutRow, int matrixOutCol, int *pLen,
            int stride, int level, extend_method extMethod);
             void wrcoef2 (double *coef, int sigInLength, double *lowRe,
            double *hiRe, int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol, int *pLen,
            int stride, int level, char *type,
            extend_method extMethod);
             void upwlev2 (double *coef, int sigInLength, double *lowRe,
            double *hiRe,
            int filterLen, int *pLen, int matrixRow, int matrixCol,
            double *approx, int approxLen, double *newCoef,
            int newCoefLen, int *newLenMatrix, int lenMatrixRow,
            int lenMatrixCol, int stride, extend_method extMethod);
             void upcoef2 (double *matrixIn, int matrixInRow,
            int matrixInCol, double *lowRe, double *hiRe,
            int filterLen, double *matrixOut,
            int matrixOutRow, int matrixOutCol,
            int matrixOutDefaultRow,
            int matrixOutDefaultCol,
            int step, char *type);//, extend_method extMethod);

             void dwt3d_tran(double *mat3DIn, int row1, int col1, int sli1,
            double *mat3DOut, int row2, int col2, int sli2);

             void dwt3d_line_forward(double *mat3DIn, int row1, int col1, int sli1,
            double *mat3DOutApp, double *mat3DOutDet,
            int row2, int col2, int sli2,
            double *loDe, double *hiDe, int filterLen,
            extend_method extMethod);

             void dwt3d_tran_z(double *mat3DIn, int row1, int col1, int sli1,
            double *mat3DOut, int row2, int col2, int sli2);

             void dwt3d_tran_z_inv(double *mat3DIn, int row1, int col1, int sli1,
            double *mat3DOut, int row2, int col2, int sli2);

             void dwt3d_combine(double *mat1, double *mat2, double *mat3,
            double *mat4, double *mat5, double *mat6,
            double *mat7, double *mat8, int rowIn,
            int colIn, int sliIn, double *matOut,
            int rowOut, int colOut, int sliOut);
             void dwt3d_line_reverse(double *mat3DInApp, double *mat3DInDet,
            int row1, int col1, int sli1,
            double *mat3DOut, int row2, int col2,
            int sli2,
            double *loDe, double *hiDe, int filterLen);


             void dwt3d_split(double *matIn, int rowIn, int colIn, int sliIn,
            double *mat1, double *mat2, double *mat3,
            double *mat4, double *mat5, double *mat6,
            double *mat7, double *mat8, int rowOut,
            int colOut, int sliOut);
             void dwt3(double *mat3DIn, int row, int col, int sli,
            double *mat3DOut, int row2, int col2, int sli2,
            int r, int c, int s, double *Lo1, double *Hi1,
            double *Lo2, double *Hi2, double *Lo3, double *Hi3,
            int fLen1, int fLen2, int fLen3, extend_method extMethod);







             void cwt_fun_parser(char *wname, int *ind);
             void cwt_len_cal (int sigInLength, int scale, int *sigOutLength, double *delta);
             void full_range_scalef (char *wname, double *f, int sigOutLength);
             void scale_real (double *f, int sigInLength, double delta, double *fout, int sigOutLength);
            // void scale_complex (double *f, int sigInLength, double delta, double *fout, int sigOutLength);
             void cwt_conv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
             void cwt_iconv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
             void cwt_conv_complex (double *sigIn, int sigInLength, double *fr, double *fi, int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);
             void cwt_conv_complex_complex (double *a, double *b, int sigInLength,double *c, double *d,  int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);




             /*********************************************
             * Function Prototype
             ********************************************/
              void swt_conv(double *sigIn, int sigInLength,
             double *approx, int approxLength,
             double *detail, int detailLength,
             double *filterLow, double *filterHi,
             int filterLength);
              void  swt_out1 (double *sigIn, int sigInLength,
             double *sigOutMatrix, int rowLength,
             int colLength, double *filterLow,
             double *filterHi, int filterLength, int step);

              void swt_out2 (double *sigIn, int sigInLength,
             double *approxMatrix, double *detailMatrix,
             int rowLength, int colLength, double *filterLow,
             double *filterHi, int filterLength, int step);

              void iswt_conv (double *approx, double *detail, int sigInLength,
             double *sigOut, int sigOutLength, double *filterLow,
             double *filterHi, int filterLength);
              void iswt_conv_step (double *approx, double *detail, int sigInLength,
             double *sigOut, int sigOutLength, double *filterLow,
             double *filterHi, int filterLength, int level);

              void iswt_input1 (double *matrixIn, int rowLength, int colLength,
             double *sigOut, int sigOutLength, double *filterLow,
             double *filterHi, int filterLength);

              void iswt_input2 (double *matrixApproxIn, double *matrixDetailIn,
             int rowLength, int colLength,
             double *sigOut, int sigOutLength, double *filterLow,
             double *filterHi, int filterLength);
              void swt2_output4(double *matrixIn, int matrixInRow, int matrixInCol,
             double *matrixOutApprox, double *matrixOutColDetail,
             double *matrixOutRowDetail, double *matrixOutDetail,
             int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);
              void swt2_output4_step(double *matrixIn, int matrixInRow, int matrixInCol,
             double *matrixOutApprox, double *matrixOutColDetail,
             double *matrixOutRowDetail, double *matrixOutDetail,
             int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);
              void swt2_output1_step(double *matrixIn, int matrixInRow,
             int matrixInCol,  double *matrixOut,
             int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);
              void iswt2(double *matrixInApprox, double *matrixInColDetail,
             double *matrixInRowDetail, double *matrixInDetail,
             int matrixInRow, int matrixInCol,
             double *matrixOut, int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);
              void iswt2_input4_step(double *matrixInApprox, double *matrixInColDetail,
             double *matrixInRowDetail, double *matrixInDetail,
             int matrixInRow, int matrixInCol,
             double *matrixOut, int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);
              void iswt2_input1_step(double *matrixIn,  int matrixInRow, int matrixInCol,
             double *matrixOut, int matrixOutRow, int matrixOutCol,
             double *filterLow, double *filterHi,
             int filterLength, int step);

             /*------------------------------------------*/
             /* Utility Function                         */
             /* -----------------------------------------*/
              void matrix_tran (double *matrixIn, int matrixInRow,
             			 int matrixInCol, double *matrixOut,
             			 int matrixOutRow, int matrixOutCol);
              void wrev (const double *sigIn, int sigInLength,
             		  double *sigOut, int sigOutLength);
              void qmf_even (const double *sigIn, int sigInLength,
             		      double *sigOut, int sigOutLength);
              void qmf_odd (double *sigIn, int sigInLength,
             		     double *sigOut, int sigOutLength);
              void qmf_wrev (const double *sigIn, int sigInLength,
             		      double *sigOut, int sigOutLength);
              void verbatim_copy (const double *sigIn, int sigInLength,
             			   double *sigOut, int sigOutLength);
              void dyaddown_1D_keep_odd (double *sigIn, int sigInLength,
             			       double *sigOut, int sigOutLength);
              void dyaddown_1D_keep_even (double *sigIn, int sigInLength,
             				double *sigOut, int sigOutLength);
              void dyaddown_2D_keep_odd_row (double *matrixIn,
             				      int matrixInRow,
             				      int matrixInCol,
             				      double *matrixOut,
             				      int matrixOutRow,
             				      int matrixOutCol);
              void dyaddown_2D_keep_odd_col (double *matrixIn,
             				      int matrixInRow,
             				      int matrixInCol,
             				      double *matrixOut,
             				      int matrixOutRow,
             				      int matrixOutCol);
              void dyaddown_2D_keep_even_row (double *matrixIn,
             				       int matrixInRow,
             				       int matrixInCol,
             				       double *matrixOut,
             				       int matrixOutRow,
             				       int matrixOutCol);
              void dyaddown_2D_keep_even_col (double *matrixIn,
             				       int matrixInRow,
             				       int matrixInCol,
             				       double *matrixOut,
             				       int matrixOutRow,
             				       int matrixOutCol);
              void dyaddown_2D_keep_odd (double *matrixIn,
             				  int matrixInRow,
             				  int matrixInCol,
             				  double *matrixOut,
             				  int matrixOutRow,
             				  int matrixOutCol);
              void dyaddown_2D_keep_even (double *matrixIn,
             				   int matrixInRow,
             				   int matrixInCol,
             				   double *matrixOut,
             				   int matrixOutRow,
             				   int matrixOutCol);
              void dyadup_1D_feed_odd (double *sigIn, int sigInLength,
             				double *sigOut, int sigOutLength);
              void dyadup_1D_feed_even (double *sigIn, int sigInLength,
             				 double *sigOut, int sigOutLength);
              void dyadup_2D_feed_odd_row (double *matrixIn,
             				    int matrixInRow,
             				    int matrixInCol,
             				    double *matrixOut,
             				    int matrixOutRow,
             				    int matrixOutCol);
              void dyadup_2D_feed_odd_col (double *matrixIn,
             				    int matrixInRow,
             				    int matrixInCol,
             				    double *matrixOut,
             				    int matrixOutRow,
             				    int matrixOutCol);
              void dyadup_2D_feed_even_row (double *matrixIn,
             				     int matrixInRow,
             				     int matrixInCol,
             				     double *matrixOut,
             				     int matrixOutRow,
             				     int matrixOutCol);
              void dyadup_2D_feed_even_col (double *matrixIn,
             				     int matrixInRow,
             				     int matrixInCol,
             				     double *matrixOut,
             				     int matrixOutRow,
             				     int matrixOutCol);
              void dyadup_2D_feed_odd (double *matrixIn,
             				int matrixInRow,
             				int matrixInCol,
             				double *matrixOut,
             				int matrixOutRow,
             				int matrixOutCol);
              void dyadup_2D_feed_even (double *matrixIn,
             				 int matrixInRow,
             				 int matrixInCol,
             				 double *matrixOut,
             				 int matrixOutRow,
             				 int matrixOutCol);
              void extend_method_parse (char *mode, extend_method *extMethod);
              void wextend_1D_center (double *sigIn, int sigInLength,
             			       double *sigOut, int sigOutLength,
             			       extend_method method);
              void wextend_1D_left (double *sigIn, int sigInLength,
             			       double *sigOut, int sigOutLength,
             			       extend_method method);
              void wextend_1D_right (double *sigIn, int sigInLength,
             			       double *sigOut, int sigOutLength,
             			       extend_method method);
              void wextend_2D (double *matrixIn, int matrixInRow,
             			int matrixInCol, double *matrixOut,
             			int matrixOutRow, int matrixOutCol,
             			extend_method extMethod, char *rowOpt,
             			char *colOpt);
              void wextend_2D_row (double *matrixIn, int matrixInRow,
             			    int matrixInCol, double *matrixOut,
             			    int matrixOutRow, int matrixOutCol,
             			    extend_method extMethod, char *Opt);
              void wextend_2D_col (double *matrixIn, int matrixInRow,
             			    int matrixInCol, double *matrixOut,
             			    int matrixOutRow, int matrixOutCol,
             			    extend_method extMethod, char *Opt);
              void wkeep_1D_center (double *sigIn, int sigInLength,
             			     double *sigOut, int sigOutLength);
              void wkeep_1D_left (double *sigIn, int sigInLength,
             			   double *sigOut, int sigOutLength);
              void wkeep_1D_right (double *sigIn, int sigInLength,
             			    double *sigOut, int sigOutLength);
              void wkeep_1D_index (double *sigIn, int sigInLength,
             			    double *sigOut, int sigOutLength,
             			    int first);
              void wkeep_2D_center (double *matrixIn, int matrixInRow,
             			     int matrixInCol, double *matrixOut,
             			     int matrixOutRow, int matrixOutCol);
              void wkeep_2D_index (double *matrixIn, int matrixInRow,
             			    int matrixInCol, double *matrixOut,
             			    int matrixOutRow, int matrixOutCol,
             			    int rowFirst, int colFirst);
              void conv (double *sigIn, int sigInLength,
             		  double *sigOut, int sigOutLength,
             		  double *fiter, int filterLength);
              void i_conv (double *sigIn, int sigInLength,
             		  double *sigOut, int sigOutLength,
             		  double *fiter, int filterLength);
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
               void wcodemat_matrix (double *matrixIn, int matrixInRow, int matrixInCol,
              	                   double *matrixOut, int matrixOutRow, int matrixOutCol,
              					   int minv, int maxv, int abso);
               void wcodemat_matrix_col (double *matrixIn, int matrixInRow, int matrixInCol,
              	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
              						   int minv, int maxv, int abso);
               void wcodemat_matrix_row (double *matrixIn, int matrixInRow, int matrixInCol,
              	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
              						   int minv, int maxv, int abso);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /** _SWTLIB_H_   **/
