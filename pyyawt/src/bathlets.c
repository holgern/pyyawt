/*
 * -------------------------------------------------------------------------
 * bathlets.c -- Bathlets wavelets coefficents.
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2007  Roger Liu
 * Copyright (C) 20010-2012  Holger Nahrstaedt
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

#include "swtlib.h"

/*********************************************
 * Local Variable (Filter Coefficent)
 ********************************************/



static const double bath4_0[4] = {
 	           0.48296291, 0.83651630,
		       0.22414387, -0.12940952};

static const double bath4_1[4] = {
			   0.49123165, 0.83422448,
			   0.21587513, -0.12711776};

static const double bath4_2[4] = {
			   0.50021916, 0.83155877,
			   0.20688763, -0.12445201};

static const double bath4_3[4] = {
	           0.50981884, 0.82850718,
               0.19728794, -0.12100415};

static const double bath4_4[4] = {
 	           0.51995636, 0.82505202,
			   0.18715315, -0.11794524};

static const double bath4_5[4] = {
			   0.53035917, 0.82124957,
               0.17674761, -0.11414279};

static const double bath4_6[4] = {
			   0.54081214, 0.81716331,
               0.16629464, -0.11005653};

static const double bath4_7[4] = {
	           0.55103446, 0.81290209,
               0.15607232, -0.10579531};

static const double bath4_8[4] = {
 	           0.56081569, 0.80857244,
			   0.14629108, -0.10146566};

static const double bath4_9[4] = {
			   0.57003174, 0.80426075,
               0.13707504, -0.09715397};

static const double bath4_10[4] = {
			   0.57860227, 0.80004303,
               0.12850451, -0.09293625};

static const double bath4_11[4] = {
	           0.58648244, 0.79598311,
               0.12062434, -0.08887633};

static const double bath4_12[4] = {
 	           0.59368629, 0.79211489,
               0.11342048, -0.08500810};

static const double bath4_13[4] = {
			   0.60026062, 0.78845055,
               0.10684615, -0.08134377};

static const double bath4_14[4] = {
			   0.60623768, 0.78500518,
               0.10086910, -0.07789882};

static const double bath4_15[4] = {
	           0.61167674, 0.78177335,
               0.09543005, -0.07466658};

static const double bath6_0[6] = {
 	           0.44841478, 0.77552168,
               0.39119223, -0.14502791,
               -0.13250022, 0.07661301};

static const double bath6_1[6] = {
			   0.45893430, 0.79242211,
               0.35372518, -0.14644658,
               -0.10555270, 0.06113125};

static const double bath6_2[6] = {
			   0.47502467, 0.80303539,
               0.31496201, -0.14495509,
               -0.08287990, 0.04902648};

static const double bath6_3[6] = {
	           0.49303553, 0.80842688,
               0.27886696, -0.14083708,
               -0.06479571, 0.03951698};

static const double bath6_4[6] = {
 	           0.51065493, 0.81006904,
               0.24732487, -0.13503181,
               -0.05087302, 0.03206956};

static const double bath6_5[6] = {
			   -0.08338125, 0.30427515,
               0.84925993, 0.41893704,
               -0.05877191, -0.01610541};

static const double bath6_6[6] = {
			   -0.08279787, 0.30367125,
               0.84917267, 0.41959533,
               -0.05926801, -0.01615980};

static const double bath6_7[6] = {
	           -0.08233775, 0.30257581,
               0.84901846, 0.42074239,
               -0.05957393, -0.01621142};

static const double bath6_8[6] = {
 	           -0.08195956, 0.30100677,
               0.84879593, 0.42236347,
               -0.05972960, -0.01626346};

static const double bath6_9[6] = {
			   -0.08163489, 0.29893998,
               0.84849656, 0.42448474,
               -0.05975489, -0.01631794};

static const double bath6_10[6] = {
			   -0.08134182, 0.29631466,
               0.84810454, 0.42716837,
               -0.05965593, -0.01637625};

static const double bath6_11[6] = {
	           -0.08106069, 0.29303272,
               0.84759506, 0.43051332,
               -0.05942758, -0.01643926};

static const double bath6_12[6] = {
 	           -0.08077071, 0.28895499,
               0.84693143, 0.43465896,
               -0.05905393, -0.01650717};

static const double bath6_13[6] = {
			   -0.08044643, 0.28389820,
               0.84606093, 0.43978754,
               -0.05850772, -0.01657896};

static const double bath6_14[6] = {
			   -0.08005385, 0.27764083,
               0.84491061, 0.44611735,
               -0.05774997, -0.01665140};

static const double bath6_15[6] = {
	           -0.07954693, 0.26995392,
               0.84338655, 0.45387025,
               -0.05673284, -0.01671738};

/*********************************************
 * Global Function
 ********************************************/

void
bathlets_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{

//   double *pFilterCoef;



  switch (member)
    {
    case 40:
//       pFilterCoef = bath4_0;
	  pWaveStruct->length = 4;
        wrev(bath4_0, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_0, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 41:
//       pFilterCoef = bath4_1;
	  pWaveStruct->length = 4;
        wrev(bath4_1, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_1, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 42:
//       pFilterCoef = bath4_2;
	  pWaveStruct->length = 4;
        wrev(bath4_2, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_2, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 43:
//       pFilterCoef = bath4_3;
	  pWaveStruct->length = 4;
        wrev(bath4_3, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_3, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 44:
//       pFilterCoef = bath4_4;
	  pWaveStruct->length = 4;
        wrev(bath4_4, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_4, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 45:
//       pFilterCoef = bath4_5;
	  pWaveStruct->length = 4;
        wrev(bath4_5, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_5, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 46:
//       pFilterCoef = bath4_6;
	  pWaveStruct->length = 4;
        wrev(bath4_6, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_6, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 47:
//       pFilterCoef = bath4_7;
	  pWaveStruct->length = 4;
        wrev(bath4_7, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_7, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 48:
//       pFilterCoef = bath4_8;
	  pWaveStruct->length = 4;
        wrev(bath4_8, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_8, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 49:
//       pFilterCoef = bath4_9;
	  pWaveStruct->length = 4;
        wrev(bath4_9, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_9, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 410:
//       pFilterCoef = bath4_10;
	  pWaveStruct->length = 4;
        wrev(bath4_10, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_10, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 411:
//       pFilterCoef = bath4_11;
	  pWaveStruct->length = 4;
        wrev(bath4_11, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_11, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 412:
//       pFilterCoef = bath4_12;
	  pWaveStruct->length = 4;
        wrev(bath4_12, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_12, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 413:
//       pFilterCoef = bath4_13;
	  pWaveStruct->length = 4;
        wrev(bath4_13, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_13, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 414:
//       pFilterCoef = bath4_14;
	  pWaveStruct->length = 4;
        wrev(bath4_14, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_14, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 415:
//       pFilterCoef = bath4_15;
	  pWaveStruct->length = 4;
        wrev(bath4_15, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath4_15, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 60:
//       pFilterCoef = bath6_0;
	  pWaveStruct->length = 6;
        wrev(bath6_0, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_0, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 61:
//       pFilterCoef = bath6_1;
	  pWaveStruct->length = 6;
        wrev(bath6_1, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_1, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 62:
//       pFilterCoef = bath6_2;
	  pWaveStruct->length = 6;
        wrev(bath6_2, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_2, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 63:
//       pFilterCoef = bath6_3;
	  pWaveStruct->length = 6;
        wrev(bath6_3, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_3, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 64:
//       pFilterCoef = bath6_4;
	  pWaveStruct->length = 6;
        wrev(bath6_4, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_4, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 65:
//       pFilterCoef = bath6_5;
	  pWaveStruct->length = 6;
      wrev(bath6_5, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_5, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 66:
//       pFilterCoef = bath6_6;
	  pWaveStruct->length = 6;
        wrev(bath6_6, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_6, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 67:
//       pFilterCoef = bath6_7;
	  pWaveStruct->length = 6;
        wrev(bath6_7, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_7, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 68:
//       pFilterCoef = bath6_8;
	  pWaveStruct->length = 6;
        wrev(bath6_8, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_8, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 69:
//       pFilterCoef = bath6_9;
	  pWaveStruct->length = 6;
        wrev(bath6_9, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_9, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 610:
//       pFilterCoef = bath6_10;
	  pWaveStruct->length = 6;
        wrev(bath6_10, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_10, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 611:
//       pFilterCoef = bath6_11;
	  pWaveStruct->length = 6;
        wrev(bath6_11, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_11, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 612:
//       pFilterCoef = bath6_12;
	  pWaveStruct->length = 6;
        wrev(bath6_12, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_12, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 613:
//       pFilterCoef = bath6_13;
	  pWaveStruct->length = 6;
        wrev(bath6_13, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_13, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 614:
//       pFilterCoef = bath6_13;
	  pWaveStruct->length = 6;
        wrev(bath6_13, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_13, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 615:
//       pFilterCoef = bath6_15;
	  pWaveStruct->length = 6;
        wrev(bath6_15, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(bath6_15, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    default:
      printf("db%d is not available!\n",member);
      exit(0);
    }

//   wrev(pFilterCoef, pWaveStruct->length,
//        LowDecomFilCoef, pWaveStruct->length);
//   qmf_wrev(pFilterCoef, pWaveStruct->length,
// 	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
bathlets_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{

  double *pFilterCoef;

  switch (member)
    {
    case 40:
//       pFilterCoef = bath4_0;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_0, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_0, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 41:
//       pFilterCoef = bath4_1;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_1, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_1, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 42:
//       pFilterCoef = bath4_2;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_2, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_2, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 43:
//       pFilterCoef = bath4_3;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_3, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_3, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 44:
//       pFilterCoef = bath4_4;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_4, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_4, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 45:
//       pFilterCoef = bath4_5;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_5, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_5, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 46:
//       pFilterCoef = bath4_6;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_6, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_6, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 47:
//       pFilterCoef = bath4_6;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_6, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_6, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 48:
//       pFilterCoef = bath4_8;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_8, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_8, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 49:
//       pFilterCoef = bath4_9;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_9, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_9, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 410:
//       pFilterCoef = bath4_10;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_10, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_10, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 411:
//       pFilterCoef = bath4_11;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_11, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_11, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 412:
//       pFilterCoef = bath4_12;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_12, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_12, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 413:
//       pFilterCoef = bath4_13;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_13, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_13, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 414:
//       pFilterCoef = bath4_14;
	  pWaveStruct->length = 4;
      verbatim_copy(bath4_14, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_14, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 415:
//       pFilterCoef = bath4_15;
	  pWaveStruct->length = 4;
        verbatim_copy(bath4_15, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath4_15, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 60:
//       pFilterCoef = bath6_0;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_0, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_0, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 61:
//       pFilterCoef = bath6_1;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_1, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_1, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 62:
//       pFilterCoef = bath6_2;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_2, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_2, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 63:
//       pFilterCoef = bath6_3;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_3, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_3, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 64:
//       pFilterCoef = bath6_4;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_4, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_4, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 65:
//       pFilterCoef = bath6_5;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_5, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_5, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 66:
//       pFilterCoef = bath6_6;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_6, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_6, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 67:
//       pFilterCoef = bath6_7;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_7, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_7, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 68:
//       pFilterCoef = bath6_8;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_8, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_8, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 69:
//       pFilterCoef = bath6_9;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_9, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_9, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 610:
//       pFilterCoef = bath6_10;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_10, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_10, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 611:
//       pFilterCoef = bath6_11;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_11, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_11, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 612:
//       pFilterCoef = bath6_12;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_12, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_12, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 613:
//       pFilterCoef = bath6_13;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_13, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_13, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 614:
//       pFilterCoef = bath6_14;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_14, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_14, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    case 615:
//       pFilterCoef = bath6_15;
	  pWaveStruct->length = 6;
        verbatim_copy(bath6_15, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(bath6_15, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
      break;
    default:
      printf("db%d is not available!\n",member);
      exit(0);
    }

//   verbatim_copy(pFilterCoef, pWaveStruct->length,
// 		LowReconFilCoef, pWaveStruct->length);
//   qmf_even(pFilterCoef, pWaveStruct->length,
//       HiReconFilCoef, pWaveStruct->length );
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
