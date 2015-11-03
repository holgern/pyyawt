/*
 * -------------------------------------------------------------------------
 * legendre.c -- Legendre wavelets coefficents.
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

static const double legd1[2] = {
	0.7071068,   0.7071068};

static const double legd2[4] = {
	0.4419417,    0.2651650,    0.2651650,    0.4419417};

static const double legd3[6] = {
    0.3480291,   0.1933495,    0.1657282,    0.1657282,    0.1933495,   0.3480291};

static const double legd4[8] = {
	0.2962391,    0.1595133,    0.1305109,    0.1208434,    0.1208434,
    0.1305109,    0.1595133,    0.2962391};

static const double legd5[10] = {
     0.2622950,    0.1388621,    0.1110897,    0.0996958,    0.0951642,
	 0.0951642,    0.0996958,    0.1110897,    0.1388621,    0.2622950};

static const double legd6[14] = {
    0.2191762,    0.1139717,    0.0891952,    0.0778688,    0.0717213,
	0.0683462,    0.0668274,    0.0668274,    0.0683462,    0.0717213,
    0.0778688,    0.0891952,    0.1139717,    0.2191762};

static const double legd7[16] = {
    0.2043035,    0.1056743,    0.0821911,    0.0712323,    0.0650382,
	0.0613217,    0.0591701,    0.0581756,    0.0581756,    0.0591701,
    0.0613217,    0.0650382,    0.0712323,    0.0821911,    0.1056743,
    0.2043035};

static const double legd8[18] = {
    0.1920980,    0.0989595,    0.0766138,    0.0660464,    0.059931,
	0.0560954,    0.0536565,    0.0521964,    0.0515097,    0.0515097,
    0.0521964,    0.0536565,    0.0560954,    0.059931,    0.0660464,
    0.0766138,    0.0989595,    0.1920980};

static const double legd9[20] = {
    0.1818140,    0.0933640,    0.0720236,    0.0618385,    0.0559826,
	0.0520021,    0.0494341,    0.0477392,    0.0467014,    0.0462072,
    0.0462072,    0.0467014,    0.0477392,    0.0494341,    0.0520021,
    0.0559826,    0.0618385,    0.0720236,    0.0933640,    0.1818140};

/*********************************************
 * Global Function
 ********************************************/

void
legendre_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;
  //double sum;
  //int count;


  switch (member)
    {
    case 1:
//       pFilterCoef = legd1;
	  pWaveStruct->length = 2;
	    wrev(legd1, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd1, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = legd2;
	  pWaveStruct->length = 4;
	    wrev(legd2, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd2, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 3:
//       pFilterCoef = legd3;
	  pWaveStruct->length = 6;
	    wrev(legd3, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd3, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 4:
//       pFilterCoef = legd4;
	  pWaveStruct->length = 8;
	    wrev(legd4, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd4, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 5:
//       pFilterCoef = legd5;
	  pWaveStruct->length = 10;
	    wrev(legd5, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd5, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 6:
//       pFilterCoef = legd6;
	  pWaveStruct->length = 14;
	    wrev(legd6, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd6, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 7:
//       pFilterCoef = legd7;
	  pWaveStruct->length = 16;
	    wrev(legd7, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd7, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 8:
//       pFilterCoef = legd8;
	  pWaveStruct->length = 18;
	    wrev(legd8, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd8, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 9:
//       pFilterCoef = legd9;
	  pWaveStruct->length = 20;
	    wrev(legd9, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd9, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    default:
      printf("legd%d is not available!\n",member);
      exit(0);
    }

//   wrev(pFilterCoef, pWaveStruct->length,
//        LowDecomFilCoef, pWaveStruct->length);
//   qmf_wrev(pFilterCoef, pWaveStruct->length,
// 	   HiDecomFilCoef, pWaveStruct->length);
  //sum = 0;
  /*for (count = 0; count < pWaveStruct->length; count++)
    LowDecomFilCoef[count] *= sqrt(2.0);

  for (count = 0; count < pWaveStruct->length; count++)
    HiDecomFilCoef[count] *= sqrt(2.0);*/
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
legendre_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;
  //double sum;
  //int count;


  switch (member)
    {
    case 1:
//       pFilterCoef = legd1;
	  pWaveStruct->length = 2;
	    verbatim_copy(legd1, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd1, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = legd2;
	  pWaveStruct->length = 4;
	    verbatim_copy(legd2, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd2, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 3:
//       pFilterCoef = legd3;
	  pWaveStruct->length = 6;
	    verbatim_copy(legd3, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd3, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 4:
//       pFilterCoef = legd4;
	  pWaveStruct->length = 8;
	    verbatim_copy(legd4, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd4, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 5:
//       pFilterCoef = legd5;
	  pWaveStruct->length = 10;
	    verbatim_copy(legd5, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd5, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 6:
//       pFilterCoef = legd6;
	  pWaveStruct->length = 14;
	    verbatim_copy(legd6, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd6, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 7:
//       pFilterCoef = legd7;
	  pWaveStruct->length = 16;
	    verbatim_copy(legd7, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd7, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 8:
//       pFilterCoef = legd8;
	  pWaveStruct->length = 18;
	    verbatim_copy(legd8, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd8, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 9:
//       pFilterCoef = legd9;
	  pWaveStruct->length = 20;
	    verbatim_copy(legd9, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd9, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    default:
      printf("legd%d is not available!\n",member);
      exit(0);
    }

//   verbatim_copy(pFilterCoef, pWaveStruct->length,
// 		LowReconFilCoef, pWaveStruct->length);
//   qmf_even(pFilterCoef, pWaveStruct->length,
//       HiReconFilCoef, pWaveStruct->length);
  /*for (count = 0; count < pWaveStruct->length; count++)
    LowDecomFilCoef[count] *= sqrt(2.0);

  for (count = 0; count < pWaveStruct->length; count++)
    HiDecomFilCoef[count] *= sqrt(2.0);*/
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
