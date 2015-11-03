/*
 * -------------------------------------------------------------------------
 * farras.c -- Farras filter coefficents.
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2008  Roger Liu
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

static const double fa1[10] = { 0,
				0.01122679215254 ,
				0.01122679215254 ,
				-0.08838834764832 ,
				0.08838834764832 ,
				0.69587998903400 ,
				0.69587998903400 ,
				0.08838834764832 ,
				-0.08838834764832,
				0};

static const double fa2[10] = { 0,
				0,
				-0.08838834764832,
				0.08838834764832,
				0.69587998903400,
				0.69587998903400,
				0.08838834764832,
				-0.08838834764832,
				0.01122679215254,
				0.01122679215254};

/*********************************************
 * Global Function
 ********************************************/

void
farras_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
  int i;
//   double *pFilterCoef;

  pWaveStruct->length = 10;

  switch (member)
    {
    case 1:
//       pFilterCoef = fa1;
        wrev(fa1, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(fa1, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = fa2;
        wrev(fa2, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(fa2, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    default:
      printf("fa%d is not available!\n",member);
      exit(0);
    }

//   wrev(pFilterCoef, pWaveStruct->length,
//        LowDecomFilCoef, pWaveStruct->length);
//   qmf_wrev(pFilterCoef, pWaveStruct->length,
// 	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  for(i=0;i<10;i++)
    HiDecomFilCoef[i] *= -1;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
farras_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
  int i;
//   double *pFilterCoef;

  pWaveStruct->length = 10;

  switch (member)
    {
    case 1:
//       pFilterCoef = fa1;
        verbatim_copy(fa1, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(fa1, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = fa2;
        verbatim_copy(fa2, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(fa2, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    default:
      printf("fa%d is not available!\n",member);
      exit(0);
    }

//   verbatim_copy(pFilterCoef, pWaveStruct->length,
// 		LowReconFilCoef, pWaveStruct->length);
//   qmf_even(pFilterCoef, pWaveStruct->length,
//       HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  for(i=0;i<10;i++)
    HiReconFilCoef[i] *= -1;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
