/*
 * -------------------------------------------------------------------------
 * kingsbury.c -- Kingsbury Q-filter coefficents.
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

static const double ksq1[10] = { 0,
				 0,
				 -0.11430184000000 ,
				 0,
				 0.58751830000000 ,
				 0.76027237000000 ,
				 0.23389032000000 ,
				 -0.08832942000000 ,
				 0,
				 0.03516384000000};

static const double ksq2[10] = { 0.03516384000000,
				 0,
				 -0.08832942000000,
				 0.23389032000000,
				 0.76027237000000,
				 0.58751830000000,
				 0,
				 -0.11430184000000,
				 0,
				 0};


/*********************************************
 * Global Function
 ********************************************/

void
kingsburyq_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
  int i;
//   double *pFilterCoef;

  pWaveStruct->length = 10;

  switch (member)
    {
    case 1:
//       pFilterCoef = ksq1;
        wrev(ksq1, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(ksq1, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = ksq2;
        wrev(ksq2, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(ksq2, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    default:
      printf("ksq%d is not available!\n",member);
      exit(0);
    }

//   wrev(pFilterCoef, pWaveStruct->length,
//        LowDecomFilCoef, pWaveStruct->length);
//   qmf_wrev(pFilterCoef, pWaveStruct->length,
// 	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  if (member==1)
    {
      for(i=0;i<10;i++)
	HiDecomFilCoef[i] *= -1;
    }
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
kingsburyq_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
  int i;
//   double *pFilterCoef;

  pWaveStruct->length = 10;

  switch (member)
    {
    case 1:
//       pFilterCoef = ksq1;
      verbatim_copy(ksq1, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(ksq1, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 2:
//       pFilterCoef = ksq2;
            verbatim_copy(ksq2, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(ksq2, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    default:
      printf("ksq%d is not available!\n",member);
      exit(0);
    }

//   verbatim_copy(pFilterCoef, pWaveStruct->length,
// 		LowReconFilCoef, pWaveStruct->length);
//   qmf_even(pFilterCoef, pWaveStruct->length,
//       HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  if (member==1)
    {
      for(i=0;i<10;i++)
	HiReconFilCoef[i] *= -1;
    }
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
