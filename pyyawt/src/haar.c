/*
 * -------------------------------------------------------------------------
 * haar.c -- Haar wavelets coefficients.
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
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

static const double haar[2] = {
     7.071067811865475244008443621048490392848359376884740365883398e-01,
 7.071067811865475244008443621048490392848359376884740365883398e-01
 //0.70710678, 0.70710678
};

/*********************************************
 * Global Function
 ********************************************/

void
haar_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = haar;

  pWaveStruct->length = 2;

  wrev(haar, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(haar, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}

void
haar_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = haar;
  pWaveStruct->length = 2;

  verbatim_copy(haar, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(haar, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
