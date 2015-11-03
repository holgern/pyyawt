/*
 * -------------------------------------------------------------------------
 * beylkin.c -- Beylkin wavelets coefficents.
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

static const double beylkin[18] = {
	0.099305765374,	0.424215360813,	0.699825214057,
	0.449718251149, -0.110927598348,-0.264497231446,
	0.026900308804, 0.155538731877, -0.017520746267,
	-0.088543630623, 0.019679866044, 0.042916387274,
	-0.017460408696, -0.014365807969, 0.010040411845,
	0.001484234782, -0.002736031626, 0.000640485329};

/*********************************************
 * Global Function
 ********************************************/

void
beylkin_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = beylkin;

  pWaveStruct->length = 18;

  wrev(beylkin, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(beylkin, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}

void
beylkin_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = beylkin;
  pWaveStruct->length = 18;

  verbatim_copy(beylkin, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(beylkin, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
