/*
 * -------------------------------------------------------------------------
 * vaidyanathan.c -- Vaidyanathan wavelets coefficents.
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

static const double vaidyanathan[24] = {
   -0.000062906118,	0.000343631905,	-0.000453956620,
   -0.000944897136,	0.002843834547,	0.000708137504,
   -0.008839103409,	0.003153847056,	0.019687215010,
   -0.014853448005,	-0.035470398607, 0.038742619293,
	0.055892523691,	-0.077709750902, -0.083928884366,
	0.131971661417,	0.135084227129,	-0.194450471766,
	-0.263494802488, 0.201612161775, 0.635601059872,
	0.572797793211,	0.250184129505,	0.045799334111};

/*********************************************
 * Global Function
 ********************************************/

void
vaidyanathan_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = vaidyanathan;

  pWaveStruct->length = 24;

  wrev(vaidyanathan, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(vaidyanathan, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}

void
vaidyanathan_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = vaidyanathan;
  pWaveStruct->length = 24;

  verbatim_copy(vaidyanathan, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(vaidyanathan, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
