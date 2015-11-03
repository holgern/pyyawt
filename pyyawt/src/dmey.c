/*
 * -------------------------------------------------------------------------
 * dmey.c -- Beylkin wavelets coefficents.
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


static const double dmey[62] = {
-1.0099999569414229e-012, 8.519459636796214e-009,
-1.111944952595278e-008, -1.0798819539621958e-008,
6.0669757413511352e-008, -1.0866516536735883e-007,
8.2006806503864813e-008, 1.1783004497663934e-007,
-5.5063405652522782e-007, 1.1307947017916706e-006,
-1.4895492164971559e-006, 7.367572885903746e-007,
3.2054419133447798e-006, -1.6312699734552807e-005,
6.5543059305751491e-005, -0.00060115023435160925,
-0.002704672124643725, 0.0022025341009110021,
0.006045814097323304, -0.0063877183184971563,
-0.011061496392513451, 0.015270015130934803,
0.017423434103729693, -0.032130793990211758,
-0.024348745906078023, 0.063739024322801596,
0.030655091960824263, -0.13284520043622938,
-0.035087555656258346, 0.44459300275757724,
0.74458559231880628, 0.44459300275757724,
-0.035087555656258346, -0.13284520043622938,
0.030655091960824263, 0.063739024322801596,
-0.024348745906078023, -0.032130793990211758,
0.017423434103729693, 0.015270015130934803,
-0.011061496392513451, -0.0063877183184971563,
0.006045814097323304, 0.0022025341009110021,
-0.002704672124643725, -0.00060115023435160925,
6.5543059305751491e-005, -1.6312699734552807e-005,
3.2054419133447798e-006, 7.367572885903746e-007,
-1.4895492164971559e-006, 1.1307947017916706e-006,
-5.5063405652522782e-007, 1.1783004497663934e-007,
8.2006806503864813e-008, -1.0866516536735883e-007,
6.0669757413511352e-008, -1.0798819539621958e-008,
-1.111944952595278e-008, 8.519459636796214e-009,
-1.0099999569414229e-012, 0.0};

/*********************************************
 * Global Function
 ********************************************/

void
dmey_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = dmey;

  pWaveStruct->length = 62;

  wrev(dmey, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(dmey, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}

void
dmey_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;

//   pFilterCoef = dmey;
  pWaveStruct->length = 62;

  verbatim_copy(dmey, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(dmey, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
