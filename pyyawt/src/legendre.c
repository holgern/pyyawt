/*
 * -------------------------------------------------------------------------
 * legendre.c -- Legendre wavelets coefficents.
 * PYYAWT - yet another python wavelet toolbox
 * Copyright (C) 2005-2007  Roger Liu
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

#include "swtlib.h"

/*********************************************
 * Local Variable (Filter Coefficent)
 * 2*N = nu+1
 * from sympy import Symbol, Rational, binomial
 * for k in np.arange(nu+1): print("%.32f, "%(np.sqrt(2)/2**(2*nu)*binomial(2*k, k)*binomial(2*nu-2*k,nu-k)))
 ********************************************/

static const double legd1[2] = {
0.70710678118654757273731092936941, 
0.70710678118654757273731092936941
};

static const double legd2[4] = {
0.44194173824159221908303152304143, 
0.26516504294495535365427940632799, 
0.26516504294495535365427940632799, 
0.44194173824159221908303152304143
};

static const double legd3[6] = {
0.34802911886525389473234781689825, 
0.19334951048069659584882629133062, 
0.16572815184059708215613682114054, 
0.16572815184059708215613682114054, 
0.19334951048069659584882629133062, 
0.34802911886525389473234781689825
};

static const double legd4[8] = {
0.29623907141506727880297944466292, 
0.15951334614657469712639681347355, 
0.13051091957447022440241823915130, 
0.12084344405043537240551643208164, 
0.12084344405043537240551643208164, 
0.13051091957447019664684262352239, 
0.15951334614657469712639681347355, 
0.29623907141506727880297944466292
};

static const double legd5[10] = {
0.26229501114875747314059140080644, 
0.13886206472581280602085485043062, 
0.11108965178065023649001119565582, 
0.09969584134160917876510410451374, 
0.09516421218971785056517376233387, 
0.09516421218971786444296157014833, 
0.09969584134160919264289191232820, 
0.11108965178065022261222338784137, 
0.13886206472581280602085485043062, 
0.26229501114875747314059140080644};

static const double legd6[12] = {
0.23785388510989599608613787040667, 
0.12459013029565980945623238085318, 
0.09836062918078407324440348702410, 
0.08678879045363299682414037761191, 
0.08100287109005746555290272681304, 
0.07851047505651723157349408666050, 
0.07851047505651723157349408666050, 
0.08100287109005746555290272681304, 
0.08678879045363299682414037761191, 
0.09836062918078405936661567920964, 
0.12459013029565980945623238085318, 
0.23785388510989599608613787040667};

static const double legd7[14] = {
  0.21917625631120224438674881639599, 
0.11397165328182516652599787221334, 
0.08919520691621100894064255726335, 
0.07786883143478738611431566596366, 
0.07172129211098837719973175808263, 
0.06834617248223598395817646178330, 
0.06682736864929740561169779766715, 
0.06682736864929740561169779766715, 
0.06834617248223598395817646178330, 
0.07172129211098837719973175808263, 
0.07786883143478738611431566596366, 
0.08919520691621099506285474944889, 
0.11397165328182516652599787221334, 
0.21917625631120224438674881639599};

static const double legd8[16] = {
0.20430358177579924228162155941391, 
0.10567426643575822398091190734704, 
0.08219109611670084858392471005573, 
0.07123228330114073081347214611014, 
0.06503817170973719141713331737265, 
0.06132170475489506777933002013015, 
0.05917006599156541396533626198107, 
0.05817561110095086779336881477320, 
0.05817561110095087473226271868043, 
0.05917006599156540702644235807384, 
0.06132170475489506084043611622292, 
0.06503817170973719141713331737265, 
0.07123228330114073081347214611014, 
0.08219109611670083470613690224127, 
0.10567426643575822398091190734704, 
0.20430358177579924228162155941391
};

static const double legd9[18] = {
0.19209794499691418279141430502932, 
0.09895954742265275494439435988170, 
0.07661384316592471932505503673383, 
0.06604641652234888304917603818467, 
0.05993100758509436326582076048908, 
0.05609542309964832967894565740608, 
0.05365649166053318430691376761388, 
0.05219645107113091153649975240114, 
0.05150965566230025077798515553695, 
0.05150965566230024383909125162972, 
0.05219645107113091847539365630837, 
0.05365649166053317736801986370665, 
0.05609542309964832274005175349885, 
0.05993100758509436326582076048908, 
0.06604641652234888304917603818467, 
0.07661384316592471932505503673383, 
0.09895954742265275494439435988170, 
0.19209794499691418279141430502932
};

/*********************************************
 * Global Function
 ********************************************/

void
legendre_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;
  //double sum;
  int count;


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
	  pWaveStruct->length = 12;
	    wrev(legd6, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd6, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 7:
//       pFilterCoef = legd7;
	  pWaveStruct->length = 14;
	    wrev(legd7, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd7, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 8:
//       pFilterCoef = legd8;
	  pWaveStruct->length = 16;
	    wrev(legd8, pWaveStruct->length,
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(legd8, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
      break;
    case 9:
//       pFilterCoef = legd9;
	  pWaveStruct->length = 18;
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
  //for (count = 0; count < pWaveStruct->length; count++)
  //  LowDecomFilCoef[count] *= sqrt(2.0);

  //for (count = 0; count < pWaveStruct->length; count++)
  //  HiDecomFilCoef[count] *= sqrt(2.0);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
legendre_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
//   double *pFilterCoef;
  //double sum;
  int count;


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
	  pWaveStruct->length = 12;
	    verbatim_copy(legd6, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd6, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 7:
//       pFilterCoef = legd7;
	  pWaveStruct->length = 14;
	    verbatim_copy(legd7, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd7, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 8:
//       pFilterCoef = legd8;
	  pWaveStruct->length = 16;
	    verbatim_copy(legd8, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(legd8, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length);
      break;
    case 9:
//       pFilterCoef = legd9;
	  pWaveStruct->length = 18;
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
  //for (count = 0; count < pWaveStruct->length; count++)
  //  LowReconFilCoef[count] *= sqrt(2.0);

  //for (count = 0; count < pWaveStruct->length; count++)
  //  HiReconFilCoef[count] *= sqrt(2.0);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
