/*
 * -------------------------------------------------------------------------
 * cowt.c -- dual tree complex wavelet transform
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

void
cowavedec (double *sigIn, int sigInLength, double *sigOutR,
	   double *sigOutI, int sigOutLength,
	   double *lowDTree1S1, double *hiDTree1S1,
	   double *lowDTree2S1, double *hiDTree2S1,
	   double *lowDTree1S2, double *hiDTree1S2,
	   double *lowDTree2S2, double *hiDTree2S2,
	   int filterLen, int *waveDecLengthArray,
	   int lengthArrayLength, int stride, extend_method extMethod)
{
  int ln1;
  double *sT1Low, *sT2Low;

  ln1 = waveDecLengthArray[lengthArrayLength-2];
  //printf("%d\n",filterLen);

  sT1Low = malloc(ln1*sizeof(double));
  sT2Low = malloc(ln1*sizeof(double));
  dwt_neo (sigIn, sigInLength, lowDTree1S1, hiDTree1S1, filterLen,
	   sT1Low, sigOutR+sigOutLength-ln1, ln1, extMethod);
  dwt_neo (sigIn, sigInLength, lowDTree2S1, hiDTree2S1, filterLen,
	   sT2Low, sigOutI+sigOutLength-ln1, ln1, extMethod);
  //printf("after dwt\n");

  wavedec (sT1Low, ln1, sigOutR, sigOutLength-ln1,
     lowDTree1S2, hiDTree1S2, filterLen, waveDecLengthArray,
     lengthArrayLength-1, stride-1, extMethod);
  free(sT1Low);
  wavedec (sT2Low, ln1, sigOutI, sigOutLength-ln1,
     lowDTree2S2, hiDTree2S2, filterLen, waveDecLengthArray,
     lengthArrayLength-1, stride-1, extMethod);
  free(sT2Low);

  return;
}

void
cowaverec (double *sigInR, double *sigInI, int sigInLength,
	   double *sigOut, int sigOutLength,
	   double *lowRTree1S1, double *hiRTree1S1,
	   double *lowRTree2S1, double *hiRTree2S1,
	   double *lowRTree1S2, double *hiRTree1S2,
	   double *lowRTree2S2, double *hiRTree2S2,
	   int filterLen, int *waveDecLengthArray,
	   int lengthArrayLength, int stride,
	   extend_method extMethod)
{
  int ln1, i;
  double *sT1Low, *sT2Low, *sigOut1, *sigOut2;

  ln1 = waveDecLengthArray[lengthArrayLength-2];

  sT1Low = malloc(ln1*sizeof(double));
  sT2Low = malloc(ln1*sizeof(double));

  waverec (sigInR, sigInLength-ln1, sT1Low, ln1,
	   lowRTree1S2, hiRTree1S2, filterLen,
	   waveDecLengthArray, lengthArrayLength-1,
	   stride-1, extMethod);
  waverec (sigInI, sigInLength-ln1, sT2Low, ln1,
	   lowRTree2S2, hiRTree2S2, filterLen,
	   waveDecLengthArray, lengthArrayLength-1,
	   stride-1, extMethod);

  sigOut1 = malloc(sigOutLength*sizeof(double));
  sigOut2 = malloc(sigOutLength*sizeof(double));

  idwt_neo (sT1Low, sigInR+sigInLength-ln1, ln1, lowRTree1S1,
	    hiRTree1S1, filterLen, sigOut1, sigOutLength);
  free(sT1Low);
  idwt_neo (sT2Low, sigInI+sigInLength-ln1, ln1, lowRTree2S1,
	    hiRTree2S1, filterLen, sigOut2, sigOutLength);
  free(sT2Low);

  for(i=0;i<sigOutLength;i++)
    sigOut[i] = (sigOut1[i] + sigOut2[i]) / 2;
  free(sigOut1);
  free(sigOut2);
  return;
}

void
copmd (double *matrixInR, double *matrixInI, int sigInLength,
      int InRow, int InCol, double *matrixOutR, double *matrixOutI)
{
  int i;
  verbatim_copy(matrixInR,InRow*InCol,matrixOutR,InRow*InCol);
  verbatim_copy(matrixInI,InRow*InCol,matrixOutI,InRow*InCol);
  for(i=InRow*InCol;i<sigInLength;i++)
    {
      matrixOutR[i] = (matrixInR[i]+matrixInI[i])/2;
      matrixOutI[i] = (matrixInR[i]-matrixInI[i])/2;
    }
  return;
}

void
copmr (double *matrixInR, double *matrixInI, int sigInLength,
      int InRow, int InCol, double *matrixOutR, double *matrixOutI)
{
  int i;
  verbatim_copy(matrixInR,InRow*InCol,matrixOutR,InRow*InCol);
  verbatim_copy(matrixInI,InRow*InCol,matrixOutI,InRow*InCol);
  for(i=InRow*InCol;i<sigInLength;i++)
    {
      matrixOutR[i] = matrixInR[i]+matrixInI[i];
      matrixOutI[i] = matrixInR[i]-matrixInI[i];
    }
  return;
}


void
cowavedec2 (double *matrixIn, int matrixInRow, int matrixInCol,
	    double *lowDTree1S1, double *hiDTree1S1,
	    double *lowDTree1S2, double *hiDTree1S2,
	    int filterLen, int *pLen, double *coef,
	    int sigOutLength, int stride, extend_method extMethod)
{
  int level, total;
  int *pLen1, *pLen2;
  double *ar;


  /* first stage */
  level = 1;
  pLen1 = malloc((level+2)*2*sizeof(int));
  matrix_wavedec_len_cal (matrixInRow, matrixInCol, level,
			  filterLen, pLen1);
  wave_mem_cal (pLen1, level, &total);
  wavedec2 (matrixIn, matrixInRow, matrixInCol,
	    lowDTree1S1, hiDTree1S1, filterLen,
	    pLen1, coef+sigOutLength-4*pLen1[0]*pLen1[1],
	    total, level, extMethod);
  ar = malloc(pLen1[0]*pLen1[1]*sizeof(double));
  verbatim_copy(coef+sigOutLength-4*pLen1[0]*pLen1[1],
		pLen1[0]*pLen1[1], ar, pLen1[0]*pLen1[1]);
  /* further stage */
  pLen2 = malloc((stride+1)*2*sizeof(int));
  matrix_wavedec_len_cal (pLen1[0], pLen1[1], stride-1,
			  filterLen, pLen2);
  wave_mem_cal (pLen2, stride-1, &total);
  wavedec2 (ar, pLen1[0], pLen1[1], lowDTree1S2, hiDTree1S2,
	    filterLen, pLen2, coef, total, stride-1, extMethod);

  free(ar);
  free(pLen1);
  free(pLen2);

  return;
}


void
cowaverec2 (double *coef, int sigInLength,
	    double *lowRTree1S1, double *hiRTree1S1,
	    double *lowRTree1S2, double *hiRTree1S2,
	    int filterLen, double *matrixOut, int matrixOutRow,
	    int matrixOutCol, int *pLen, int stride,
	    extend_method extMethod)
{
  int level, total;
  int *pLen1, *pLen2;
  double *ar;


  level = 1;
  pLen1 = malloc((level+2)*2*sizeof(int));
  matrix_wavedec_len_cal (matrixOutRow, matrixOutCol, level,
			  filterLen, pLen1);

  pLen2 = malloc((stride+1)*2*sizeof(int));
  matrix_wavedec_len_cal (pLen1[0], pLen1[1], stride-1,
			  filterLen, pLen2);
  wave_mem_cal (pLen2, stride-1, &total);
  ar = malloc(pLen1[0]*pLen1[1]*sizeof(double));
  waverec2 (coef, total, lowRTree1S2, hiRTree1S2,
	    filterLen, ar, pLen1[0], pLen1[1],
	    pLen2, stride-1, extMethod);

  verbatim_copy(ar, pLen1[0]*pLen1[1],
		coef+sigInLength-4*pLen1[0]*pLen1[1],
		pLen1[0]*pLen1[1]);

  waverec2 (coef+sigInLength-4*pLen1[0]*pLen1[1],
	    4*pLen1[0]*pLen1[1], lowRTree1S1, hiRTree1S1,
	    filterLen, matrixOut, matrixOutRow, matrixOutCol,
	    pLen1, level, extMethod);
  free(ar);
  free(pLen1);
  free(pLen2);

  return;
}

void
cowavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
	     double *lowDTree1S1R, double *hiDTree1S1R,
	     double *lowDTree1S1C, double *hiDTree1S1C,
	     double *lowDTree1S2R, double *hiDTree1S2R,
	     double *lowDTree1S2C, double *hiDTree1S2C,
	     int filterLen, int *pLen, double *coef,
	     int sigOutLength, int stride, extend_method extMethod)
{
  int level, total;
  int *pLen1, *pLen2;
  double *ar;


  /* first stage */
  level = 1;
  pLen1 = malloc((level+2)*2*sizeof(int));
  matrix_wavedec_len_cal (matrixInRow, matrixInCol, level,
			  filterLen, pLen1);
  wave_mem_cal (pLen1, level, &total);
  wavedec2a (matrixIn, matrixInRow, matrixInCol,
	     lowDTree1S1R, hiDTree1S1R,
	     lowDTree1S1C, hiDTree1S1C,filterLen,
	     pLen1, coef+sigOutLength-4*pLen1[0]*pLen1[1],
	     total, level, extMethod);
  ar = malloc(pLen1[0]*pLen1[1]*sizeof(double));
  verbatim_copy(coef+sigOutLength-4*pLen1[0]*pLen1[1],
		pLen1[0]*pLen1[1], ar, pLen1[0]*pLen1[1]);
  /* further stage */
  pLen2 = malloc((stride+1)*2*sizeof(int));
  matrix_wavedec_len_cal (pLen1[0], pLen1[1], stride-1,
			  filterLen, pLen2);
  wave_mem_cal (pLen2, stride-1, &total);
  wavedec2a (ar, pLen1[0], pLen1[1], lowDTree1S2R, hiDTree1S2R,
	     lowDTree1S2C, hiDTree1S2C,
	     filterLen, pLen2, coef, total, stride-1, extMethod);

  free(ar);
  free(pLen1);
  free(pLen2);

  return;
}

void
cowaverec2a (double *coef, int sigInLength,
	     double *lowRTree1S1R, double *hiRTree1S1R,
	     double *lowRTree1S1C, double *hiRTree1S1C,
	     double *lowRTree1S2R, double *hiRTree1S2R,
	     double *lowRTree1S2C, double *hiRTree1S2C,
	     int filterLen, double *matrixOut, int matrixOutRow,
	     int matrixOutCol, int *pLen, int stride,
	     extend_method extMethod)
{
  int level, total;
  int *pLen1, *pLen2;
  double *ar;


  level = 1;
  pLen1 = malloc((level+2)*2*sizeof(int));
  matrix_wavedec_len_cal (matrixOutRow, matrixOutCol, level,
			  filterLen, pLen1);

  pLen2 = malloc((stride+1)*2*sizeof(int));
  matrix_wavedec_len_cal (pLen1[0], pLen1[1], stride-1,
			  filterLen, pLen2);
  wave_mem_cal (pLen2, stride-1, &total);
  ar = malloc(pLen1[0]*pLen1[1]*sizeof(double));
  waverec2a (coef, total, lowRTree1S2R, hiRTree1S2R,
	    lowRTree1S2C, hiRTree1S2C,
	    filterLen, ar, pLen1[0], pLen1[1],
	    pLen2, stride-1, extMethod);

  verbatim_copy(ar, pLen1[0]*pLen1[1],
		coef+sigInLength-4*pLen1[0]*pLen1[1],
		pLen1[0]*pLen1[1]);

  waverec2a (coef+sigInLength-4*pLen1[0]*pLen1[1],
	    4*pLen1[0]*pLen1[1], lowRTree1S1R, hiRTree1S1R,
	    lowRTree1S1C, hiRTree1S1C,
	    filterLen, matrixOut, matrixOutRow, matrixOutCol,
	    pLen1, level, extMethod);
  free(ar);
  free(pLen1);
  free(pLen2);

  return;
}
