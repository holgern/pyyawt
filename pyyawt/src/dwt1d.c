/*
 * -------------------------------------------------------------------------
 * dwt1d.c -- 1-D signal decomposition and reconstruction
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
/*-------------------------------------------
 * wavelet name parser
 *-----------------------------------------*/

 extend_method getdwtMode(){
   return dwtMode;
 }

void
filter_clear ()
{
  int count;
  for(count=0;count<30;count++)
    {
      LowDecomFilCoef[count] = 0;
      LowReconFilCoef[count] = 0;
      HiDecomFilCoef[count] = 0;
      HiReconFilCoef[count] = 0;
    }
  return;
}

/*-------------------------------------------
 * Orthfilter Group
 *-----------------------------------------*/

void
orth_filt_group (double *filterIn, int sigInLength,
		 double *filterLowRec, double *filterLowDec,
		 double *filterHiRec, double *filterHiDec)
{
  int count;

  for (count = 0; count < sigInLength; count++)
    filterLowRec[count] = filterIn[count];
  wrev (filterLowRec, sigInLength, filterLowDec, sigInLength);
  qmf_even (filterLowRec, sigInLength, filterHiRec, sigInLength);
  wrev (filterHiRec, sigInLength, filterHiDec, sigInLength);
  return;
}

void
bior_filt_group (double *f1, int sigInLength1,
		 double *f2, int sigInLength2,
		 double *lowDecom, int sigOutLength1,
		 double *hiDecom, int sigOutLength2,
		 double *lowRecon, int sigOutLength3,
		 double *hiRecon, int sigOutLength4)
{
  verbatim_copy (f2, sigInLength2, lowRecon, sigOutLength3);
  wrev (f1, sigInLength1, lowDecom, sigOutLength1);
  qmf_even (f1, sigInLength1, hiRecon, sigOutLength4);
  //qmf_odd (f1, sigInLength1, hiRecon, sigOutLength4);
  qmf_wrev (f2, sigInLength2, hiDecom, sigOutLength2);
  return;
}

/*-------------------------------------------
 * wavelet name parser
 *-----------------------------------------*/

void
wavelet_parser (char *wname, int *family, int *member)
{
  int count;

  *family = NOT_DEFINED;
  *member = NOT_DEFINED;
  for(count=0;count<waveletIdentityNum;count++)
    {
      if (strcmp(wname,wi[count].wname) == 0)
	{
	  *family = wi[count].family;
	  *member = wi[count].member;
	  break;
	}
    }
  return;
}

void
wavelet_fun_parser (char *wname, int *ii)
{
  int count;

  *ii = -1;
  for(count=0;count<waveletIdentityNum;count++)
    {
      if (strcmp(wname,wi[count].wname) == 0)
	{
	  *ii = count;
	  break;
	}
    }
  return;
}

/*void
wave_len_validate (int sigInLen, int waveLength, int *lev, int *val)
{
  int n;
  int m;
  int n1;
  int m1;
  float di;

  *val = 0;
  di = (float) sigInLen / (float) waveLength;
  if (di < 1)
    {
      *lev = 0;
      *val = 0;
      return;
    }
  else
    {
      n = (int) floor (log (di) / log ((float) 2));
      m = (int) ceil (log (di) / log ((float) 2));
      if ((((long) 1 << n) * waveLength == sigInLen)
	  || (((long) 1 << m) * waveLength == sigInLen))
	*lev = m + 1;
      else
	*lev = n + 1;
      *val = 1;
      //	  n1 = (int) floor (log (waveLength) / log ((float) 2));
      //m1 = (int) ceil (log (waveLength) / log ((float) 2));
      //if (n1 != m1)
      //    	  *lev = (*lev) - 1;
      if (di == 2)
	*lev -= 1;
      return;
    }
    }*/

void
wave_len_validate (int sigInLen, int waveLength, int *lev, int *val)
{
  int n;
  int m;

  *val = 0;
  if (sigInLen < 2*waveLength)
    {
      *lev = 0;
      *val = 0;
      return;
    }
  else
    {
      *val = 1;
      *lev = 0;
      n = sigInLen;
      do
	{
	  m = (int)floor((n + waveLength - 1)/2);
	  *lev += 1;
	  n = m;
	}
      while (m>=2*waveLength);
      return;
    }
}



void
dwt_write (char *mode, int *errCode)
{
  int count;
  *errCode = UNKNOWN_INPUT_ERR;

  for (count=0;count<extensionIdentityNum;count++)
    {
      if (strcmp(mode,ei[count].extMethodName) == 0)
	{
	  dwtMode = ei[count].extMethod;
	  *errCode = SUCCESS;
	  break;
	}
    }
  return;
}

void
dwt_parse(char **str1)
{
  int count;
  for (count=0;count<extensionIdentityNum;count++)
    {
      if (ei[count].extMethod == dwtMode)
	{
	  strcpy(*str1,ei[count].extMethodName);
	  break;
	}
    }
  return;
}

void
dwt_nex (double *sigIn, int sigInLength, double *lowDe,
	 double *hiDe, int filterLen, double *approx,
	 double *detail, int sigOutLength)
{
  //int sigInLengthTemp,
  int sigOutLengthTemp, sigOutLengthPre;
  double *approxTemp, *approxPre;
  double *detailTemp, *detailPre;

  //sigInLengthTemp = sigInLength + 2 * (filterLen - 1);
  //sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  //wextend_1D_center (sigIn, sigInLength, sigInTemp,
  //     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLength + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, approxTemp,
	sigOutLengthTemp, lowDe, filterLen);
  sigOutLengthPre = sigOutLengthTemp/2;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 approxPre, sigOutLengthPre);
  wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);

  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, detailTemp,
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 detailPre, sigOutLengthPre);
  wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);

  free(detailPre);
  free(detailTemp);
  //free(sigInTemp);
  return;
}


void
dwt (double *sigIn, int sigInLength, double *lowDe,
     double *hiDe, int filterLen, double *approx,
     double *detail, int sigOutLength, extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, sigOutLengthPre;
  double *sigInTemp, *approxTemp, *approxPre;
  double *detailTemp, *detailPre;

  sigInLengthTemp = sigInLength + 2 * (filterLen - 1);
  sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  wextend_1D_center (sigIn, sigInLength, sigInTemp,
		     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, approxTemp,
	sigOutLengthTemp, lowDe, filterLen);
  sigOutLengthPre = sigOutLengthTemp/2;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 approxPre, sigOutLengthPre);
  wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);

  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, detailTemp,
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 detailPre, sigOutLengthPre);
  wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);

  free(detailPre);
  free(detailTemp);
  free(sigInTemp);
  return;
}

void
dwt_neo (double *sigIn, int sigInLength, double *lowDe,
     double *hiDe, int filterLen, double *approx,
     double *detail, int sigOutLength, extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, sigOutLengthPre;
  double *sigInTemp, *approxTemp, *approxPre;
  double *detailTemp, *detailPre;

  sigInLengthTemp = sigInLength + 2 * filterLen;
  if ((extMethod == PER)&&(sigInLength%2 != 0))
	  sigInLengthTemp = sigInLength + 1 + 2 * filterLen;
  //if ((extMethod == PER)&&(sigInLength%2 == 0))
	//  sigInLengthTemp = sigInLength + 2 * filterLen;
  sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  wextend_1D_center (sigIn, sigInLength, sigInTemp,
		     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, approxTemp,
	sigOutLengthTemp, lowDe, filterLen);
  //sigOutLengthPre = sigOutLengthTemp/2;
  sigOutLengthPre = sigInLength + filterLen - 1;
  if ((extMethod==PER)&&(sigInLength%2 == 0))
	  sigOutLengthPre = sigInLength;
  if ((extMethod==PER)&&(sigInLength%2 !=0))
	  sigOutLengthPre = sigInLength + 1;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  wkeep_1D_center (approxTemp, sigOutLengthTemp, approxPre, sigOutLengthPre);
  dyaddown_1D_keep_even (approxPre, sigOutLengthPre,approx,sigOutLength);
  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 //approxPre, sigOutLengthPre);
  //wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);

  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, detailTemp,
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  wkeep_1D_center (detailTemp, sigOutLengthTemp, detailPre, sigOutLengthPre);
  dyaddown_1D_keep_even (detailPre, sigOutLengthPre,detail,sigOutLength);
  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 //detailPre, sigOutLengthPre);
  //wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);

  free(detailPre);
  free(detailTemp);
  free(sigInTemp);
  return;
}


void
dwt_conv (double *sigIn, int sigInLength, double *lowDe,
     double *hiDe, int filterLen, double *approx,
     double *detail, int sigOutLength)
{

  //int sigOutLengthTemp;
  //double *approxTemp;
  //double *detailTemp;


  //sigOutLengthTemp = sigInLength + filterLen - 1;
  //approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  //conv (sigIn, sigInLength, approxTemp,
	//sigOutLengthTemp, lowDe, filterLen);
  conv (sigIn, sigInLength, approx,
	sigOutLength, lowDe, filterLen);

  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,approx,sigOutLength);

  //free(approxTemp);
  //detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  //conv (sigIn, sigInLength, detailTemp,
	//sigOutLengthTemp, hiDe, filterLen);
  conv (sigIn, sigInLength, detail,
	sigOutLength, hiDe, filterLen);

  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,detail,sigOutLength);

  //free(detailTemp);

  return;
}

void
dwt_no_extension (double *sigIn, int sigInLength, double *lowDe,
     double *hiDe, int filterLen, double *approx,
     double *detail, int sigOutLength)
{
  //int sigInLengthTemp,
  int sigOutLengthTemp;
  //sigOutLengthPre;
  //double *sigInTemp,
  double *approxTemp;//, *approxPre;
  double *detailTemp;//, *detailPre;

  //sigInLengthTemp = sigInLength + 2 * filterLen;
  //if ((extMethod == PER)||(sigInLength%2 != 0))
	//  sigInLengthTemp = sigInLength + 1 + 2 * filterLen;
  //sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  //wextend_1D_center (sigIn, sigInLength, sigInTemp,
	//	     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLength + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, approxTemp,
	sigOutLengthTemp, lowDe, filterLen);
  //sigOutLengthPre = sigOutLengthTemp/2;
  //sigOutLengthPre = sigInLength + filterLen - 1;
  //approxPre = malloc (sigOutLengthPre * sizeof(double));
  //wkeep_1D_center (approxTemp, sigOutLengthTemp, approxPre, sigOutLengthPre);
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,approx,sigOutLength);
  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 //approxPre, sigOutLengthPre);
  //wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);

  //free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, detailTemp,
	sigOutLengthTemp, hiDe, filterLen);
  //detailPre = malloc (sigOutLengthPre * sizeof(double));
  //wkeep_1D_center (detailTemp, sigOutLengthTemp, detailPre, sigOutLengthPre);
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,detail,sigOutLength);
  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 //detailPre, sigOutLengthPre);
  //wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);

  //free(detailPre);
  free(detailTemp);
  //free(sigInTemp);
  return;
}


void
idwt_complete_ex (double *approx, double *detail, int sigInLength,
		  double *lowRe, double *hiRe, int filterLen,
		  double *sigOut, int sigOutLength,
		  extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, count, ind, sigInLen;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre, *approxEx, *detailEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  approxEx = malloc(sigInLen*sizeof(double));
  detailEx = malloc(sigInLen*sizeof(double));

  wextend_1D_center (approx, sigInLength, approxEx, sigInLen,
		     extMethod);
  wextend_1D_center (detail, sigInLength, detailEx, sigInLen,
		     extMethod);

  sigInLengthTemp = 2 * sigInLen - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approxEx, sigInLen,
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_odd (detailEx, sigInLen,
		      detailTemp, sigInLengthTemp);
  free(approxEx);
  free(detailEx);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  //wkeep_1D_center (sigOutPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  wkeep_1D_index (sigOutPre, sigOutLengthTemp,
  	  sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}


void
idwt_complete (double *approx, double *detail, int sigInLength,
	       double *lowRe, double *hiRe, int filterLen,
	       double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp, count, ind;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre;

  //printf("enter idwt!\n");
  sigInLengthTemp = 2 * sigInLength - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approx, sigInLength,
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_odd (detail, sigInLength,
		      detailTemp, sigInLengthTemp);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  //wkeep_1D_center (sigOutPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (sigOutPre, sigOutLengthTemp,
  	  sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}



void
idwt_neo (double *approx, double *detail, int sigInLength,
	       double *lowRe, double *hiRe, int filterLen,
	       double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp, count;//, ind;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre;

  //printf("enter idwt!\n");
  sigInLengthTemp = 2 * sigInLength + 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (approx, sigInLength,
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_even (detail, sigInLength,
		      detailTemp, sigInLengthTemp);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  wkeep_1D_center (sigOutPre, sigOutLengthTemp,
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (sigOutPre, sigOutLengthTemp,
  	  //sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}


void
idwt_approx_ex (double *approx, int sigInLength,
		double *lowRe, int filterLen,
		double *sigOut, int sigOutLength,
		extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, ind, sigInLen;
  double *approxTemp, *approxPre, *approxEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  approxEx = malloc(sigInLen * sizeof(double));
  wextend_1D_center (approx, sigInLength, approxEx, sigInLen,
		     extMethod);

  sigInLengthTemp = 2 * sigInLen - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approxEx, sigInLen,
		      approxTemp, sigInLengthTemp);
  free(approxEx);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  //wkeep_1D_center (approxPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (approxPre, sigOutLengthTemp,
  	  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}



void
idwt_approx (double *approx, int sigInLength,
	     double *lowRe, int filterLen,
	     double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp, ind;
  double *approxTemp, *approxPre;

  sigInLengthTemp = 2 * sigInLength - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approx, sigInLength,
		       approxTemp, sigInLengthTemp);

  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  //wkeep_1D_center (approxPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (approxPre, sigOutLengthTemp,
  	  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}

void
idwt_approx_neo (double *approx, int sigInLength,
	     double *lowRe, int filterLen,
	     double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp;//, ind;
  double *approxTemp, *approxPre;

  sigInLengthTemp = 2 * sigInLength + 1;
  approxTemp = (double*)malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (approx, sigInLength,
		       approxTemp, sigInLengthTemp);

  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = (double*)malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre,
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  wkeep_1D_center (approxPre, sigOutLengthTemp,
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (approxPre, sigOutLengthTemp,
  	//  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}

void
idwt_detail (double *detail, int sigInLength,
	     double *hiRe, int filterLen,
	     double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp, ind;
  double *detailTemp, *detailPre;

  sigInLengthTemp = 2 * sigInLength - 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (detail, sigInLength,
		       detailTemp, sigInLengthTemp);

  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  //wkeep_1D_center (detailPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (detailPre, sigOutLengthTemp,
		  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}

void
idwt_detail_ex (double *detail, int sigInLength,
		double *hiRe, int filterLen,
		double *sigOut, int sigOutLength,
		extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, ind, sigInLen;
  double *detailTemp, *detailPre, *detailEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  detailEx = malloc(sigInLen * sizeof(double));
  wextend_1D_center (detail, sigInLength, detailEx, sigInLen,
		     extMethod);

  sigInLengthTemp = 2 * sigInLen - 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (detailEx, sigInLen,
		      detailTemp, sigInLengthTemp);
  free(detailEx);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  //wkeep_1D_center (detailPre, sigOutLengthTemp,
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (detailPre, sigOutLengthTemp,
		  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}


void
idwt_detail_neo (double *detail, int sigInLength,
	     double *hiRe, int filterLen,
	     double *sigOut, int sigOutLength)
{
  int sigInLengthTemp, sigOutLengthTemp;//, ind;
  double *detailTemp, *detailPre;

  sigInLengthTemp = 2 * sigInLength + 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (detail, sigInLength,
		       detailTemp, sigInLengthTemp);

  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre,
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  wkeep_1D_center (detailPre, sigOutLengthTemp,
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (detailPre, sigOutLengthTemp,
	//	  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}


void
wave_dec_len_cal (int filterLen, int sigLength,
		  int stride, int *waveDecLengthArray)
{
  int count = 0;
  int calLen;
  waveDecLengthArray[stride + 1] = sigLength;
  if (dwtMode!=PER)
    {
      calLen = sigLength;
      for (count = 0; count < stride; count++)
	{
	  calLen += (filterLen - 1);
	  waveDecLengthArray[stride-count]=(int)(floor(calLen/2));
	  calLen = *(waveDecLengthArray + stride - count);
	}
      waveDecLengthArray[0] = waveDecLengthArray[1];
    }
  else
    {
      for (count=stride; count > 0; count--)
	waveDecLengthArray[count] =
	  (int)ceil(((double)(waveDecLengthArray[count+1]))/2.0);
     waveDecLengthArray[0] = waveDecLengthArray[1];
    }
  return;
}

void
upcoef_len_cal (int sigInLength, int filterLen, int stride,
		int *sigOutLength, int *sigOutLengthDefault)
{
  int count;
  *sigOutLength = sigInLength;
  *sigOutLengthDefault = sigInLength;
//   if ((2*(*sigOutLength) - filterLen + 2)<0){ //this was implemented for cwt
      for(count=0;count<stride;count++)
      {
	// original version
	*sigOutLengthDefault = 2*(*sigOutLengthDefault) + filterLen - 1;
	*sigOutLength = 2*(*sigOutLength) + filterLen - 2;

      }
//     } else { //works with dwt
//       for(count=0;count<stride;count++)
// 	{
//
//     //version 1.14 - but does not work with cwt
// 	  *sigOutLengthDefault = 2*(*sigOutLengthDefault) + filterLen - 1;
// 	  *sigOutLength = 2*(*sigOutLength) - filterLen + 2;
//
// 	}
//      }
  return;
}

void
upcoef (double *sigIn, int sigInLength, double *lowRe,double *hiRe,
	int filterLen, double *sigOut, int sigOutLength,
	int defaultLength, char *coefType, int step)
{
  int count, sigInLengthTemp, leng;
  double *sigInTemp, *sigOutTemp;
  //version 1.14 fow dwt - but does not work with cwt
//    sigInLengthTemp = 2 * sigInLength - filterLen + 2;

//   // works with wavefun, cwt
      sigInLengthTemp = 2 * sigInLength + filterLen - 2;



  //sigInLengthTemp = 2 * sigInLength + filterLen - 1;
  sigInTemp = (double *) malloc(defaultLength*sizeof(double));

  if (strcmp(coefType,"a")==0)
  {
	  //sciprint("recognized\n");
// 	  printf("sigInLength %d, filterLen%d, sigInLengthTemp %d\n",sigInLength,filterLen,sigInLengthTemp);
	  idwt_approx_neo (sigIn, sigInLength, lowRe, filterLen,
		 sigInTemp, sigInLengthTemp);
// 	  sciprint("recognized\n");
  }
  else
    idwt_detail_neo (sigIn, sigInLength, hiRe, filterLen,
		 sigInTemp, sigInLengthTemp);

  if (step > 1)
    {
      sigOutTemp = (double *) malloc(defaultLength*sizeof(double));
      for(count=0;count<defaultLength;count++)
	sigOutTemp[count] = 0;
      leng = sigInLengthTemp;
      for(count=0;count<(step-1);count++)
	{
	  //printf("leng %d, filterLen%d, leng*2-filterLen+2 %d\n",leng,filterLen,leng*2-filterLen+2);

	  // original version
	  idwt_approx_neo (sigInTemp, leng, lowRe, filterLen,
	               sigOutTemp, leng*2+filterLen-2);
	  leng = leng*2+filterLen-2;

// 	  //version 1.14 - but does not work with cwt
// 	  idwt_approx_neo (sigInTemp, leng, lowRe, filterLen,
// 		       sigOutTemp, leng*2-filterLen+2);
// 	  //sciprint("ok\n");
// 	  leng = leng*2-filterLen+2;

	  verbatim_copy (sigOutTemp, leng, sigInTemp, leng);
	}
      sigInLengthTemp = leng;
      free(sigOutTemp);
    }


  wkeep_1D_center (sigInTemp, sigInLengthTemp, sigOut, sigOutLength);
  free(sigInTemp);
  return;
}

void
wavedec (double *sigIn, int sigInLength, double *sigOut,
	 int sigOutLength, double *lowDe, double *hiDe,
	 int filterLen, int *waveDecLengthArray,
	 int lengthArrayLengh, int stride, extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  //int filterLen;
  double *pApprox, *pDetail, *pApproxTemp;
  double *pSig;

  pSig = sigIn;
  sigInLen = sigInLength;

  pApprox = malloc (sigInLength * sizeof (double));
  pApproxTemp = malloc (sigInLength * sizeof (double));
  for (count = 0; count < sigInLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = sigOutLength - waveDecLengthArray[stride];
  pDetail = sigOut + pos;
  for (count = 0; count < stride; count++)
    {
      dwt_neo (pSig, sigInLen, lowDe, hiDe, filterLen, pApprox,
	   pDetail, waveDecLengthArray[stride - count], extMethod);
      for (countt = 0; countt < waveDecLengthArray[stride - count]; countt++)
	pApproxTemp[countt] = pApprox[countt];
      pSig = pApproxTemp;
      sigInLen = waveDecLengthArray[stride - count];
      pos = pos - waveDecLengthArray[stride - count - 1];
      pDetail = (sigOut + pos);
    }
  for (count = 0; count < sigInLen; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}


void
waverec (double *sigIn, int sigInLength, double *sigOut,
	 int sigOutLength, double *lowRe, double *hiRe,
	 int filterLen, int *waveDecLengthArray,
	 int lengthArraylength, int stride,
	 extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  double *pApprox, *pDetail, *pApproxTemp;

  sigInLen = waveDecLengthArray[1];

  pApprox = malloc (sigOutLength * sizeof (double));
  pApproxTemp = malloc (sigOutLength * sizeof (double));
  for (count = 0; count < sigOutLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = waveDecLengthArray[0];
  pDetail = sigIn + pos;
  for (count = 0; count < waveDecLengthArray[1]; count++)
    pApprox[count] = sigIn[count];

  for (count = 0; count < stride; count++)
    {
      idwt_neo (pApprox, pDetail, sigInLen, lowRe, hiRe,
		     filterLen, pApproxTemp,
		     waveDecLengthArray[count + 2]);
      for (countt = 0; countt < waveDecLengthArray[count + 2]; countt++)
	pApprox[countt] = pApproxTemp[countt];
      sigInLen = waveDecLengthArray[count + 2];
      pos += waveDecLengthArray[count + 1];
      pDetail = sigIn + pos;
    }
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}

void
wenergy (double *coef, int coefLen, int *lenArray, int arrayLen,
	 double *aE, int aELen, double *dE, int dELen)
{
  int count, countt, *pos;
  double ene;

  ene = 0;
  for(count=0;count<coefLen;count++)
    ene += coef[count]*coef[count];
  *aE = 0;
  for(count=0;count<lenArray[0];count++)
    *aE += coef[count]*coef[count];
  *aE = (*aE)*100/ene;

  pos = malloc (dELen*sizeof(int));
  for(count=0;count<dELen;count++)
    pos[count] = 0;
  pos[0] = lenArray[0];
  for(count=1;count<dELen;count++)
    pos[count] += (lenArray[count] + pos[count-1]);
  for(count=0;count<dELen;count++)
    {
      dE[count] = 0;
      for(countt=0;countt<lenArray[count+1];countt++)
	dE[count] += coef[pos[count]+countt]*coef[pos[count]+countt];
      dE[count] = dE[count]*100/ene;
    }
  free(pos);
  return;
}

void
detcoef (double *sigIn, int sigInLength, int *waveDecLengthArray,
	 int arrayLen, double *sigOut, int sigOutLength,
	 int stride, int level)
{
   int leng=0;
  int startCount=0, count=0;
//  printf ("level %d, leng %d count %d stride %d waveDecLengthArray %d\n",level,leng,count,stride,waveDecLengthArray[stride - level+1]);
  if (level != 0)
    {

      for (count = 0; count < level; count++){
           //sciprint("tmp %d",tmp);
			     leng +=waveDecLengthArray[stride - count];;
			printf("");

		  }
  }
    //printf ("sigInLength %d, leng %d\n",sigInLength,leng);
  startCount = sigInLength - leng;
  //printf ("level %d, leng %d startCount %d sigOutLength %d,sigInLength %d \n",level,leng,startCount,sigOutLength,sigInLength);
  for (count = startCount; count <= (startCount + sigOutLength - 1); count++){
    sigOut[count - startCount] = sigIn[count];
  }

  return;
}

void
appcoef (double *sigIn, int sigInLength, double *sigOut,
	 int sigOutLength, double *lowRe, double *hiRe,
	 int filterLen, int *waveDecLengthArray,
	 int lengthArraylength, int stride, int level,
	 extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  double *pApprox, *pDetail, *pApproxTemp;

  if (level == stride)
    {
      for (count = 0; count < waveDecLengthArray[stride - level + 1]; count++)
	sigOut[count] = sigIn[count];
      return;
    }

  sigInLen = waveDecLengthArray[1];

  pApprox = malloc (sigOutLength * sizeof (double));
  pApproxTemp = malloc (sigOutLength * sizeof (double));
  for (count = 0; count < sigOutLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = waveDecLengthArray[0];
  pDetail = sigIn + pos;
  for (count = 0; count < waveDecLengthArray[1]; count++)
    pApprox[count] = sigIn[count];

  for (count = 0; count < (stride - level); count++)
    {
      idwt_neo (pApprox, pDetail, sigInLen, lowRe, hiRe,
		     filterLen, pApproxTemp,
		     waveDecLengthArray[count + 2]);
      for (countt = 0; countt < waveDecLengthArray[count + 2]; countt++)
	pApprox[countt] = pApproxTemp[countt];
      sigInLen = waveDecLengthArray[count + 2];
      pos += waveDecLengthArray[count + 1];
      pDetail = sigIn + pos;
    }
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}

void
wrcoef (double *sigIn, int sigInLength, double *lowRe, double *hiRe,
	int filterLen, int *waveDecLengthArray, int arrayLen,
	double *sigOut, int sigOutLength, char *coefType,
	int stride, int level, extend_method extMethod)
{

  int count = 0;
  int startCount, endCount, leng;
  double *sigOutTemp;

  sigOutTemp = malloc (sigInLength * sizeof (double));

  if (level != 0)
    {
      leng = 0;
      for (count = 0; count < level; count++)
	leng += waveDecLengthArray[stride - count];
    }

  if (strcmp (coefType, "d") == 0)
    {
      for (count = 0; count < sigInLength; count++)
	sigOutTemp[count] = 0;
      if (level != 0)
	{
	  startCount = sigInLength - leng;
	  endCount = startCount + waveDecLengthArray[stride - level + 1] - 1;
	  for (count = startCount; count <= endCount; count++)
	    sigOutTemp[count] = sigIn[count];
	}
    }
  else
    {
      for (count = 0; count < sigInLength; count++)
	sigOutTemp[count] = sigIn[count];
      if (level != 0)
	{
	  endCount = sigInLength - 1;
	  startCount = endCount - leng + 1;
	  for (count = startCount; count <= endCount; count++)
	    sigOutTemp[count] = 0;
	}
    }
  waverec (sigOutTemp, sigInLength, sigOut, sigOutLength,
	   lowRe, hiRe, filterLen, waveDecLengthArray, arrayLen,
	   stride, extMethod);
  //waverec (sigInLength, sigOutTemp, sigOutLength, sigOut,
  //   stride, waveDecLengthArray, waveType, mode);
  free (sigOutTemp);
  return;
}

void
upwlev (double *coefArray, int coefLen, int *waveDecLengthArray,
	int arrayLen, double *lowRe, double *hiRe, int filterLen,
	double *newCoefArray, int newCoefLen, int *newLenArray,
	int newArrayLen, double *approx, int approxLen,
	int stride, extend_method extMethod)
{
  int count, pos1;
  char c='a';
  double *app, *det;

  //printf("enter upwlev!\n");
  for(count=0;count<approxLen;count++)
    approx[count]=coefArray[count];
  for(count=arrayLen-1;count>1;count--)
    newLenArray[count-1]=waveDecLengthArray[count];
  newLenArray[0]=newLenArray[1];

  pos1 = waveDecLengthArray[0] + waveDecLengthArray[1];
  for(count=coefLen-1;count>=pos1;count--)
    newCoefArray[count-coefLen+newCoefLen]=coefArray[count];

  app = malloc(waveDecLengthArray[1]*sizeof(double));
  det = malloc(waveDecLengthArray[1]*sizeof(double));
  for(count=0;count<waveDecLengthArray[1];count++)
    {
      app[count]=coefArray[count];
      det[count]=coefArray[count+waveDecLengthArray[1]];
    }
  idwt_neo (app, det, waveDecLengthArray[1], lowRe, hiRe,
		 filterLen, newCoefArray, waveDecLengthArray[2]);
  free(app);
  free(det);
  return;
}
