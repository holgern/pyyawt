/*
 * -------------------------------------------------------------------------
 * dwt2d.c -- 2-D signal decomposition and reconstruction
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
dwt2D (double *matrixIn, int matrixInRow, int matrixInCol,
       double *matrixOutApprox, double *matrixOutColDetail,
       double *matrixOutRowDetail, double *matrixOutDetail,
       int matrixOutRow, int matrixOutCol, double *lowDe,
       double *hiDe, int filterLen, extend_method extMethod)
{
  int row, col;
  int filterOutLength;
  int matrixOutRowM, matrixOutColM, matrixInRowM, matrixInColM;
  //double temp;
  double *matrixInPre, *matrixInRo;
  double *pMatrixOutApproxPre, *pMatrixOutDetailPre;
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxM, *matrixOutColDetailM;
  double *matrixOutRowDetailM, *matrixOutDetailM;
  char c='b';

  /* extension */
  //matrixInRowM = matrixInRow + 2 * (filterLen - 1);
  //matrixInColM = matrixInCol + 2 * (filterLen - 1);
  matrixInRowM = matrixInRow + 2 * filterLen;
  matrixInColM = matrixInCol + 2 * filterLen;
  if ((extMethod==PER)&&(matrixInRow%2!=0))
     matrixInRowM = matrixInRow + 2 * filterLen + 1;
  if ((extMethod==PER)&&(matrixInCol%2!=0))
     matrixInColM = matrixInCol + 2 * filterLen + 1;

  matrixInRo = malloc (matrixInRowM * matrixInColM * sizeof (double));
  matrixInPre = malloc(matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixIn, matrixInRow, matrixInCol, matrixInRo,
	      matrixInRowM, matrixInColM, extMethod, &c, &c);
  matrix_tran (matrixInRo, matrixInColM, matrixInRowM,
	       matrixInPre, matrixInRowM, matrixInColM);
  free (matrixInRo);

  /* Row DWT  */

  filterOutLength = matrixInColM + filterLen - 1;

  matrixOutColM = filterOutLength/2;
  filterOutLength = matrixInRowM + filterLen - 1;
  //if ((extMethod==PER)&&(matrixInRow%2!=0))
	//  filterOutLength =
  matrixOutRowM = filterOutLength/2;

  matrixOutApproxPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  matrixOutDetailPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  filterOutLength = matrixInColM + filterLen - 1;
  for (row = 0; row < matrixInRowM; row++)
    {
      dwt_no_extension ((matrixInPre + row * matrixInColM), matrixInColM, lowDe,
	   hiDe, filterLen,
	   (matrixOutApproxPre + row * matrixOutColM),
	   (matrixOutDetailPre + row * matrixOutColM),
	       matrixOutColM);//, extMethod);
    }

  free (matrixInPre);

  /* transpose */
  pMatrixOutApproxPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutApproxPre, matrixInRowM, matrixOutColM,
	       pMatrixOutApproxPre, matrixInRowM, matrixOutColM);
  free (matrixOutApproxPre);

  pMatrixOutDetailPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInRowM, matrixOutColM,
	       pMatrixOutDetailPre, matrixInRowM, matrixOutColM);

  free (matrixOutDetailPre);

  /* Col DWT */
  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutApproxM = malloc (matrixOutRowM * matrixOutColM
			     * sizeof (double));
  matrixOutColDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));

  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_no_extension ((pMatrixOutApproxPre + col * matrixInRowM), matrixInRowM,
	   lowDe, hiDe, filterLen,
	   (matrixOutApproxM + col * matrixOutRowM),
	   (matrixOutColDetailM + col * matrixOutRowM),
	       matrixOutRowM);//, extMethod);
    }
  free (pMatrixOutApproxPre);

  wkeep_2D_center (matrixOutApproxM, matrixOutRowM, matrixOutColM,
		   matrixOutApprox, matrixOutRow, matrixOutCol);
  free (matrixOutApproxM);
  wkeep_2D_center (matrixOutColDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutColDetail, matrixOutRow, matrixOutCol);
  free (matrixOutColDetailM);

  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutRowDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));
  matrixOutDetailM = malloc (matrixOutRowM * matrixOutColM * sizeof (double));
  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_no_extension ((pMatrixOutDetailPre + col * matrixInRowM), matrixInRowM,
	   lowDe, hiDe, filterLen,
	   (matrixOutRowDetailM + col * matrixOutRowM),
	   (matrixOutDetailM + col * matrixOutRowM),
	       matrixOutRowM);//, extMethod);
    }
  free (pMatrixOutDetailPre);
  wkeep_2D_center (matrixOutRowDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutRowDetail, matrixOutRow, matrixOutCol);
  free (matrixOutRowDetailM);
  wkeep_2D_center (matrixOutDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutDetail, matrixOutRow, matrixOutCol);
  free (matrixOutDetailM);

  return;
}


void
dwt2D_neo (double *matrixIn, int matrixInRow, int matrixInCol,
       double *matrixOutApprox, double *matrixOutColDetail,
       double *matrixOutRowDetail, double *matrixOutDetail,
       int matrixOutRow, int matrixOutCol, double *lowDe,
       double *hiDe, int filterLen, extend_method extMethod)
{
  int row, col;
  int filterOutLength;
  int matrixOutRowM, matrixOutColM, matrixInRowM, matrixInColM;
  int matrixOutRowT, matrixOutColT;
  double *matrixInPre, *matrixInRo;
  double *pMatrixOutApproxPre, *pMatrixOutDetailPre;
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxM, *matrixOutColDetailM;
  double *matrixOutRowDetailM, *matrixOutDetailM;
  double *matrixOutApproxT, *matrixOutColDetailT;
  double *matrixOutRowDetailT, *matrixOutDetailT;
  char c='b';

  /* extension */
  matrixInRowM = matrixInRow + 2 * filterLen;
  matrixInColM = matrixInCol + 2 * filterLen;
  if ((extMethod==PER)&&(matrixInRow%2!=0))
     matrixInRowM = matrixInRow + 2 * filterLen + 1;
  if ((extMethod==PER)&&(matrixInCol%2!=0))
     matrixInColM = matrixInCol + 2 * filterLen + 1;

  matrixInRo = malloc (matrixInRowM * matrixInColM * sizeof (double));
  matrixInPre = malloc(matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixIn, matrixInRow, matrixInCol, matrixInRo,
	      matrixInRowM, matrixInColM, extMethod, &c, &c);
  matrix_tran (matrixInRo, matrixInColM, matrixInRowM,
	       matrixInPre, matrixInRowM, matrixInColM);
  free (matrixInRo);

  /* Row DWT  */

  filterOutLength = matrixInColM + filterLen - 1;
  matrixOutColM = filterOutLength;
  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutRowM = filterOutLength;

  matrixOutApproxPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  matrixOutDetailPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  filterOutLength = matrixInColM + filterLen - 1;
  for (row = 0; row < matrixInRowM; row++)
    {
      dwt_conv ((matrixInPre + row * matrixInColM), matrixInColM, lowDe,
	   hiDe, filterLen,
	   (matrixOutApproxPre + row * matrixOutColM),
	   (matrixOutDetailPre + row * matrixOutColM),
	       matrixOutColM);
    }

  free (matrixInPre);

  /* transpose */
  pMatrixOutApproxPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutApproxPre, matrixInRowM, matrixOutColM,
	       pMatrixOutApproxPre, matrixInRowM, matrixOutColM);
  free (matrixOutApproxPre);

  pMatrixOutDetailPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInRowM, matrixOutColM,
	       pMatrixOutDetailPre, matrixInRowM, matrixOutColM);

  free (matrixOutDetailPre);

  /* Col DWT */
  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutApproxM = malloc (matrixOutRowM * matrixOutColM
			     * sizeof (double));
  matrixOutColDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));

  //printf("extmethod %d",extMethod);
  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_conv ((pMatrixOutApproxPre + col * matrixInRowM), matrixInRowM,
	   lowDe, hiDe, filterLen,
	   (matrixOutApproxM + col * matrixOutRowM),
	   (matrixOutColDetailM + col * matrixOutRowM),
	       matrixOutRowM);//, extMethod);
    }
  free (pMatrixOutApproxPre);

  matrixOutRowT = matrixInRow + filterLen - 1;
  matrixOutColT = matrixInCol + filterLen - 1;
  if ((extMethod==PER)&&(matrixInRow%2!=0))
     matrixOutRowT = matrixInRow + 1;
  if ((extMethod==PER)&&(matrixInCol%2!=0))
     matrixOutColT = matrixInCol + 1;
  if ((extMethod==PER)&&(matrixInRow%2==0))
     matrixOutRowT = matrixInRow;
  if ((extMethod==PER)&&(matrixInCol%2==0))
     matrixOutColT = matrixInCol;
  matrixOutApproxT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));
  matrixOutColDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));

  wkeep_2D_center (matrixOutApproxM, matrixOutRowM, matrixOutColM,
		   matrixOutApproxT, matrixOutRowT, matrixOutColT);
  free (matrixOutApproxM);
  dyaddown_2D_keep_even (matrixOutApproxT, matrixOutRowT,
		       matrixOutColT, matrixOutApprox,
		       matrixOutRow, matrixOutCol);
  free(matrixOutApproxT);

  wkeep_2D_center (matrixOutColDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutColDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutColDetailM);
  dyaddown_2D_keep_even (matrixOutColDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutColDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutColDetailT);




  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutRowDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));
  matrixOutDetailM = malloc (matrixOutRowM * matrixOutColM * sizeof (double));
  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_conv ((pMatrixOutDetailPre + col * matrixInRowM), matrixInRowM,
	   lowDe, hiDe, filterLen,
	   (matrixOutRowDetailM + col * matrixOutRowM),
	   (matrixOutDetailM + col * matrixOutRowM),
	       matrixOutRowM);
    }
  free (pMatrixOutDetailPre);
  matrixOutRowDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));
  matrixOutDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));


  wkeep_2D_center (matrixOutRowDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutRowDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutRowDetailM);
  dyaddown_2D_keep_even (matrixOutRowDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutRowDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutRowDetailT);

  wkeep_2D_center (matrixOutDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutDetailM);
  dyaddown_2D_keep_even (matrixOutDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutDetailT);

  return;
}


void
idwt2D (double *matrixInApprox, double *matrixInColDetail,
	double *matrixInRowDetail, double *matrixInDetail,
	int matrixInRow, int matrixInCol, double *lowRe,
	double *hiRe, int filterLen, double *matrixOut,
	int matrixOutRow, int matrixOutCol, extend_method extMethod)
{
  int row, col;
  int filterOutLength;
  int matrixInRowM, matrixInColM;
  char c='b';
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxTemp, *matrixOutDetailTemp;
  double *matrixOutPre;
  double *matrixInApproxM, *matrixInColDetailM;
  double *matrixInRowDetailM, *matrixInDetailM;

  matrixInRowM = matrixInRow + 2 * (filterLen - 1);
  matrixInColM = matrixInCol + 2 * (filterLen - 1);

  /* Approx extension */
  matrixInApproxM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixInApprox, matrixInRow, matrixInCol,
	    matrixInApproxM, matrixInRowM, matrixInColM,
	    extMethod, &c, &c);

  matrixInColDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixInColDetail, matrixInRow, matrixInCol,
	    matrixInColDetailM, matrixInRowM, matrixInColM,
	    extMethod, &c, &c);

  matrixInRowDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixInRowDetail, matrixInRow, matrixInCol,
	      matrixInRowDetailM, matrixInRowM, matrixInColM,
	      extMethod, &c, &c);

  matrixInDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));

  wextend_2D (matrixInDetail, matrixInRow, matrixInCol,
	      matrixInDetailM, matrixInRowM, matrixInColM,
	      extMethod, &c, &c);

  filterOutLength = 2 * matrixInRowM + filterLen - 1;

  /* Approx Calculation */
  matrixOutApproxPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrixOutApproxTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInApproxM + col * matrixInRowM),
		   (matrixInColDetailM + col * matrixInRowM),
		   matrixInRowM, lowRe, hiRe, filterLen,
		   (matrixOutApproxPre + col * matrixOutRow),
		   matrixOutRow);
   matrix_tran (matrixOutApproxPre, matrixInColM, matrixOutRow,
	       matrixOutApproxTemp, matrixInColM, matrixOutRow);
  free (matrixOutApproxPre);
  free (matrixInApproxM);
  free (matrixInColDetailM);

  /* Detail Calculation */
  matrixOutDetailPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInRowDetailM + col * matrixInRowM),
		   (matrixInDetailM + col * matrixInRowM),
		   matrixInRowM, lowRe, hiRe, filterLen,
		   (matrixOutDetailPre + col * matrixOutRow),
		   matrixOutRow);
  matrixOutDetailTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInColM, matrixOutRow,
	       matrixOutDetailTemp, matrixInColM, matrixOutRow);
  free (matrixOutDetailPre);
  free (matrixInRowDetailM);
  free (matrixInDetailM);

  /* Final Inverse Transform */
  filterOutLength = 2 * matrixInColM + filterLen - 1;
  matrixOutPre = malloc (matrixOutRow * matrixOutCol * sizeof (double));
  for (row = 0; row < matrixOutRow; row++)
    idwt_neo ((matrixOutApproxTemp + row * matrixInColM),
		   (matrixOutDetailTemp + row * matrixInColM),
		   matrixInColM, lowRe, hiRe, filterLen,
		   (matrixOutPre + row * matrixOutCol),
		   matrixOutCol);
  free (matrixOutApproxTemp);
  free (matrixOutDetailTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixOutCol, matrixOut,
	       matrixOutRow, matrixOutCol);
  free (matrixOutPre);

  return;
}


void
idwt2D_neo (double *matrixInApprox, double *matrixInColDetail,
	double *matrixInRowDetail, double *matrixInDetail,
	int matrixInRow, int matrixInCol, double *lowRe,
	double *hiRe, int filterLen, double *matrixOut,
	int matrixOutRow, int matrixOutCol)
{
  int row, col;
  int filterOutLength;
  int matrixInRowM, matrixInColM;
  char c='b';
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxTemp, *matrixOutDetailTemp;
  double *matrixOutPre;
  //double *matrixInApproxM, *matrixInColDetailM;
  //double *matrixInRowDetailM, *matrixInDetailM;

  matrixInRowM = matrixInRow;// + 2 * (filterLen - 1);
  matrixInColM = matrixInCol;// + 2 * (filterLen - 1);

  /* Approx extension */
  //matrixInApproxM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  //wextend_2D (matrixInApprox, matrixInRow, matrixInCol,
	//    matrixInApproxM, matrixInRowM, matrixInColM,
	  //  extMethod, &c, &c);

  //matrixInColDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  //wextend_2D (matrixInColDetail, matrixInRow, matrixInCol,
	//    matrixInColDetailM, matrixInRowM, matrixInColM,
	  //  extMethod, &c, &c);

  //matrixInRowDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));
  //wextend_2D (matrixInRowDetail, matrixInRow, matrixInCol,
	//      matrixInRowDetailM, matrixInRowM, matrixInColM,
	  //    extMethod, &c, &c);

  //matrixInDetailM = malloc (matrixInRowM * matrixInColM * sizeof (double));

  //wextend_2D (matrixInDetail, matrixInRow, matrixInCol,
	//      matrixInDetailM, matrixInRowM, matrixInColM,
	  //    extMethod, &c, &c);

  filterOutLength = 2 * matrixInRowM + filterLen - 1;

  /* Approx Calculation */
  matrixOutApproxPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrixOutApproxTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInApprox + col * matrixInRowM),
		   (matrixInColDetail + col * matrixInRowM),
		   matrixInRowM, lowRe, hiRe, filterLen,
		   (matrixOutApproxPre + col * matrixOutRow),
		   matrixOutRow);
   matrix_tran (matrixOutApproxPre, matrixInColM, matrixOutRow,
	       matrixOutApproxTemp, matrixInColM, matrixOutRow);
  free (matrixOutApproxPre);
  //free (matrixInApproxM);
  //free (matrixInColDetailM);

  /* Detail Calculation */
  matrixOutDetailPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInRowDetail + col * matrixInRowM),
		   (matrixInDetail + col * matrixInRowM),
		   matrixInRowM, lowRe, hiRe, filterLen,
		   (matrixOutDetailPre + col * matrixOutRow),
		   matrixOutRow);
  matrixOutDetailTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInColM, matrixOutRow,
	       matrixOutDetailTemp, matrixInColM, matrixOutRow);
  free (matrixOutDetailPre);
  //free (matrixInRowDetailM);
  //free (matrixInDetailM);

  /* Final Inverse Transform */
  filterOutLength = 2 * matrixInColM + filterLen - 1;
  matrixOutPre = malloc (matrixOutRow * matrixOutCol * sizeof (double));
  for (row = 0; row < matrixOutRow; row++)
    idwt_neo ((matrixOutApproxTemp + row * matrixInColM),
		   (matrixOutDetailTemp + row * matrixInColM),
		   matrixInColM, lowRe, hiRe, filterLen,
		   (matrixOutPre + row * matrixOutCol),
		   matrixOutCol);
  free (matrixOutApproxTemp);
  free (matrixOutDetailTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixOutCol, matrixOut,
	       matrixOutRow, matrixOutCol);
  free (matrixOutPre);

  return;
}


void
wave_mem_cal (int *pLen, int stride, int *total)
{
  int count = 0;
  (*total) = 4 * (pLen[2]) * (pLen[3]);
  for (count = 2; count < (stride+1); count++)
    (*total) += 3 * (pLen[2 * count]) * (pLen[2 * count + 1]);
  return;
}

void
matrix_wavedec_len_cal (int matrixInRow, int matrixInCol, int stride,
			int filterLen, int *pLen)
{
  int row;
  int filterOutLength;

  pLen[(stride+1) * 2] = matrixInRow;
  pLen[(stride+1) * 2 + 1] = matrixInCol;

  if (getdwtMode()!=PER)
  {
     for (row = stride; row > 0; row--)
       {
         filterOutLength = pLen[(row + 1) * 2] + filterLen - 1;
         pLen[row*2] = filterOutLength/2;
         filterOutLength = pLen[(row + 1) * 2 + 1] + filterLen - 1;
         pLen[row*2+1] = filterOutLength/2;
       }
     pLen[0] = pLen[2];
     pLen[1] = pLen[3];
  }
  else
  {
     for (row = stride; row > 0; row--)
       {
         //filterOutLength = pLen[(row + 1) * 2] + filterLen - 1;
         pLen[row*2] = (int)ceil(((double)(pLen[(row + 1) * 2])) / 2.0);
         //filterOutLength = pLen[(row + 1) * 2 + 1] + filterLen - 1;
         pLen[row*2+1] = (int)ceil(((double)(pLen[(row + 1) * 2 + 1])) / 2.0);
       }
     pLen[0] = pLen[2];
     pLen[1] = pLen[3];

  }
  return;
}

void
matrix_locate (int stride, int *pLen, int *pH, int *pV, int *pD)
{
  int count;
  int area, acre;

  /* *pH = (*pLen)*(*(pLen+1)); */
  pH[0] = pLen[0] * pLen[1];
  /* *pV = 2 * (*pH); */
  pV[0] = 2 * pH[0];
  /* *pD = 3 * (*pH); */
  pD[0] = 3 * pH[0];

  for (count = 1; count < stride; count++)
    {
      /* area = (*(pLen+(count-1)*2)) * (*(pLen+(count-1)*2+1)); */
      area = pLen[2 * count] * pLen[2 * count + 1];
      /* acre = (*(pLen+count*2)) * (*(pLen+count*2+1)); */
      acre = pLen[2 * (count + 1)] * pLen[2 * (count + 1) + 1];
      /* *(pH+count) = *(pH+count-1) + 2 * area; */
      pH[count] = pH[count - 1] + 3 * area;
      /* *(pV+count) = *(pV+count-1) + area + acre; */
      pV[count] = pV[count - 1] + 2 * area + acre;
      /* *(pD+count) = *(pD+count-1) + 2 * acre; */
      pD[count] = pD[count - 1] + area + 2 * acre;
    }
  /* *total = *(pD+stride-1)+ acre - 1; */

  return;
}


void
wavedec2 (double *matrixIn, int matrixInRow, int matrixInCol,
	  double *lowDe, double *hiDe, int filterLen,
	  int *pLen, double *coef, int sigOutLength, int stride,
	  extend_method extMethod)
{
  int count, row, col;
  int *pH, *pV, *pD;
  double *matrixInTemp;
  double *matrixOutApproxTemp;
  //printf("waverec2 ok1");
  matrixInTemp = malloc ((pLen[(stride+1) * 2]) *(pLen[(stride+1) * 2 + 1]) * sizeof (double));
  matrixOutApproxTemp = malloc ((pLen[stride * 2]) *(pLen[stride * 2 + 1]) * sizeof (double));
  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for (row = 0; row < pLen[(stride+1) * 2]; row++)
    {
      for (col = 0; col < pLen[(stride+1) * 2 + 1]; col++)
	matrixInTemp[col + row * (pLen[(stride+1) * 2 + 1])] =
	  matrixIn[col + row * (pLen[(stride+1) * 2 + 1])];
    }
    //printf("waverec2 ok2");
  for (count = stride - 1; count >= 0; count--)
    {
      dwt2D_neo (matrixInTemp, pLen[2 * (count + 2)],
	     pLen[2 * (count + 2) + 1], matrixOutApproxTemp,
	     (pH[count] + coef), (pV[count] + coef),
	     (pD[count] + coef), pLen[2 * (count + 1)],
	     pLen[2 * (count + 1) + 1], lowDe,
	     hiDe, filterLen, extMethod);
      //swtdwt2Dex (matrixInTemp, pLen[2 * (count + 1)],
      //	  pLen[2 * (count + 1) + 1], matrixOutApproxTemp,
      //	  (pH[count] + coef), (pV[count] + coef),
      //	  (pD[count] + coef), waveType, mode);
      for (row = 0; row < pLen[2 * (count+1)]; row++)
	{
	  for (col = 0; col < pLen[2 * (count+1) + 1]; col++)
	    matrixInTemp[col + row * (pLen[2 * (count+1) + 1])] =
	      matrixOutApproxTemp[col + row * (pLen[2 * (count+1) + 1])];
       // printf("");
	}
    }
  free (matrixInTemp);
  free (pH);
  free (pV);
  free (pD);
 // printf("waverec2 ok3");
  for (row = 0; row < pLen[0]; row++)
    {
      for (col = 0; col < pLen[1]; col++)
	coef[col + row * (pLen[1])] =
	  matrixOutApproxTemp[col + row * (pLen[1])];
    }
  free (matrixOutApproxTemp);
  return;
}

void
waverec2 (double *coef, int sigInLength, double *lowRe, double *hiRe,
	  int filterLen, double *matrixOut, int matrixOutRow,
	  int matrixOutCol, int *pLen, int stride,
	  extend_method extMethod)
{
  int count, row, col;
  int *pH, *pV, *pD;
  double *matrixOutTemp;
  double *matrixInApproxTemp;

  matrixOutTemp = malloc ((pLen[(stride+1) * 2]) *
			  (pLen[(stride+1) * 2 + 1]) * sizeof (double));
  matrixInApproxTemp = malloc ((pLen[(stride+1) * 2]) *
			       (pLen[(stride+1) * 2 + 1]) * sizeof (double));

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for (row = 0; row < pLen[0]; row++)
    {
      for (col = 0; col < pLen[1]; col++)
	matrixInApproxTemp[col + row * (pLen[1])] =
	  coef[col + row * (pLen[1])];
    }

  for (count = 0; count < stride; count++)
    {
      //idwt2D (matrixInApproxTemp, (coef + pH[count]),
	    //  (coef + pV[count]), (coef + pD[count]),
	      //pLen[(count+1) * 2], pLen[(count+1) * 2 + 1],
	      //lowRe, hiRe, filterLen, matrixOutTemp,
	      //pLen[(count + 2) * 2], pLen[(count + 2) * 2 + 1],
	      //extMethod);
	  idwt2D_neo (matrixInApproxTemp, (coef + pH[count]),
	      (coef + pV[count]), (coef + pD[count]),
	      pLen[(count+1) * 2], pLen[(count+1) * 2 + 1],
	      lowRe, hiRe, filterLen, matrixOutTemp,
	      pLen[(count + 2) * 2], pLen[(count + 2) * 2 + 1]);
      //swtidwt2Dex (matrixInApproxTemp, (coef + pH[count]),
      //	   (coef + pV[count]), (coef + pD[count]),
      //	   pLen[count * 2], pLen[count * 2 + 1], matrixOutTemp,
      //	   pLen[(count + 1) * 2], pLen[(count + 1) * 2 + 1], waveType,
      //	   mode);
      for (row = 0; row < pLen[(count + 2) * 2]; row++)
	{
	  for (col = 0; col < pLen[(count + 2) * 2 + 1]; col++)
	    matrixInApproxTemp[col + row * (pLen[(count + 2) * 2 + 1])] =
	      matrixOutTemp[col + row * (pLen[(count + 2) * 2 + 1])];
	}
    }

  for (row = 0; row < pLen[(stride+1) * 2]; row++)
    {
      for (col = 0; col < pLen[(stride+1) * 2 + 1]; col++)
	matrixOut[col + row * (pLen[(stride+1) * 2 + 1])] =
	  matrixInApproxTemp[col + row * (pLen[(stride+1) * 2 + 1])];
    }


  free (pH);
  free (pV);
  free (pD);
  free (matrixInApproxTemp);
  free (matrixOutTemp);
  return;
}

void
wenergy_2output (double *coef, int sigInLength, int *pLen,
		 double *ae, double *de, int deLength, int stride)
{
  int count, countt, startp, endp;
  int *pH, *pV, *pD;
  double ene;

  ene = 0;
  for(count=0;count<sigInLength;count++)
    ene += coef[count]*coef[count];

  *ae = 0;
  for(count=0;count<pLen[0]*pLen[1];count++)
    *ae += coef[count]*coef[count];
  *ae = (*ae)*100/ene;

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for(count=0;count<stride;count++)
    {
      startp = pH[count];
      endp = pH[count] + 3 * pLen[(count+1)*2] * pLen[(count+1)*2+1];
      de[count] = 0;
      for(countt=startp;countt<endp;countt++)
	de[count] += coef[countt]*coef[countt];
      de[count] = de[count]*100/ene;
    }
  free(pH);
  free(pV);
  free(pD);
  return;
}


void
wenergy_4output (double *coef, int sigInLength, int *pLen,
		 double *ae, double *he, double *ve, double *de,
		 int deLength, int stride)
{
  int count, countt, *pH, *pV, *pD;
  int starth, startv, startd, endh, endv, endd;
  double ene;
  ene = 0;
  for(count=0;count<sigInLength;count++)
    ene += coef[count]*coef[count];

  *ae = 0;
  for(count=0;count<pLen[0]*pLen[1];count++)
    *ae += coef[count]*coef[count];
  *ae = (*ae)*100/ene;

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for(count=0;count<stride;count++)
    {
      /* horizontal */
      starth = pH[count];
      endh = pH[count] + pLen[(count+1)*2] * pLen[(count+1)*2+1];
      he[count] = 0;
      for(countt=starth;countt<endh;countt++)
	he[count] += coef[countt]*coef[countt];
      he[count] = he[count]*100/ene;
      /* vertical */
      startv = pV[count];
      endv = pV[count] + pLen[(count+1)*2] * pLen[(count+1)*2+1];
      ve[count] = 0;
      for(countt=startv;countt<endv;countt++)
	ve[count] += coef[countt]*coef[countt];
      ve[count] = ve[count]*100/ene;
      /* detail */
      startd = pD[count];
      endd = pD[count] + pLen[(count+1)*2] * pLen[(count+1)*2+1];
      de[count] = 0;
      for(countt=startd;countt<endd;countt++)
	de[count] += coef[countt]*coef[countt];
      de[count] = de[count]*100/ene;
    }
  free(pH);
  free(pV);
  free(pD);

  return;
}


void
detcoef2 (double *coef, int sigInLength, double *coefOut,
	  int sigOutLength, int *pLen, int stride, int level,
	  char *coefType)
{
  int row, col, sta;
  int *pH, *pV, *pD;

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);
  if (strcmp (coefType, "h") == 0)
    {
      sta = pH[stride - level];
    }
  if (strcmp (coefType, "v") == 0)
    {
      sta = pV[stride - level];
    }
  if (strcmp (coefType, "d") == 0)
    {
      sta = pD[stride - level];
    }

  for (row = 0; row < pLen[(stride - level + 1) * 2]; row++)
    {
      for (col = 0; col < pLen[(stride - level + 1) * 2 + 1]; col++)
	coefOut[col + row * (pLen[(stride - level + 1) * 2 + 1])] =
	  coef[sta + col + row * (pLen[(stride - level + 1) * 2 + 1])];
    }
  free (pH);
  free (pV);
  free (pD);

  return;
}


void
appcoef2 (double *coef, int sigInLength, double *lowRe,
	  double *hiRe, int filterLen, double *coefOut,
	  int matrixOutRow, int matrixOutCol, int *pLen,
	  int stride, int level, extend_method extMethod)
{
  int count;

  if (level == stride)
    {
      for (count = 0; count < ((pLen[0]) * (pLen[1])); count++)
	coefOut[count] = coef[count];
    }
  else
    {
      waverec2 (coef, sigInLength, lowRe, hiRe, filterLen,
		coefOut, matrixOutRow, matrixOutCol,
		pLen, stride-level, extMethod);
    }

  return;
}


void
wrcoef2 (double *coef, int sigInLength, double *lowRe,
	 double *hiRe, int filterLen, double *matrixOut,
	 int matrixOutRow, int matrixOutCol, int *pLen,
	 int stride, int level, char *type, extend_method extMethod)
{
  int count, total, sta, si;
  double *coefTemp;
  int *pH, *pV, *pD;

  wave_mem_cal (pLen, stride, &total);
  coefTemp = malloc (total * sizeof (double));

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for (count = 0; count < total; count++)
    coefTemp[count] = 0;

  if (strcmp (type, "h") == 0)
    {
      sta = pH[stride - level];
      si = (pLen[(stride - level + 1) * 2]) * (pLen[(stride - level + 1) * 2 + 1]);
    }
  if (strcmp (type, "v") == 0)
    {
      sta = pV[stride - level];
      si = (pLen[(stride - level + 1) * 2]) * (pLen[(stride - level + 1) * 2 + 1]);
    }
  if (strcmp (type, "d") == 0)
    {
      sta = pD[stride - level];
      si = (pLen[(stride - level + 1) * 2]) * (pLen[(stride - level + 1) * 2 + 1]);
    }
  if (strcmp (type, "a") == 0)
    {
      sta = 0;
      si = (pLen[0]) * (pLen[1]);
      if (level != stride)
	{
	  for (count = 1; count <= (stride - level); count++)
	    si += 3 * (pLen[(count) * 2]) * (pLen[(count) * 2 + 1]);
	}
    }

  for (count = sta; count < (si + sta); count++)
    coefTemp[count] = coef[count];

  waverec2 (coefTemp, sigInLength, lowRe, hiRe, filterLen,
	    matrixOut, matrixOutRow, matrixOutCol,
	    pLen, stride, extMethod);

  //waverec2 (coefTemp, pLen, stride, waveType, matrixOut, mode);

  free (pH);
  free (pV);
  free (pD);
  free (coefTemp);
  return;
}


void
upwlev2 (double *coef, int sigInLength, double *lowRe, double *hiRe,
	 int filterLen, int *pLen, int matrixRow, int matrixCol,
	 double *approx, int approxLen, double *newCoef,
	 int newCoefLen, int *newLenMatrix, int lenMatrixRow,
	 int lenMatrixCol, int stride, extend_method extMethod)
{
  int count, *le, *pH, *pV, *pD, row, col, pos;

  for(count=0;count<approxLen;count++)
    approx[count]=coef[count];

  le = malloc((matrixRow-1)*matrixCol*sizeof(int));
  for(count=stride+1;count>1;count--)
    {
      le[(count-1)*2] = pLen[count*2];
      le[(count-1)*2+1] = pLen[count*2+1];
    }
  le[0] = pLen[4];
  le[1] = pLen[5];

  //for(count=0;count<=stride;count++)
  //sciprint("%d %d\n",le[count*2],le[count*2+1]);
  for (col = 0; col < matrixCol; col++)
    {
      for (row = 0; row < (matrixRow-1); row++)
	{
	  newLenMatrix[row + col * (matrixRow-1)] =
	    le[col + row * matrixCol];
	}
   }

  //for (row = 0; row < matrixCol; row++)
  //{
  //for (col = 0; col < (matrixRow-1); col++)
  //  newLenMatrix[row + col * matrixCol] = le[col + row * matrixRow];
  //}
 //matrix_tran(le,matrixRow-1,matrixCol,newLenMatrix,
  //      matrixRow-1,matrixCol);
  free(le);

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  pos = 4*pH[0];
  for(count=sigInLength-1;count>=pos;count--)
    newCoef[count-sigInLength+newCoefLen] = coef[count];


  //idwt2D (coef, coef+pH[0], coef+pV[0], coef+pD[0],
	//  pLen[0], pLen[1], lowRe, hiRe, filterLen, newCoef,
	  //pLen[4], pLen[5], extMethod);
  idwt2D_neo (coef, coef+pH[0], coef+pV[0], coef+pD[0],
	  pLen[0], pLen[1], lowRe, hiRe, filterLen, newCoef,
	  pLen[4], pLen[5]);


  free(pH);
  free(pV);
  free(pD);
  return;
}


void
upcoef2 (double *matrixIn, int matrixInRow, int matrixInCol,
	 double *lowRe, double *hiRe, int filterLen,
	 double *matrixOut, int matrixOutRow, int matrixOutCol,
	 int matrixOutDefaultRow, int matrixOutDefaultCol,
	 int step, char *type)//, extend_method extMethod)
{
  double *vo, *matrixOutTemp, *matrixOutPre;
  int matrixOutTempRow, matrixOutTempCol, rowLeng, colLeng,count, count1;

//   matrixOutTempRow = 2*matrixInRow - filterLen + 2;
//   matrixOutTempCol = 2*matrixInCol - filterLen + 2;
  matrixOutTempRow = 2*matrixInRow + filterLen - 2;
  matrixOutTempCol = 2*matrixInCol + filterLen - 2;

  vo = malloc(matrixInRow*matrixInCol*sizeof(double));
  for(count = 0;count<matrixInRow*matrixInCol;count++)
    vo[count] = 0;
  matrixOutTemp = malloc(matrixOutDefaultRow*
			 matrixOutDefaultCol*sizeof(double));
  for(count=0;count<matrixOutDefaultRow*matrixOutDefaultCol;count++)
    matrixOutTemp[count] = 0;

  if (!strcmp(type,"a"))
    {
      //idwt2D (matrixIn, vo, vo, vo, matrixInRow, matrixInCol,
	    //  lowRe, hiRe, filterLen, matrixOutTemp,
	      //matrixOutTempRow, matrixOutTempCol, extMethod);
	  idwt2D_neo (matrixIn, vo, vo, vo, matrixInRow, matrixInCol,
	      lowRe, hiRe, filterLen, matrixOutTemp,
	      matrixOutTempRow, matrixOutTempCol);
    }

  if (!strcmp(type,"h"))
    {
      //idwt2D (vo, matrixIn, vo, vo, matrixInRow, matrixInCol,
	    //  lowRe, hiRe, filterLen, matrixOutTemp,
	      //matrixOutTempRow, matrixOutTempCol, extMethod);
	  idwt2D_neo (vo, matrixIn, vo, vo, matrixInRow, matrixInCol,
	      lowRe, hiRe, filterLen, matrixOutTemp,
	      matrixOutTempRow, matrixOutTempCol);
    }

  if (!strcmp(type,"v"))
    {
      //idwt2D (vo, vo, matrixIn, vo, matrixInRow, matrixInCol,
	    //  lowRe, hiRe, filterLen, matrixOutTemp,
	      //matrixOutTempRow, matrixOutTempCol, extMethod);
	  idwt2D_neo (vo, vo, matrixIn, vo, matrixInRow, matrixInCol,
	      lowRe, hiRe, filterLen, matrixOutTemp,
	      matrixOutTempRow, matrixOutTempCol);
    }

  if (!strcmp(type,"d"))
    {
      //idwt2D (vo, vo, vo, matrixIn, matrixInRow, matrixInCol,
	    //  lowRe, hiRe, filterLen, matrixOutTemp,
	      //matrixOutTempRow, matrixOutTempCol, extMethod);
	  idwt2D_neo (vo, vo, vo, matrixIn, matrixInRow, matrixInCol,
	      lowRe, hiRe, filterLen, matrixOutTemp,
	      matrixOutTempRow, matrixOutTempCol);
    }

  free(vo);

  if (step>1)
    {
      matrixOutPre = malloc(matrixOutDefaultRow*
			    matrixOutDefaultCol*sizeof(double));
      for(count=0;count<matrixOutDefaultRow*
	    matrixOutDefaultCol;count++)
	matrixOutPre[count] = 0;
      rowLeng = matrixOutTempRow;
      colLeng = matrixOutTempCol;
      for(count=0;count<step-1;count++)
	{
	  //idwt2D (matrixOutTemp, vo, vo, vo, rowLeng, colLeng,
		//  lowRe, hiRe, filterLen, matrixOutPre,
		  //rowLeng*2-filterLen+2, colLeng*2-filterLen+2,
		  //extMethod);
      vo = malloc(rowLeng*colLeng*sizeof(double));
      for(count1 = 0;count1<rowLeng*colLeng;count1++)
         vo[count1] = 0;
	  idwt2D_neo (matrixOutTemp, vo, vo, vo, rowLeng, colLeng,
		  lowRe, hiRe, filterLen, matrixOutPre,
		  rowLeng*2+filterLen-2, colLeng*2+filterLen-2);
	  rowLeng = rowLeng*2 + filterLen - 2;
	  colLeng = colLeng*2 + filterLen - 2;
// 	  idwt2D_neo (matrixOutTemp, vo, vo, vo, rowLeng, colLeng,
// 		  lowRe, hiRe, filterLen, matrixOutPre,
// 		  rowLeng*2-filterLen+2, colLeng*2-filterLen+2);
// 	  rowLeng = rowLeng*2 - filterLen + 2;
// 	  colLeng = colLeng*2 - filterLen + 2;
	  verbatim_copy (matrixOutPre, rowLeng*colLeng,
			 matrixOutTemp, rowLeng*colLeng);
	  free(vo);
	}
      matrixOutTempRow = rowLeng;
      matrixOutTempCol = colLeng;
      free(matrixOutPre);
    }

  wkeep_2D_center (matrixOutTemp, matrixOutDefaultRow,
		   matrixOutDefaultCol, matrixOut,
		   matrixOutRow, matrixOutCol);


  //free(vo);
  free(matrixOutTemp);

  return;
}

void
wavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
	   double *lowDeR, double *hiDeR, double *lowDeC,
	   double *hiDeC, int filterLen, int *pLen,
	   double *coef, int sigOutLength, int stride,
	   extend_method extMethod)
{
  int count, row, col;
  int *pH, *pV, *pD;
  double *matrixInTemp;
  double *matrixOutApproxTemp;

  matrixInTemp = malloc ((pLen[(stride+1) * 2]) *
			 (pLen[(stride+1) * 2 + 1]) * sizeof (double));
  matrixOutApproxTemp = malloc ((pLen[stride * 2]) *
				(pLen[stride * 2 + 1]) *
				sizeof (double));
  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for (row = 0; row < pLen[(stride+1) * 2]; row++)
    {
      for (col = 0; col < pLen[(stride+1) * 2 + 1]; col++)
	matrixInTemp[col + row * (pLen[(stride+1) * 2 + 1])] =
	  matrixIn[col + row * (pLen[(stride+1) * 2 + 1])];
    }

  for (count = stride - 1; count >= 0; count--)
    {
      dwt2D_neo_a (matrixInTemp, pLen[2 * (count + 2)],
		 pLen[2 * (count + 2) + 1], matrixOutApproxTemp,
		 (pH[count] + coef), (pV[count] + coef),
		 (pD[count] + coef), pLen[2 * (count + 1)],
		 pLen[2 * (count + 1) + 1], lowDeR,
		 hiDeR, lowDeC, hiDeC, filterLen, extMethod);
      for (row = 0; row < pLen[2 * (count+1)]; row++)
	{
	  for (col = 0; col < pLen[2 * (count+1) + 1]; col++)
	    matrixInTemp[col + row * (pLen[2 * (count+1) + 1])] =
	      matrixOutApproxTemp[col + row * (pLen[2 * (count+1) + 1])];
	}
    }
  free (matrixInTemp);
  free (pH);
  free (pV);
  free (pD);
  for (row = 0; row < pLen[0]; row++)
    {
      for (col = 0; col < pLen[1]; col++)
	coef[col + row * (pLen[1])] =
	  matrixOutApproxTemp[col + row * (pLen[1])];
    }
  free (matrixOutApproxTemp);
  return;
}

void
dwt2D_neo_a (double *matrixIn, int matrixInRow, int matrixInCol,
	     double *matrixOutApprox, double *matrixOutColDetail,
	     double *matrixOutRowDetail, double *matrixOutDetail,
	     int matrixOutRow, int matrixOutCol, double *lowDeR,
	     double *hiDeR, double *lowDeC, double *hiDeC,
	     int filterLen, extend_method extMethod)
{
  int row, col;
  int filterOutLength;
  int matrixOutRowM, matrixOutColM, matrixInRowM, matrixInColM;
  int matrixOutRowT, matrixOutColT;
  double *matrixInPre, *matrixInRo;
  double *pMatrixOutApproxPre, *pMatrixOutDetailPre;
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxM, *matrixOutColDetailM;
  double *matrixOutRowDetailM, *matrixOutDetailM;
  double *matrixOutApproxT, *matrixOutColDetailT;
  double *matrixOutRowDetailT, *matrixOutDetailT;
  char c='b';

  /* extension */
  matrixInRowM = matrixInRow + 2 * filterLen;
  matrixInColM = matrixInCol + 2 * filterLen;
  if ((extMethod==PER)&&(matrixInRow%2!=0))
     matrixInRowM = matrixInRow + 2 * filterLen + 1;
  if ((extMethod==PER)&&(matrixInCol%2!=0))
     matrixInColM = matrixInCol + 2 * filterLen + 1;

  matrixInRo = malloc (matrixInRowM * matrixInColM * sizeof (double));
  matrixInPre = malloc(matrixInRowM * matrixInColM * sizeof (double));
  wextend_2D (matrixIn, matrixInRow, matrixInCol, matrixInRo,
	      matrixInRowM, matrixInColM, extMethod, &c, &c);
  matrix_tran (matrixInRo, matrixInColM, matrixInRowM,
	       matrixInPre, matrixInRowM, matrixInColM);
  free (matrixInRo);

  /* Row DWT  */

  filterOutLength = matrixInColM + filterLen - 1;
  matrixOutColM = filterOutLength;
  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutRowM = filterOutLength;

  matrixOutApproxPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  matrixOutDetailPre = malloc (matrixInRowM *
			       matrixOutColM * sizeof (double));
  filterOutLength = matrixInColM + filterLen - 1;
  for (row = 0; row < matrixInRowM; row++)
    {
      dwt_conv ((matrixInPre + row * matrixInColM), matrixInColM, lowDeR,
	   hiDeR, filterLen,
	   (matrixOutApproxPre + row * matrixOutColM),
	   (matrixOutDetailPre + row * matrixOutColM),
	       matrixOutColM);
    }

  free (matrixInPre);

  /* transpose */
  pMatrixOutApproxPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutApproxPre, matrixInRowM, matrixOutColM,
	       pMatrixOutApproxPre, matrixInRowM, matrixOutColM);
  free (matrixOutApproxPre);

  pMatrixOutDetailPre = malloc (matrixInRowM *
				matrixOutColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInRowM, matrixOutColM,
	       pMatrixOutDetailPre, matrixInRowM, matrixOutColM);

  free (matrixOutDetailPre);

  /* Col DWT */
  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutApproxM = malloc (matrixOutRowM * matrixOutColM
			     * sizeof (double));
  matrixOutColDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));

  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_conv ((pMatrixOutApproxPre + col * matrixInRowM), matrixInRowM,
	   lowDeC, hiDeC, filterLen,
	   (matrixOutApproxM + col * matrixOutRowM),
	   (matrixOutColDetailM + col * matrixOutRowM),
	       matrixOutRowM);//, extMethod);
    }
  free (pMatrixOutApproxPre);

  matrixOutRowT = matrixInRow + filterLen - 1;
  matrixOutColT = matrixInCol + filterLen - 1;
  if ((extMethod==PER)&&(matrixInRow%2!=0))
     matrixOutRowT = matrixInRow + 1;
  if ((extMethod==PER)&&(matrixInCol%2!=0))
     matrixOutColT = matrixInCol + 1;
  if ((extMethod==PER)&&(matrixInRow%2==0))
     matrixOutRowT = matrixInRow;
  if ((extMethod==PER)&&(matrixInCol%2==0))
     matrixOutColT = matrixInCol;
  matrixOutApproxT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));
  matrixOutColDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));

  wkeep_2D_center (matrixOutApproxM, matrixOutRowM, matrixOutColM,
		   matrixOutApproxT, matrixOutRowT, matrixOutColT);
  free (matrixOutApproxM);
  dyaddown_2D_keep_even (matrixOutApproxT, matrixOutRowT,
		       matrixOutColT, matrixOutApprox,
		       matrixOutRow, matrixOutCol);
  free(matrixOutApproxT);

  wkeep_2D_center (matrixOutColDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutColDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutColDetailM);
  dyaddown_2D_keep_even (matrixOutColDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutColDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutColDetailT);




  filterOutLength = matrixInRowM + filterLen - 1;
  matrixOutRowDetailM = malloc (matrixOutRowM * matrixOutColM
				* sizeof (double));
  matrixOutDetailM = malloc (matrixOutRowM * matrixOutColM * sizeof (double));
  for (col = 0; col < matrixOutColM; col++)
    {
      dwt_conv ((pMatrixOutDetailPre + col * matrixInRowM), matrixInRowM,
	   lowDeC, hiDeC, filterLen,
	   (matrixOutRowDetailM + col * matrixOutRowM),
	   (matrixOutDetailM + col * matrixOutRowM),
	       matrixOutRowM);
    }
  free (pMatrixOutDetailPre);
  matrixOutRowDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));
  matrixOutDetailT = malloc (matrixOutRowT*matrixOutColT*sizeof(double));


  wkeep_2D_center (matrixOutRowDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutRowDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutRowDetailM);
  dyaddown_2D_keep_even (matrixOutRowDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutRowDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutRowDetailT);

  wkeep_2D_center (matrixOutDetailM, matrixOutRowM, matrixOutColM,
		   matrixOutDetailT, matrixOutRowT, matrixOutColT);
  free (matrixOutDetailM);
  dyaddown_2D_keep_even (matrixOutDetailT, matrixOutRowT,
		       matrixOutColT, matrixOutDetail,
		       matrixOutRow, matrixOutCol);
  free(matrixOutDetailT);

  return;
}

void
waverec2a (double *coef, int sigInLength, double *lowReR,
	   double *hiReR, double *lowReC, double *hiReC,
	  int filterLen, double *matrixOut, int matrixOutRow,
	  int matrixOutCol, int *pLen, int stride,
	  extend_method extMethod)
{
  int count, row, col;
  int *pH, *pV, *pD;
  double *matrixOutTemp;
  double *matrixInApproxTemp;

  matrixOutTemp = malloc ((pLen[(stride+1) * 2]) *
			  (pLen[(stride+1) * 2 + 1]) * sizeof (double));
  matrixInApproxTemp = malloc ((pLen[(stride+1) * 2]) *
			       (pLen[(stride+1) * 2 + 1]) * sizeof (double));

  pH = malloc (stride * sizeof (int));
  pV = malloc (stride * sizeof (int));
  pD = malloc (stride * sizeof (int));
  matrix_locate (stride, pLen, pH, pV, pD);

  for (row = 0; row < pLen[0]; row++)
    {
      for (col = 0; col < pLen[1]; col++)
	matrixInApproxTemp[col + row * (pLen[1])] =
	  coef[col + row * (pLen[1])];
    }

  for (count = 0; count < stride; count++)
    {
	  idwt2D_neo_a (matrixInApproxTemp, (coef + pH[count]),
			(coef + pV[count]), (coef + pD[count]),
			pLen[(count+1) * 2], pLen[(count+1) * 2 + 1],
			lowReR, hiReR, lowReC, hiReC, filterLen,
			matrixOutTemp, pLen[(count + 2) * 2],
			pLen[(count + 2) * 2 + 1]);

      for (row = 0; row < pLen[(count + 2) * 2]; row++)
	{
	  for (col = 0; col < pLen[(count + 2) * 2 + 1]; col++)
	    matrixInApproxTemp[col + row * (pLen[(count + 2) * 2 + 1])] =
	      matrixOutTemp[col + row * (pLen[(count + 2) * 2 + 1])];
	}
    }

  for (row = 0; row < pLen[(stride+1) * 2]; row++)
    {
      for (col = 0; col < pLen[(stride+1) * 2 + 1]; col++)
	matrixOut[col + row * (pLen[(stride+1) * 2 + 1])] =
	  matrixInApproxTemp[col + row * (pLen[(stride+1) * 2 + 1])];
    }


  free (pH);
  free (pV);
  free (pD);
  free (matrixInApproxTemp);
  free (matrixOutTemp);
  return;
}

void
idwt2D_neo_a (double *matrixInApprox, double *matrixInColDetail,
	      double *matrixInRowDetail, double *matrixInDetail,
	      int matrixInRow, int matrixInCol, double *lowReR,
	      double *hiReR, double *lowReC, double *hiReC,
	      int filterLen, double *matrixOut,
	      int matrixOutRow, int matrixOutCol)
{
  int row, col;
  int filterOutLength;
  int matrixInRowM, matrixInColM;
  char c='b';
  double *matrixOutApproxPre, *matrixOutDetailPre;
  double *matrixOutApproxTemp, *matrixOutDetailTemp;
  double *matrixOutPre;

  matrixInRowM = matrixInRow;// + 2 * (filterLen - 1);
  matrixInColM = matrixInCol;// + 2 * (filterLen - 1);


  filterOutLength = 2 * matrixInRowM + filterLen - 1;

  /* Approx Calculation */
  matrixOutApproxPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrixOutApproxTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInApprox + col * matrixInRowM),
		   (matrixInColDetail + col * matrixInRowM),
		   matrixInRowM, lowReC, hiReC, filterLen,
		   (matrixOutApproxPre + col * matrixOutRow),
		   matrixOutRow);
   matrix_tran (matrixOutApproxPre, matrixInColM, matrixOutRow,
	       matrixOutApproxTemp, matrixInColM, matrixOutRow);
  free (matrixOutApproxPre);

  /* Detail Calculation */
  matrixOutDetailPre = malloc (matrixOutRow * matrixInColM * sizeof (double));
  for (col = 0; col < matrixInColM; col++)
    idwt_neo ((matrixInRowDetail + col * matrixInRowM),
		   (matrixInDetail + col * matrixInRowM),
		   matrixInRowM, lowReC, hiReC, filterLen,
		   (matrixOutDetailPre + col * matrixOutRow),
		   matrixOutRow);
  matrixOutDetailTemp =
    malloc (matrixOutRow * matrixInColM * sizeof (double));
  matrix_tran (matrixOutDetailPre, matrixInColM, matrixOutRow,
	       matrixOutDetailTemp, matrixInColM, matrixOutRow);
  free (matrixOutDetailPre);

  /* Final Inverse Transform */
  filterOutLength = 2 * matrixInColM + filterLen - 1;
  matrixOutPre = malloc (matrixOutRow * matrixOutCol * sizeof (double));
  for (row = 0; row < matrixOutRow; row++)
    idwt_neo ((matrixOutApproxTemp + row * matrixInColM),
		   (matrixOutDetailTemp + row * matrixInColM),
		   matrixInColM, lowReR, hiReR, filterLen,
		   (matrixOutPre + row * matrixOutCol),
		   matrixOutCol);
  free (matrixOutApproxTemp);
  free (matrixOutDetailTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixOutCol, matrixOut,
	       matrixOutRow, matrixOutCol);
  free (matrixOutPre);

  return;
}
