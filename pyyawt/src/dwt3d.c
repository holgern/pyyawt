/*
 * -------------------------------------------------------------------------
 * dwt3d.c -- 3-D signal decomposition and reconstruction
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

void dwt3d_tran(double *mat3DIn, int row1, int col1, int sli1,
		double *mat3DOut, int row2, int col2, int sli2)
{
  int i;

  //printf("enter dwt3d_tran!\n");
  for(i=0;i<sli1;i++)
    {
      //printf("%d\n",i);
      matrix_tran(mat3DIn+i*row1*col1, row1, col1,
		  mat3DOut+i*row1*col1, col1, row1);
    }

  return;
}

void dwt3d_line_forward(double *mat3DIn, int row1, int col1, int sli1,
			double *mat3DOutApp, double *mat3DOutDet,
			int row2, int col2, int sli2,
			double *loDe, double *hiDe, int filterLen,
			extend_method extMethod)
{
  int i;

  for(i=0;i<row1*sli1;i++)
    dwt_neo (mat3DIn+i*col1, col1, loDe,
	     hiDe, filterLen, mat3DOutApp+i*col2,
	     mat3DOutDet+i*col2, col2, extMethod);

  return;
}

void dwt3d_line_reverse(double *mat3DInApp, double *mat3DInDet,
			int row1, int col1, int sli1,
			double *mat3DOut, int row2, int col2,
			int sli2,
			double *loDe, double *hiDe, int filterLen)

{
  int i;

  for(i=0;i<row1*sli1;i++)
    idwt_neo (mat3DInApp+i*col1, mat3DInDet+i*col1, col1,
	      loDe, hiDe, filterLen,
	      mat3DOut+i*col2, col2);
  return;
}


void dwt3d_tran_z(double *mat3DIn, int row1, int col1, int sli1,
		  double *mat3DOut, int row2, int col2, int sli2)
{
  int i,j,k;

  for(i=0;i<row1;i++)
    {
      for(j=0;j<col1;j++)
	{
	  for(k=0;k<sli1;k++)
	    mat3DOut[k+j*sli1+i*col1*sli1] =
	      mat3DIn[k*row1*col1+i*col1+j];
	}
    }

  return;
}

void dwt3d_tran_z_inv(double *mat3DIn, int row1, int col1, int sli1,
		      double *mat3DOut, int row2, int col2, int sli2)
{
  int i,j,k;

  for(i=0;i<row2;i++)
    {
      for(j=0;j<col2;j++)
	{
	  for(k=0;k<sli2;k++)
	    mat3DOut[k*row2*col2+i*col2+j] =
	      mat3DIn[k+j*sli2+i*col2*sli2];
	}
    }

  return;
}

void dwt3d_combine(double *mat1, double *mat2, double *mat3,
		   double *mat4, double *mat5, double *mat6,
		   double *mat7, double *mat8, int rowIn,
		   int colIn, int sliIn, double *matOut,
		   int rowOut, int colOut, int sliOut)
{
  int siz;

  siz = rowOut*colOut*sliOut;

  dwt3d_tran_z_inv(mat1, rowIn, colIn, sliIn, matOut,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat2, rowIn, colIn, sliIn, matOut+siz,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat3, rowIn, colIn, sliIn, matOut+siz*2,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat4, rowIn, colIn, sliIn, matOut+siz*3,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat5, rowIn, colIn, sliIn, matOut+siz*4,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat6, rowIn, colIn, sliIn, matOut+siz*5,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat7, rowIn, colIn, sliIn, matOut+siz*6,
		   rowOut, colOut, sliOut);
  dwt3d_tran_z_inv(mat8, rowIn, colIn, sliIn, matOut+siz*7,
		   rowOut, colOut, sliOut);

  return;
}

void dwt3d_split(double *matIn, int rowIn, int colIn, int sliIn,
		 double *mat1, double *mat2, double *mat3,
		 double *mat4, double *mat5, double *mat6,
		 double *mat7, double *mat8, int rowOut,
		 int colOut, int sliOut)
{
  int siz;

  siz = rowIn*colIn*sliIn;
  dwt3d_tran_z(matIn, rowIn, colIn, sliIn,
	       mat1, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz, rowIn, colIn, sliIn,
	       mat2, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*2, rowIn, colIn, sliIn,
	       mat3, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*3, rowIn, colIn, sliIn,
	       mat4, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*4, rowIn, colIn, sliIn,
	       mat5, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*5, rowIn, colIn, sliIn,
	       mat6, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*6, rowIn, colIn, sliIn,
	       mat7, rowIn, colIn, sliIn);
  dwt3d_tran_z(matIn+siz*7, rowIn, colIn, sliIn,
	       mat8, rowIn, colIn, sliIn);
  return;
}


void dwt3(double *mat3DIn, int row, int col, int sli,
	  double *mat3DOut, int row2, int col2, int sli2,
	  int r, int c, int s, double *Lo1, double *Hi1,
	  double *Lo2, double *Hi2, double *Lo3, double *Hi3,
	  int fLen1, int fLen2, int fLen3, extend_method extMethod)
{
  double *tm1;
  double *rowApp, *rowDet, *rowAppRot, *rowDetRot;
  double *rowAppColApp, *rowAppColDet, *rowDetColApp, *rowDetColDet;
  double *rowAppColAppRot, *rowAppColDetRot;
  double *rowDetColAppRot, *rowDetColDetRot;
  double *rowAppColAppSliApp, *rowAppColAppSliDet;
  double *rowAppColDetSliApp, *rowAppColDetSliDet;
  double *rowDetColAppSliApp, *rowDetColAppSliDet;
  double *rowDetColDetSliApp, *rowDetColDetSliDet;

  /* rotate input cube for transform on rows */

  tm1 = malloc(row*col*sli*sizeof(double));
  dwt3d_tran(mat3DIn,col,row,sli,tm1,row,col,sli);

  /* transform on rows and obtain 2 cubes */

  rowApp = malloc(row*c*sli*sizeof(double));
  rowDet = malloc(row*c*sli*sizeof(double));
  dwt3d_line_forward(tm1, row, col, sli,
		     rowApp, rowDet, row, c , sli,
		     Lo1, Hi1, fLen1, extMethod);
  free(tm1);

  /* rotate 2 cubes for transform on columns */

  rowAppRot = malloc(row*c*sli*sizeof(double));
  dwt3d_tran(rowApp,row,c,sli,rowAppRot,c,row,sli);
  free(rowApp);
  rowDetRot = malloc(row*c*sli*sizeof(double));
  dwt3d_tran(rowDet,row,c,sli,rowDetRot,c,row,sli);
  free(rowDet);

  /* transform on column and obtain 4 cubes */

  rowAppColApp = malloc(r*c*sli*sizeof(double));
  rowAppColDet = malloc(r*c*sli*sizeof(double));
  dwt3d_line_forward(rowAppRot, c, row, sli,
		     rowAppColApp, rowAppColDet, c, r, sli,
		     Lo2, Hi2, fLen2, extMethod);
  free(rowAppRot);

  rowDetColApp = malloc(r*c*sli*sizeof(double));
  rowDetColDet = malloc(r*c*sli*sizeof(double));
  dwt3d_line_forward(rowDetRot, c, row, sli,
		     rowDetColApp, rowDetColDet, c, r, sli,
		     Lo2, Hi2, fLen2, extMethod);
  free(rowDetRot);

  /* rotate 4 cubes for transform on slices */

  rowAppColAppRot = malloc(r*c*sli*sizeof(double));
  dwt3d_tran_z(rowAppColApp,c,r,sli,rowAppColAppRot,r,sli,c);
  free(rowAppColApp);

  rowAppColDetRot = malloc(r*c*sli*sizeof(double));
  dwt3d_tran_z(rowAppColDet,c,r,sli,rowAppColDetRot,r,sli,c);
  free(rowAppColDet);

  rowDetColAppRot = malloc(r*c*sli*sizeof(double));
  dwt3d_tran_z(rowDetColApp,c,r,sli,rowDetColAppRot,r,sli,c);
  free(rowDetColApp);

  rowDetColDetRot = malloc(r*c*sli*sizeof(double));
  dwt3d_tran_z(rowDetColDet,c,r,sli,rowDetColDetRot,r,sli,c);
  free(rowDetColDet);

  /* transform on on slices and obtain eight cubes*/

  rowAppColAppSliApp = malloc(r*c*s*sizeof(double));
  rowAppColAppSliDet = malloc(r*c*s*sizeof(double));
  dwt3d_line_forward(rowAppColAppRot, r, sli, c,
		     rowAppColAppSliApp, rowAppColAppSliDet, r, s, c,
		     Lo3, Hi3, fLen3, extMethod);
  free(rowAppColAppRot);


  rowAppColDetSliApp = malloc(r*c*s*sizeof(double));
  rowAppColDetSliDet = malloc(r*c*s*sizeof(double));
  dwt3d_line_forward(rowAppColDetRot, r, sli, c,
		     rowAppColDetSliApp, rowAppColDetSliDet, r, s, c,
		     Lo3, Hi3, fLen3, extMethod);
  free(rowAppColDetRot);

  rowDetColAppSliApp = malloc(r*c*s*sizeof(double));
  rowDetColAppSliDet = malloc(r*c*s*sizeof(double));
  dwt3d_line_forward(rowDetColAppRot, r, sli, c,
		     rowDetColAppSliApp, rowDetColAppSliDet, r, s, c,
		     Lo3, Hi3, fLen3, extMethod);
  free(rowDetColAppRot);


  rowDetColDetSliApp = malloc(r*c*s*sizeof(double));
  rowDetColDetSliDet = malloc(r*c*s*sizeof(double));
  dwt3d_line_forward(rowDetColDetRot, r, sli, c,
		     rowDetColDetSliApp, rowDetColDetSliDet, r, s, c,
		     Lo3, Hi3, fLen3, extMethod);
  free(rowDetColDetRot);

  /* combine all the cubes to the one */
  dwt3d_combine(rowAppColAppSliApp, rowAppColAppSliDet,
		rowAppColDetSliApp, rowAppColDetSliDet,
		rowDetColAppSliApp, rowDetColAppSliDet,
		rowDetColDetSliApp, rowDetColDetSliDet,
		r, s, c, mat3DOut, c, r, s);
  free(rowAppColAppSliApp);
  free(rowAppColAppSliDet);
  free(rowAppColDetSliApp);
  free(rowAppColDetSliDet);
  free(rowDetColAppSliApp);
  free(rowDetColAppSliDet);
  free(rowDetColDetSliApp);
  free(rowDetColDetSliDet);

  return;
}

void idwt3(double *mat3DIn, int row, int col, int sli,
	   double *mat3DOut, int r, int c, int s,
	   double *Lo1, double *Hi1, double *Lo2,
	   double *Hi2, double *Lo3, double *Hi3,
	   int fLen1, int fLen2, int fLen3)
{
  double *rowApp, *rowDet, *rowAppRot, *rowDetRot;
  double *rowAppColApp, *rowAppColDet, *rowDetColApp, *rowDetColDet;
  double *rowAppColAppRot, *rowAppColDetRot;
  double *rowDetColAppRot, *rowDetColDetRot;
  double *rowAppColAppSliApp, *rowAppColAppSliDet;
  double *rowAppColDetSliApp, *rowAppColDetSliDet;
  double *rowDetColAppSliApp, *rowDetColAppSliDet;
  double *rowDetColDetSliApp, *rowDetColDetSliDet;
  double *outTemp;

  /* split the result to 8 cubes, rotation also done */
  rowAppColAppSliApp = malloc(row*col*sli*sizeof(double));
  rowAppColAppSliDet = malloc(row*col*sli*sizeof(double));
  rowAppColDetSliApp = malloc(row*col*sli*sizeof(double));
  rowAppColDetSliDet = malloc(row*col*sli*sizeof(double));
  rowDetColAppSliApp = malloc(row*col*sli*sizeof(double));
  rowDetColAppSliDet = malloc(row*col*sli*sizeof(double));
  rowDetColDetSliApp = malloc(row*col*sli*sizeof(double));
  rowDetColDetSliDet = malloc(row*col*sli*sizeof(double));
  dwt3d_split(mat3DIn, col, row, sli,
	      rowAppColAppSliApp,
	      rowAppColAppSliDet,
	      rowAppColDetSliApp,
	      rowAppColDetSliDet,
	      rowDetColAppSliApp,
	      rowDetColAppSliDet,
	      rowDetColDetSliApp,
	      rowDetColDetSliDet,
	      row, sli, col);

  /* merge 8 cubes to 4 */
  rowAppColAppRot = malloc(row*s*col*sizeof(double));
  dwt3d_line_reverse(rowAppColAppSliApp, rowAppColAppSliDet,
		     row, sli, col, rowAppColAppRot,
		     row, s, col, Lo1, Hi1, fLen1);
  free(rowAppColAppSliApp);
  free(rowAppColAppSliDet);

  rowAppColDetRot = malloc(row*s*col*sizeof(double));
  dwt3d_line_reverse(rowAppColDetSliApp, rowAppColDetSliDet,
		     row, sli, col, rowAppColDetRot,
		     row, s, col, Lo1, Hi1, fLen1);
  free(rowAppColDetSliApp);
  free(rowAppColDetSliDet);

  rowDetColAppRot = malloc(row*s*col*sizeof(double));
  dwt3d_line_reverse(rowDetColAppSliApp, rowDetColAppSliDet,
		     row, sli, col, rowDetColAppRot,
		     row, s, col, Lo1, Hi1, fLen1);
  free(rowDetColAppSliApp);
  free(rowDetColAppSliDet);

  rowDetColDetRot = malloc(row*s*col*sizeof(double));
  dwt3d_line_reverse(rowDetColDetSliApp, rowDetColDetSliDet,
		     row, sli, col, rowDetColDetRot,
		     row, s, col, Lo1, Hi1, fLen1);
  free(rowDetColDetSliApp);
  free(rowDetColDetSliDet);

  /* rotate 4 cubes along z axix*/
  rowAppColApp = malloc(row*s*col*sizeof(double));
  dwt3d_tran_z_inv(rowAppColAppRot, row, s, col,
		   rowAppColApp, col, row, s);
  free(rowAppColAppRot);

  rowAppColDet = malloc(row*s*col*sizeof(double));
  dwt3d_tran_z_inv(rowAppColDetRot, row, s, col,
		   rowAppColDet, col, row, s);
  free(rowAppColDetRot);

  rowDetColApp = malloc(row*s*col*sizeof(double));
  dwt3d_tran_z_inv(rowDetColAppRot, row, s, col,
		   rowDetColApp, col, row, s);
  free(rowDetColAppRot);

  rowDetColDet = malloc(row*s*col*sizeof(double));
  dwt3d_tran_z_inv(rowDetColDetRot, row, s, col,
		   rowDetColDet, col, row, s);
  free(rowDetColDetRot);

  /* merge 4 cubes to 2 */
  rowAppRot = malloc(r*s*col*sizeof(double));
  dwt3d_line_reverse(rowAppColApp, rowAppColDet,
		     col, row, s, rowAppRot,
		     col, r, s, Lo2, Hi2, fLen2);
  free(rowAppColApp);
  free(rowAppColDet);

  rowDetRot = malloc(r*s*col*sizeof(double));
  dwt3d_line_reverse(rowDetColApp, rowDetColDet,
		     col, row, s, rowDetRot,
		     col, r, s, Lo2, Hi2, fLen2);
  free(rowDetColApp);
  free(rowDetColDet);

  /* rotate 2 cubes along xy */
  rowApp = malloc(r*s*col*sizeof(double));
  dwt3d_tran(rowAppRot, col, r, s, rowApp, r, col, s);
  free(rowAppRot);

  rowDet = malloc(r*s*col*sizeof(double));
  dwt3d_tran(rowDetRot, col, r, s, rowDet, r, col, s);
  free(rowDetRot);

  /* merge 2 cubes to one */
  outTemp = malloc(r*s*c*sizeof(double));
  dwt3d_line_reverse(rowApp, rowDet,
		     r, col, s, outTemp,
		     r, c, s, Lo3, Hi3, fLen3);
  free(rowApp);
  free(rowDet);

  /* rotate the cube to fit scliab matrix format */
  dwt3d_tran(outTemp, r, c, s, mat3DOut, c, r, s);
  free(outTemp);
  return;
}
