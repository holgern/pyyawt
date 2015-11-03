/*
 * -------------------------------------------------------------------------
 * utility.c -- utility function
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
 * Matrix Transposition Operation
 *-----------------------------------------*/

void
matrix_tran (double *matrixIn, int matrixInRow, int matrixInCol,
	     double *matrixOut, int matrixOutRow, int matrixOutCol)
{
  int row, col;

  for (col = 0; col < matrixInCol; col++)
    {
      for (row = 0; row < matrixInRow; row++)
	{
	  matrixOut[row + col * matrixInRow] =
	    matrixIn[col + row * matrixInCol];
	}
    }
  return;
}

/*-------------------------------------------
 * Flipping Operation
 *-----------------------------------------*/

void
wrev (const double *sigIn, int sigInLength,
      double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength; count++)
    sigOut[count] = sigIn[sigInLength - count - 1];
  return;
}

/*-------------------------------------------
 * Quadrature Mirror Filtering Operation
 *-----------------------------------------*/

void
qmf_odd (double *sigIn, int sigInLength,
     double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength; count++)
    {
      sigOut[count] = sigIn[sigInLength - count - 1];
      if (count % 2 == 0)
	  {
	    sigOut[count] = -1 * sigOut[count];
	  }
    }
    return;
}

void
qmf_even (const double *sigIn, int sigInLength,
     double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength; count++)
    {
      sigOut[count] = sigIn[sigInLength - count - 1];
      if (count % 2 != 0)
	  {
	    sigOut[count] = -1 * sigOut[count];
	  }
    }
  return;
}

/*-------------------------------------------
 * Flipping and QMF at the same time
 *-----------------------------------------*/
void
qmf_wrev (const double *sigIn, int sigInLength,
	  double *sigOut, int sigOutLength)
{
  int count = 0;
  double *sigOutTemp;
  sigOutTemp = malloc(sigInLength*sizeof(double));
  for (count = 0; count < sigInLength; count++)
    {
      sigOutTemp[count] = sigIn[sigInLength - count - 1];
      if (count % 2 != 0)
	{
	  sigOutTemp[count] = -1 * sigOutTemp[count];
	}
    }

  for (count = 0; count < sigInLength; count++)
    sigOut[count] = sigOutTemp[sigInLength - count - 1];
  free(sigOutTemp);
  return;
}

/*-------------------------------------------
 * Verbatim Copying
 *-----------------------------------------*/

void
verbatim_copy (const double *sigIn, int sigInLength,
	       double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength; count++)
    sigOut[count] = sigIn[count];
}

/*-------------------------------------------
 * Dyaddown
 *-----------------------------------------*/

void
dyaddown_1D_keep_odd (double *sigIn, int sigInLength,
		   double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = sigIn[count * 2];
  return;
}

void
dyaddown_1D_keep_even (double *sigIn, int sigInLength,
		    double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = sigIn[count * 2 + 1];
  return;
}

void
dyaddown_2D_keep_odd_row (double *matrixIn, int matrixInRow,
			  int matrixInCol, double *matrixOut,
			  int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixOutRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	matrixOutPre[col+row*matrixInCol]=
	  matrixInTemp[col+row*matrixInCol*2];
    }
  free (matrixInTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOut, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  return;
}

void
dyaddown_2D_keep_odd_col (double *matrixIn, int matrixInRow,
			  int matrixInCol, double *matrixOut,
			  int matrixOutRow, int matrixOutCol)
{
  int row, col;
  for(col=0;col<matrixOutCol;col++)
    {
      for(row=0;row<matrixInRow;row++)
	matrixOut[row+col*matrixInRow]=
	  matrixIn[row+col*matrixInRow*2];
    }
  return;
}

void
dyaddown_2D_keep_even_row (double *matrixIn, int matrixInRow,
			  int matrixInCol, double *matrixOut,
			  int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixOutRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	matrixOutPre[col+row*matrixInCol]=
	  matrixInTemp[col+(2*row+1)*matrixInCol];
    }
  free (matrixInTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOut, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  return;
}

void
dyaddown_2D_keep_even_col (double *matrixIn, int matrixInRow,
			  int matrixInCol, double *matrixOut,
			  int matrixOutRow, int matrixOutCol)
{
  int row, col;
  for(col=0;col<matrixOutCol;col++)
    {
      for(row=0;row<matrixInRow;row++)
	matrixOut[row+col*matrixInRow]=
	  matrixIn[row+(2*col+1)*matrixInRow];
    }
  return;
}


void
dyaddown_2D_keep_odd (double *matrixIn, int matrixInRow,
		      int matrixInCol, double *matrixOut,
		      int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixOutRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	matrixOutPre[col+row*matrixInCol]=
	  matrixInTemp[col+2*row*matrixInCol];
    }
  free (matrixInTemp);

  matrixOutTemp = malloc (matrixOutRow * matrixInCol * sizeof (double));
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOutTemp, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  for(col=0;col<matrixOutCol;col++)
    {
      for(row=0;row<matrixOutRow;row++)
	matrixOut[row+col*matrixOutRow]=
	  matrixOutTemp[row+2*col*matrixOutRow];
    }
  free(matrixOutTemp);
  return;
}

void
dyaddown_2D_keep_even (double *matrixIn, int matrixInRow,
		       int matrixInCol, double *matrixOut,
		       int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixOutRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	matrixOutPre[col+row*matrixInCol]=
	  matrixInTemp[col+(2*row+1)*matrixInCol];
    }
  free (matrixInTemp);

  matrixOutTemp = malloc (matrixOutRow * matrixInCol * sizeof (double));
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOutTemp, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  for(col=0;col<matrixOutCol;col++)
    {
      for(row=0;row<matrixOutRow;row++)
	matrixOut[row+col*matrixOutRow]=
	  matrixOutTemp[row+(2*col+1)*matrixOutRow];
    }
  free(matrixOutTemp);
  return;
}

/*-------------------------------------------
 * Dyadup
 *-----------------------------------------*/

void
dyadup_1D_feed_odd (double *sigIn, int sigInLength,
		 double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength - 1; count++)
    {
      sigOut[count * 2 + 1] = 0;
      sigOut[count * 2] = sigIn[count];
    }
  sigOut[sigOutLength-1] = sigIn[sigInLength-1];
  return;
}

void
dyadup_1D_feed_even (double *sigIn, int sigInLength,
		  double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigInLength; count++)
    {
      //printf("count %d\n",count);
      sigOut[count * 2] = 0;
      sigOut[count * 2 + 1] = sigIn[count];
    }
  sigOut[sigOutLength-1] = 0;
  return;
}

void
dyadup_2D_feed_odd_row (double *matrixIn, int matrixInRow,
			int matrixInCol, double *matrixOut,
			int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixInRow-1;row++)
    {
      for(col=0;col<matrixInCol;col++)
	{
	  matrixOutPre[col+2*row*matrixInCol]=
	    matrixInTemp[col+row*matrixInCol];
	  matrixOutPre[col+(2*row+1)*matrixInCol] = 0;
	}
    }
  for(col=0;col<matrixInCol;col++)
    matrixOutPre[col+(matrixOutRow-1)*matrixInCol] =
      matrixInTemp[col+(matrixInRow-1)*matrixInCol];
  free (matrixInTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOut, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  return;
}

void
dyadup_2D_feed_odd_col (double *matrixIn, int matrixInRow,
			int matrixInCol, double *matrixOut,
			int matrixOutRow, int matrixOutCol)
{
  int row, col;
  for(col=0;col<matrixInCol-1;col++)
    {
      for(row=0;row<matrixInRow;row++)
	{
	  matrixOut[row+2*col*matrixInRow]=
	    matrixIn[row+col*matrixInRow];
	  matrixOut[row+(2*col+1)*matrixInRow] = 0;
	}
    }
  for(row=0;row<matrixInRow;row++)
    matrixOut[row+(matrixOutCol-1)*matrixInRow] =
      matrixIn[row+(matrixInCol-1)*matrixInRow];
  return;
}

void
dyadup_2D_feed_even_row (double *matrixIn, int matrixInRow,
			int matrixInCol, double *matrixOut,
			int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixInRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	{
	  matrixOutPre[col+2*row*matrixInCol]= 0;
	  matrixOutPre[col+(2*row+1)*matrixInCol]
	    = matrixInTemp[col+row*matrixInCol];
	}
    }
  for(col=0;col<matrixInCol;col++)
    matrixOutPre[col+(matrixOutRow-1)*matrixInCol] = 0;
  free (matrixInTemp);
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOut, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  return;
}

void
dyadup_2D_feed_even_col (double *matrixIn, int matrixInRow,
			int matrixInCol, double *matrixOut,
			int matrixOutRow, int matrixOutCol)
{
  int row, col;
  for(col=0;col<matrixInCol;col++)
    {
      for(row=0;row<matrixInRow;row++)
	{
	  matrixOut[row+2*col*matrixInRow] = 0;
	  matrixOut[row+(2*col+1)*matrixInRow] =
	    matrixIn[row+col*matrixInRow];
	}
    }
  for(row=0;row<matrixOutRow;row++)
    matrixOut[row+(matrixOutCol-1)*matrixOutRow]=0;
  return;
}


void
dyadup_2D_feed_odd (double *matrixIn, int matrixInRow,
		    int matrixInCol, double *matrixOut,
		    int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixInRow-1;row++)
    {
      for(col=0;col<matrixInCol;col++)
	{
	  matrixOutPre[col+2*row*matrixInCol] =
	    matrixInTemp[col+row*matrixInCol];
	  matrixOutPre[col+(2*row+1)*matrixInCol] = 0;
	}
    }
  for(col=0;col<matrixInCol;col++)
    matrixOutPre[col+(matrixOutRow-1)*matrixInCol] =
      matrixInTemp[col+(matrixInRow-1)*matrixInCol];
  free (matrixInTemp);

  matrixOutTemp = malloc (matrixOutRow * matrixInCol * sizeof (double));
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOutTemp, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  for(col=0;col<matrixInCol-1;col++)
    {
      for(row=0;row<matrixOutRow;row++)
	{
	  matrixOut[row+2*col*matrixOutRow]=
	    matrixOutTemp[row+col*matrixOutRow];
	  matrixOut[row+(2*col+1)*matrixOutRow] = 0;
	}
    }
  for(row=0;row<matrixOutRow;row++)
    matrixOut[row+(matrixOutCol-1)*matrixOutRow] =
      matrixOutTemp[row+(matrixInCol-1)*matrixOutRow];
  free(matrixOutTemp);
  return;
}

void
dyadup_2D_feed_even (double *matrixIn, int matrixInRow,
		     int matrixInCol, double *matrixOut,
		     int matrixOutRow, int matrixOutCol)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixOutRow * matrixInCol * sizeof (double));
  for(row=0;row<matrixInRow;row++)
    {
      for(col=0;col<matrixInCol;col++)
	{
	  matrixOutPre[col+(2*row+1)*matrixInCol]=
	    matrixInTemp[col+row*matrixInCol];
	  matrixOutPre[col+2*row*matrixInCol] = 0;
	}
    }
  for(col=0;col<matrixInCol;col++)
    matrixOutPre[col+(matrixOutRow-1)*matrixInCol] = 0;
  free (matrixInTemp);

  matrixOutTemp = malloc (matrixOutRow * matrixInCol * sizeof (double));
  matrix_tran (matrixOutPre, matrixOutRow, matrixInCol,
	       matrixOutTemp, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  for(col=0;col<matrixInCol;col++)
    {
      for(row=0;row<matrixOutRow;row++)
	{
	  matrixOut[row+(2*col+1)*matrixOutRow]=
	    matrixOutTemp[row+col*matrixOutRow];
	  matrixOut[row+2*col*matrixOutRow] = 0;
	}
    }
  for(row=0;row<matrixOutRow;row++)
    matrixOut[row+(matrixOutCol-1)*matrixOutRow]=0;
  free(matrixOutTemp);
  return;
}


/*-------------------------------------------
 * Signal Extending
 *-----------------------------------------*/

void
extend_method_parse (char *mode, extend_method *extMethod)
{
  int count;

  for (count=0;count<extensionIdentityNum;count++)
    {
      if (strcmp(mode,ei[count].extMethodName) == 0)
	{
	  *extMethod = ei[count].extMethod;
	  break;
	}
    }
  return;
}


void
wextend_1D_center (double *sigIn, int sigInLength,
	 double *sigOut, int sigOutLength,
	 extend_method extMethod)
{
  int count = 0;
  int addLength = 0;

  addLength = (sigOutLength - sigInLength) >> 1;
  for (count = 0; count < addLength; count++)
    {
      sigOut[count] = 0;
      sigOut[count + sigInLength + addLength] = 0;
    }

  for (count = 0; count < sigInLength; count++)
    {
      sigOut[count + addLength] = sigIn[count];
    }

  switch (extMethod) {
  case ZPD: break;
  case SYMH:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[addLength - count - 1];
	  sigOut[count + sigInLength + addLength] =
	    sigIn[sigInLength - count - 1];
	}
      break;
    }
  case SYMW:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[addLength - count];
	  sigOut[count + sigInLength + addLength] =
	    sigIn[sigInLength - count - 2];
	}
      break;
    }
  case ASYMH:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[addLength - count - 1] * (-1);
	  sigOut[count + sigInLength + addLength] =
	    sigIn[sigInLength - count - 1] * (-1);
	}
      break;
    }
  case ASYMW:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[addLength - count] * (-1);
	  sigOut[count + sigInLength + addLength] =
	    sigIn[sigInLength - count - 2] * (-1);
	}
      break;
    }
  case SP0:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[0];
	  sigOut[count + sigInLength + addLength] =
	    sigIn[sigInLength - 1];
	}
      break;
    }
  case SP1:
    {
      for (count = (addLength - 1); count >= 0; count--)
	{
		sigOut[count] = sigIn[0]-(sigIn[1]-sigIn[0])*(addLength-count);
		sigOut[sigInLength + 2 * addLength - count - 1] =
			sigIn[sigInLength - 1] - (sigIn[sigInLength-2] - sigIn[sigInLength-1])*(addLength-count);
	 }
      break;
    }
  case PPD:
    {
      for (count = 0; count < addLength; count++)
	{
	  sigOut[count] = sigIn[sigInLength - addLength + count];
	  sigOut[count + sigInLength + addLength] = sigIn[count];
	}
      break;
    }
  case PER:
    {
      if (sigInLength%2 == 0)
	{
	  for (count = 0; count < addLength; count++)
	    {
	      sigOut[count] = sigIn[sigInLength - addLength + count];
	      sigOut[count + sigInLength + addLength] = sigIn[count];
	    }
	}
      else
	{
	  sigOut[addLength + sigInLength] = sigIn[sigInLength - 1];
	  for (count = 0; count < addLength; count++)
	    {
	      sigOut[count] =
		     sigOut[sigInLength + 1 + count];
		  sigOut[count + sigInLength + addLength + 1] =
		     sigIn[count];
	    }
	}
      break;
    }
  default: break;
  }
  return;
}


void
wextend_1D_left (double *sigIn, int sigInLength,
	 double *sigOut, int sigOutLength,
	 extend_method extMethod)
{
  int count = 0;
  int addLength = 0;

  addLength = (sigOutLength - sigInLength);
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = 0;

  for (count = 0; count < sigInLength; count++)
    {
      sigOut[count + addLength] = sigIn[count];
    }

  switch (extMethod) {
  case ZPD: break;
  case SYMH:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[addLength - count - 1];
      break;
    }
  case SYMW:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[addLength - count];
      break;
    }
  case ASYMH:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[addLength - count - 1] * (-1);
      break;
    }
  case ASYMW:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[addLength - count] * (-1);
      break;
    }
  case SP0:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[0];
      break;
    }
  case SP1:
    {
      for (count = addLength - 1; count >= 0; count--)
         sigOut[count] = sigIn[0]-(sigIn[1]-sigIn[0])*(addLength-count);


      break;
    }
  case PPD:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count] = sigIn[sigInLength - addLength + count];
      break;
    }
  case PER:
    {
      if (sigInLength%2 == 0)
	{
	  for (count = 0; count < addLength; count++)
	    sigOut[count] = sigIn[sigInLength - addLength + count];
	}
      else
	{
	  for (count = 0; count < sigInLength; count++)
	    {
	      sigOut[count + addLength - 1] = sigIn[count];
	    }
	  sigOut[sigOutLength-1] = sigOut[sigOutLength-2];
	  for (count = 0; count < (addLength - 1); count++)
	    sigOut[count] =
	      sigOut[sigInLength + count + 1];
	}
      break;
    }
  default: break;
  }
  return;
}

void
wextend_1D_right (double *sigIn, int sigInLength,
	 double *sigOut, int sigOutLength,
	 extend_method extMethod)
{
  int count = 0;
  int addLength = 0;

  addLength = (sigOutLength - sigInLength);
  for (count = 0; count < addLength; count++)
    sigOut[count + sigInLength] = 0;

  for (count = 0; count < sigInLength; count++)
    sigOut[count] = sigIn[count];

  switch (extMethod) {
  case ZPD: break;
  case SYMH:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count + sigInLength] = sigIn[sigInLength - count - 1];
      break;
    }
  case SYMW:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count + sigInLength] = sigIn[sigInLength - count - 2];
      break;
    }
  case ASYMH:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count + sigInLength] =
	  sigIn[sigInLength - count - 1] * (-1);
      break;
    }
  case ASYMW:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count + sigInLength] =
	  sigIn[sigInLength - count - 2] * (-1);
      break;
    }
  case SP0:
    {
      for (count = 0; count < addLength; count++)
	  sigOut[count + sigInLength] = sigIn[sigInLength - 1];
      break;
    }
  case SP1:
    {
      for (count = 0; count < addLength; count++)
	  	sigOut[count+sigInLength] = sigIn[sigInLength-1] - (sigIn[sigInLength-2]-sigIn[sigInLength-1])*(count+1);

      break;
    }
  case PPD:
    {
      for (count = 0; count < addLength; count++)
	sigOut[count + sigInLength] = sigIn[count];
      break;
    }
  case PER:
    {
      if (sigInLength%2 == 0)
	{
	  for (count = 0; count < addLength; count++)
	    sigOut[count + sigInLength] = sigIn[count];
	}
      else
	{
	  sigOut[sigInLength] = sigOut[sigInLength - 1];
	  for(count = 0; count < (addLength - 1); count++)
	    {
	      sigOut[sigInLength + count + 1] = sigOut[count];
	    }
	}
      break;
    }
  default: break;
  }
  return;
}

void
wextend_2D (double *matrixIn, int matrixInRow, int matrixInCol,
	    double *matrixOut, int matrixOutRow, int matrixOutCol,
	    extend_method extMethod, char *rowOpt, char *colOpt)
{
  int row, col;
  double *matrixInTemp, *matrixOutPre, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixInRow * matrixOutCol * sizeof (double));

  for (row = 0; row < matrixInRow; row++)
    {
      if (*rowOpt=='b')
	wextend_1D_center((matrixInTemp + row * matrixInCol),
			  matrixInCol,
			  (matrixOutPre + row * matrixOutCol),
			  matrixOutCol, extMethod);
      if (*rowOpt=='l')
	wextend_1D_left((matrixInTemp + row * matrixInCol),
			matrixInCol,
			(matrixOutPre + row * matrixOutCol),
			matrixOutCol, extMethod);
      if (*rowOpt=='r')
	wextend_1D_right((matrixInTemp + row * matrixInCol),
			 matrixInCol,
			 (matrixOutPre + row * matrixOutCol),
			 matrixOutCol, extMethod);
    }
  free (matrixInTemp);

  matrixOutTemp = malloc (matrixInRow * matrixOutCol * sizeof (double));
  matrix_tran (matrixOutPre, matrixInRow, matrixOutCol,
	       matrixOutTemp, matrixInRow, matrixOutCol);
  free (matrixOutPre);

  for (col = 0; col < matrixOutCol; col++)
    {
      if (*colOpt=='b')
	wextend_1D_center((matrixOutTemp + col * matrixInRow),
			  matrixInRow,
			  (matrixOut + col * matrixOutRow),
			  matrixOutRow, extMethod);
      if (*colOpt=='l')
	wextend_1D_left((matrixOutTemp + col * matrixInRow),
			matrixInRow,
			(matrixOut + col * matrixOutRow),
			matrixOutRow, extMethod);
      if (*colOpt=='r')
	wextend_1D_right((matrixOutTemp + col * matrixInRow),
			 matrixInRow,
			 (matrixOut + col * matrixOutRow),
			 matrixOutRow, extMethod);
    }
  free (matrixOutTemp);
  return;
}

void
wextend_2D_col (double *matrixIn, int matrixInRow, int matrixInCol,
		double *matrixOut, int matrixOutRow, int matrixOutCol,
		extend_method extMethod, char *Opt)
{
  int row;//, col;
  double *matrixInTemp, *matrixOutPre;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixOutCol, matrixInRow);

  matrixOutPre = malloc (matrixInRow * matrixOutCol * sizeof (double));
  for (row = 0; row < matrixInRow; row++)
    {
      if (strcmp(Opt,"b") == 0)
	wextend_1D_center((matrixInTemp + row * matrixInCol),
			  matrixInCol,
			  (matrixOutPre + row * matrixOutCol),
			  matrixOutCol, extMethod);
      if (strcmp(Opt,"l") == 0)
	wextend_1D_left((matrixInTemp + row * matrixInCol),
			matrixInCol,
			(matrixOutPre + row * matrixOutCol),
			matrixOutCol, extMethod);
      if (strcmp(Opt,"r") == 0)
	wextend_1D_right((matrixInTemp + row * matrixInCol),
			 matrixInCol,
			 (matrixOutPre + row * matrixOutCol),
			 matrixOutCol, extMethod);
    }
  free (matrixInTemp);
  matrix_tran (matrixOutPre, matrixInRow, matrixOutCol,
	       matrixOut, matrixInRow, matrixOutCol);
  free (matrixOutPre);
  return;
}

void
wextend_2D_row (double *matrixIn, int matrixInRow, int matrixInCol,
		double *matrixOut, int matrixOutRow, int matrixOutCol,
		extend_method extMethod, char *Opt)
{
  int col;
  //int row, col;
  for (col = 0; col < matrixInCol; col++)
    {
      if (strcmp(Opt,"b") == 0)
	wextend_1D_center((matrixIn + col * matrixInRow),
			  matrixInRow,
			  (matrixOut + col * matrixOutRow),
			  matrixOutRow, extMethod);
      if (strcmp(Opt,"l") == 0)
	wextend_1D_left((matrixIn + col * matrixInRow),
			matrixInRow,
			(matrixOut + col * matrixOutRow),
			matrixOutRow, extMethod);
      if (strcmp(Opt,"r") == 0)
	wextend_1D_right((matrixIn + col * matrixInRow),
			 matrixInRow,
			 (matrixOut + col * matrixOutRow),
			 matrixOutRow, extMethod);
    }
  return;
}


/*-------------------------------------------
 * Signal Extraction
 *-----------------------------------------*/
void
wkeep_1D_center (double *sigIn, int sigInLength,
		 double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
	sigOut[count] = sigIn[(sigInLength -
			       sigOutLength) / 2 + count];
  return;
}

void
wkeep_1D_left (double *sigIn, int sigInLength,
	       double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = sigIn[count];
  return;
}

void
wkeep_1D_right (double *sigIn, int sigInLength,
		double *sigOut, int sigOutLength)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = sigIn[sigInLength - sigOutLength + count];
  return;
}

void
wkeep_1D_index (double *sigIn, int sigInLength,
		double *sigOut, int sigOutLength, int first)
{
  int count = 0;
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = sigIn[count + first - 1];
  return;
}

void
wkeep_2D_center (double *matrixIn, int matrixInRow,
		 int matrixInCol, double *matrixOut,
		 int matrixOutRow, int matrixOutCol)
{
  int row, col, startRow, startCol;
  double *matrixInTemp, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrixOutTemp = malloc (matrixOutRow * matrixOutCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixInCol, matrixInRow);

  startRow = (matrixInRow - matrixOutRow) / 2;
  startCol = (matrixInCol - matrixOutCol) / 2;

  for (row = startRow; row < startRow + matrixOutRow; row++)
    {
      for (col = startCol; col < startCol + matrixOutCol; col++)
	matrixOutTemp[col - startCol + (row - startRow) * matrixOutCol] =
	  matrixInTemp[col + row * matrixInCol];
    }
  matrix_tran (matrixOutTemp, matrixOutRow, matrixOutCol,
	       matrixOut, matrixOutRow, matrixOutCol);

  free (matrixInTemp);
  free (matrixOutTemp);
  return;
}

void
wkeep_2D_index (double *matrixIn, int matrixInRow,
		int matrixInCol, double *matrixOut,
		int matrixOutRow, int matrixOutCol,
		int rowFirst, int colFirst)
{
  int row, col, startRow, startCol;
  double *matrixInTemp, *matrixOutTemp;

  matrixInTemp = malloc (matrixInRow * matrixInCol * sizeof (double));
  matrixOutTemp = malloc (matrixOutRow * matrixOutCol * sizeof (double));
  matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixInTemp, matrixInCol, matrixInRow);

  startRow = rowFirst - 1;
  startCol = colFirst - 1;

  for (row = startRow; row < startRow + matrixOutRow; row++)
    {
      for (col = startCol; col < startCol + matrixOutCol; col++)
	matrixOutTemp[col - startCol + (row - startRow) * matrixOutCol] =
	  matrixInTemp[col + row * matrixInCol];
    }
  matrix_tran (matrixOutTemp, matrixOutRow, matrixOutCol,
	       matrixOut, matrixOutRow, matrixOutCol);
  free (matrixInTemp);
  free (matrixOutTemp);
  return;
}

/*-------------------------------------------
 * Convolution
 *-----------------------------------------*/
void
conv (double *sigIn, int sigInLength,
      double *sigOut, int sigOutLength,
      double *filter, int filterLength)
{
  int count = 0;
  int bufferLength = 0;
  int conCount = 0;
  double *pBuf;

  bufferLength = sigInLength + 2 * (filterLength - 1);
  pBuf = malloc (bufferLength * sizeof (double));

  /* initialize the buffer */
  for (count = 0; count < filterLength - 1; count++)
     {
      pBuf[count] = 0;		/* head */
      pBuf[count + sigInLength + filterLength - 1] = 0;	/* tail */
    }
  for (count = 0; count < sigInLength; count++)
    pBuf[count + filterLength - 1] = sigIn[count];

  /* convolution */
  for (count = 0; count < sigOutLength; count++)
    {
      sigOut[count] = 0;
      for (conCount = filterLength - 1; conCount >= 0; conCount--)
	   {
	     sigOut[count] += filter[conCount] *
	     pBuf[count + filterLength - conCount - 1];
	   }
    }
  free (pBuf);
  return;
}

/*-------------------------------------------
 * Periodic Convolution
 *-----------------------------------------*/

void
i_conv(double *sigIn, int sigInLength,
	  double *sigOut, int sigOutLength,
	  double *filter, int filterLength)
{
  int count = 0;
  int bufferLength = 0;
  int outLength = 0;
  double *pBuf;
  double *pOutBuf;

  bufferLength = 2 * sigInLength;
  pBuf = malloc(bufferLength*sizeof(double));
  for(count=0;count<sigInLength;count++)
  {
	  pBuf[count] = sigIn[count];
	  pBuf[count+sigInLength] = sigIn[count];
  }

  outLength = filterLength + 2 * sigInLength - 1;
  pOutBuf = malloc(outLength*sizeof(double));
  conv(pBuf,bufferLength,pOutBuf,outLength,filter,filterLength);

  free(pBuf);

  for(count=0;count<sigOutLength;count++)
    sigOut[count] = pOutBuf[filterLength+count];
  free(pOutBuf);
  return;
}

/*-------------------------------------------
 * EXP2 Convolution
 *-----------------------------------------*/

 void swt_exp2(int lev, int *outputV)
 {
	 int count;
	 *outputV = 1;
	 if (lev>0)
	 {
		 for (count=0;count<lev;count++)
			 *outputV = *outputV * 2;
	 }
	 return;
 }

 /*-------------------------------------------
 * abs
 *-----------------------------------------*/
 double swt_abs(double sigIn)
 {
     if (sigIn>=0)
		 return sigIn;
	 else
		 return (-1*sigIn);
 }

/*-------------------------------------------
 * min
 *-----------------------------------------*/

 void swt_min(double *sigIn, int sigInLength, double *sigMin)
 {
	 int count;

	 *sigMin = sigIn[0];
	 for(count=1;count<sigInLength;count++)
	 {
         if (sigIn[count] < *sigMin)
			 *sigMin = sigIn[count];
	 }
	 return;
 }

 void swt_min_abs(double *sigIn, int sigInLength, double *sigMin)
 {
	 int count;

	 *sigMin = swt_abs(sigIn[0]);
	 for(count=1;count<sigInLength;count++)
	 {
         if (swt_abs(sigIn[count]) < *sigMin)
			 *sigMin = swt_abs(sigIn[count]);
	 }
	 return;
 }


 /*-------------------------------------------
 * max
 *-----------------------------------------*/

 void swt_max(double *sigIn, int sigInLength, double *sigMax)
 {
	 int count;

	 *sigMax = sigIn[0];
	 for(count=1;count<sigInLength;count++)
	 {
         if (sigIn[count] > *sigMax)
			 *sigMax = sigIn[count];
	 }
	 return;
 }

 void swt_max_abs(double *sigIn, int sigInLength, double *sigMax)
 {
	 int count;

	 *sigMax = swt_abs(sigIn[0]);
	 for(count=1;count<sigInLength;count++)
	 {
         if (swt_abs(sigIn[count]) > *sigMax)
			 *sigMax = swt_abs(sigIn[count]);
	 }
	 return;
 }


/*-------------------------------------------
 * wcodemat
 *-----------------------------------------*/

 void wcodemat(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv)
 {
     int count;
     double sigMin, sigMax;

	 swt_max(sigIn, sigInLength, &sigMax);
	 swt_min(sigIn, sigInLength, &sigMin);

	 for(count=0;count<sigInLength;count++)
     	 sigOut[count] = ((sigIn[count]-sigMin)/(sigMax-sigMin))*(maxv-minv) + minv;
	 return;
 }

 void wcodematd(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, double minv, double maxv)
 {
     int count;
     double sigMin, sigMax;

	 swt_max(sigIn, sigInLength, &sigMax);
	 swt_min(sigIn, sigInLength, &sigMin);

	 for(count=0;count<sigInLength;count++)
     	 sigOut[count] = ((sigIn[count]-sigMin)/(sigMax-sigMin))*(maxv-minv) + minv;
	 return;
 }


 void wcodemat_abs(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv)
 {
     int count;
     double sigMin, sigMax;

	 swt_max_abs(sigIn, sigInLength, &sigMax);
	 swt_min_abs(sigIn, sigInLength, &sigMin);

	 for(count=0;count<sigInLength;count++)
     	 sigOut[count] = ((swt_abs(sigIn[count])-sigMin)/(sigMax-sigMin))*(maxv-minv) + minv;
	 return;
 }

 void wcodemat_matrix_row (double *matrixIn, int matrixInRow, int matrixInCol,
	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
						   int minv, int maxv, int abso)
 {
	 int count;
	 double *matrixOutTemp, *matrixOutPre;

	 matrixOutTemp = malloc(matrixInRow*matrixInCol*sizeof(double));
	 matrixOutPre = malloc(matrixInRow*matrixInCol*sizeof(double));

	 matrix_tran(matrixIn,matrixInCol,matrixInRow,
		         matrixOutTemp, matrixInRow, matrixInCol);
	 for(count=0;count<matrixInRow;count++)
	 {
		 if (abso)
             wcodemat_abs(matrixOutTemp+count*matrixInCol,matrixInCol,matrixOutPre+count*matrixInCol,matrixInCol,minv,maxv);
		 else
			 wcodemat(matrixOutTemp+count*matrixInCol,matrixInCol,matrixOutPre+count*matrixInCol,matrixInCol,minv,maxv);
	 }
	 free(matrixOutTemp);
	 //matrixOutM = malloc(matrixInRow*matrixInCol*sizeof(double));
	 matrix_tran(matrixOutPre,matrixInRow,matrixInCol,matrixOut,matrixInRow,matrixInCol);
	 free(matrixOutPre);

	 //for(count=0;count<matrixInRow*matrixInCol;count++)
	//	 matrixOut[count] = (int)(floor(matrixOutM[count]));
	 //free(matrixOutM);

	 return;
 }

 void wcodemat_matrix_col (double *matrixIn, int matrixInRow, int matrixInCol,
	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
						   int minv, int maxv, int abso)
 {
	 int count;
	 //double *matrixOutPre;

	 //matrixOutPre = malloc(matrixInRow*matrixInCol*sizeof(double));
	 for(count=0;count<matrixInCol;count++)
	 {
		 if (abso)
             wcodemat_abs(matrixIn+count*matrixInRow,matrixInRow,matrixOut+count*matrixInRow,matrixInRow,minv,maxv);
		 else
			 wcodemat(matrixIn+count*matrixInRow,matrixInRow,matrixOut+count*matrixInRow,matrixInRow,minv,maxv);
	 }
     //for(count=0;count<matrixInRow*matrixInCol;count++)
	//	 matrixOut[count] = (int)(floor(matrixOutPre[count]));
	 //free(matrixOutPre);

	 return;
 }

 void wcodemat_matrix (double *matrixIn, int matrixInRow, int matrixInCol,
	                   double *matrixOut, int matrixOutRow, int matrixOutCol,
					   int minv, int maxv, int abso)
 {
     //int count;
	 //double *matrixOutPre;

	 //matrixOutPre = malloc(matrixInRow*matrixInCol*sizeof(double));
	 //for(count=0;count<matrixInRow*matrixInCol;count++)
	 //{
		 if (abso)
             wcodemat_abs(matrixIn,matrixInRow*matrixInCol,matrixOut,matrixInRow*matrixInCol,minv,maxv);
		 else
			 wcodemat(matrixIn,matrixInRow*matrixInCol,matrixOut,matrixInRow*matrixInCol,minv,maxv);
	 //}
     //for(count=0;count<matrixInRow*matrixInCol;count++)
	//	 matrixOut[count] = (int)(floor(matrixOutPre[count]));
	 //free(matrixOutPre);
	 return;
 }

 void linspace(double lb, double ub, int n, double *sigOut, int sigOutLength)
 {
	 int count;
	 for(count=0;count<n;count++)
	 	 sigOut[count] = lb + count*(ub-lb)/(n-1);
	 return;
 }



/*-------------------------------------------
 * ocumsum
 *-----------------------------------------*/
/*void ocumsum (double *sigIn, int sigInLength)
{
	int i, j;
	for (i=1;i<sigInLength;i++)
	{
		sigIn[i]+=sigIn[i-1];
	}
	return;
}*/
