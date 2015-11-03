/*
 * -------------------------------------------------------------------------
 * swt.c -- 1-D and 2-D discrete stationary wavelet transform
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

 void swt_conv(double *sigIn, int sigInLength,
	           double *approx, int approxLength,
			   double *detail, int detailLength,
			   double *filterLow, double *filterHi,
			   int filterLength)
 {
	 i_conv(sigIn,sigInLength,approx,sigInLength,filterLow,filterLength);
     i_conv(sigIn,sigInLength,detail,sigInLength,filterHi,filterLength);
	 return;
 }

 void swt_out1 (double *sigIn, int sigInLength,
	            double *sigOutMatrix, int rowLength, int colLength,
				double *filterLow, double *filterHi, int filterLength, int step)
 {
     int row, count, *filterLen;
	 double *approxTemp,*approxOut,*sigOutTemp;
	 double *filterLowTemp, *filterHiTemp, *filterLowM, *filterHiM;

	 filterLen = malloc(step*sizeof(int));
	 approxTemp = malloc(colLength * sizeof(double));
	 approxOut = malloc(colLength*sizeof(double));
	 sigOutTemp = malloc(rowLength*colLength*sizeof(double));

	 filterLen[0] = filterLength;
	 if (step>1)
	 {
		 for(count=1;count<step;count++)
			 filterLen[count] = filterLen[count-1] * 2;
	 }

     filterLowTemp = malloc((filterLen[step-1]+1)*sizeof(double));
	 filterHiTemp = malloc((filterLen[step-1]+1)*sizeof(double));
	 filterLowM = malloc((filterLen[step-1]+1)*sizeof(double));
     filterHiM = malloc((filterLen[step-1]+1)*sizeof(double));
	 for(count=0;count<(filterLen[step-1]+1);count++)
	 {
		 filterLowM[count] = 0;
		 filterHiM[count] = 0;
	 }

	 verbatim_copy(filterLow,filterLength,filterLowTemp,filterLen[0]);
	 verbatim_copy(filterHi,filterLength,filterHiTemp,filterLen[0]);
	 verbatim_copy(sigIn,sigInLength,approxTemp,sigInLength);
     for (row=0;row<step;row++)
	 {
		 swt_conv(approxTemp,sigInLength,
			      approxOut,sigInLength,
				  sigOutTemp+row*colLength,sigInLength,
			      filterLowTemp,filterHiTemp,filterLen[row]);
		 verbatim_copy(approxOut,sigInLength,approxTemp,sigInLength);
		 if (row<(step-1))
			 {
				 dyadup_1D_feed_even (filterLowTemp, filterLen[row],
		                             filterLowM, filterLen[row]*2+1);
		         dyadup_1D_feed_even (filterHiTemp, filterLen[row],
		                             filterHiM, filterLen[row]*2+1);
			     verbatim_copy(filterLowM+1,filterLen[row]*2,filterLowTemp,filterLen[row]*2);
			     verbatim_copy(filterHiM+1,filterLen[row]*2,filterHiTemp,filterLen[row]*2);
		      }
	  }
	  verbatim_copy(approxTemp,sigInLength,sigOutTemp+step*colLength,sigInLength);
	  free(approxTemp);
	  free(approxOut);
	  free(filterLowM);
	  free(filterHiM);
	  free(filterLowTemp);
	  free(filterHiTemp);
	  free(filterLen);
      matrix_tran (sigOutTemp, rowLength, colLength,
	               sigOutMatrix, colLength, rowLength);
	  free(sigOutTemp);
	  return;
 }

 void swt_out2 (double *sigIn, int sigInLength,
	            double *approxMatrix, double *detailMatrix,
				int rowLength, int colLength, double *filterLow,
				double *filterHi, int filterLength, int step)
 {
     int row, count, *filterLen;
	 double *approxMatrixTemp, *detailMatrixTemp, *approxTemp;
	 double *filterLowTemp, *filterHiTemp, *filterLowM, *filterHiM;


     filterLen = malloc(step*sizeof(int));
	 approxTemp = malloc(colLength * sizeof(double));
	 approxMatrixTemp = malloc(rowLength*colLength*sizeof(double));
	 detailMatrixTemp = malloc(rowLength*colLength*sizeof(double));

	 filterLen[0] = filterLength;
	 if (step>1)
	 {
		 for(count=1;count<step;count++)
			 filterLen[count] = filterLen[count-1] * 2;
	 }

	 filterLowTemp = malloc((filterLen[step-1])*sizeof(double));
	 filterHiTemp = malloc((filterLen[step-1])*sizeof(double));
	 filterLowM = malloc((filterLen[step-1])*sizeof(double));
     filterHiM = malloc((filterLen[step-1])*sizeof(double));
	 for(count=0;count<(filterLen[step-1]);count++)
	 {
		 filterLowM[count] = 0;
		 filterHiM[count] = 0;
	 }

	 verbatim_copy(filterLow,filterLength,filterLowTemp,filterLen[0]);
	 verbatim_copy(filterHi,filterLength,filterHiTemp,filterLen[0]);
	 verbatim_copy(sigIn,sigInLength,approxTemp,sigInLength);
     for (row=0;row<step;row++)
	 {
		 swt_conv(approxTemp,sigInLength,
			      approxMatrixTemp+row*colLength,sigInLength,
				  detailMatrixTemp+row*colLength,sigInLength,
			      filterLowTemp,filterHiTemp,filterLen[row]);
		 verbatim_copy(approxMatrixTemp+row*colLength,sigInLength,approxTemp,sigInLength);
		 if (row<(step-1))
			 {
				 for(count=0;count<filterLen[row];count++)
				 {
					 filterLowM[2*count] = filterLowTemp[count];
					 filterLowM[2*count+1] = 0;
                     filterHiM[2*count] = filterHiTemp[count];
					 filterHiM[2*count+1] = 0;
				 }

			     verbatim_copy(filterLowM,filterLen[row]*2,filterLowTemp,filterLen[row]*2);
			     verbatim_copy(filterHiM,filterLen[row]*2,filterHiTemp,filterLen[row]*2);
		      }
	 }
	 free(approxTemp);
	 matrix_tran (approxMatrixTemp, rowLength, colLength,
	              approxMatrix, colLength, rowLength);
	 matrix_tran (detailMatrixTemp, rowLength, colLength,
	              detailMatrix, colLength, rowLength);
	 free(approxMatrixTemp);
     free(detailMatrixTemp);
	 free(filterLowM);
	 free(filterHiM);
	 free(filterLowTemp);
	 free(filterHiTemp);
	 free(filterLen);
	 return;
 }

 void swt2_output4(double *matrixIn, int matrixInRow, int matrixInCol,
	               double *matrixOutApprox, double *matrixOutColDetail,
				   double *matrixOutRowDetail, double *matrixOutDetail,
				   int matrixOutRow, int matrixOutCol,
				   double *filterLow, double *filterHi,
				   int filterLength, int step)
 {
	 int row, col, count;
	 int filterLen, twt;
	 double *matrixOutApproxPre, *matrixOutDetailPre;
	 double *filterLowM, *filterHiM;
	 double *pMatrixOutApproxPre, *pMatrixOutDetailPre;

	 //matrixInPre = malloc(matrixInRow*matrixInCol*sizeof(double));
	 //matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       //matrixInPre, matrixInRow, matrixInCol);

	 matrixOutApproxPre = malloc (matrixInRow *
			       matrixOutCol * sizeof (double));
     matrixOutDetailPre = malloc (matrixInRow *
			       matrixOutCol * sizeof (double));

	 swt_exp2(step,&twt);
	 filterLen = filterLength*twt;
     filterLowM = malloc(filterLen*sizeof(double));
     filterHiM = malloc(filterLen*sizeof(double));
	 for(count=0;count<filterLen;count++)
	 {
		 if (count%twt==0)
		 {
			 filterLowM[count] = filterLow[count/twt];
			 filterHiM[count] = filterHi[count/twt];
		 }
		 else
		 {
			 filterLowM[count] = 0;
			 filterHiM[count] = 0;
		 }
	 }


	 for (row=0;row<matrixInRow;row++)
	 {
		 swt_conv(matrixIn+row*matrixInCol,matrixInCol,
			      matrixOutApproxPre+row*matrixInCol,matrixInCol,
				  matrixOutDetailPre+row*matrixInCol,matrixInCol,
			      filterLowM,filterHiM,filterLen);

	 }

	 //free (matrixInPre);
	 pMatrixOutApproxPre = malloc (matrixInRow *
				matrixOutCol * sizeof (double));
     matrix_tran (matrixOutApproxPre, matrixInRow, matrixInCol,
	       pMatrixOutApproxPre, matrixInRow, matrixInCol);
     free (matrixOutApproxPre);

     pMatrixOutDetailPre = malloc (matrixInRow *
				matrixOutCol * sizeof (double));
     matrix_tran (matrixOutDetailPre, matrixInRow, matrixInCol,
	       pMatrixOutDetailPre, matrixInRow, matrixInCol);

     free (matrixOutDetailPre);

	 for (col = 0; col < matrixInCol; col++)
     {
		swt_conv(pMatrixOutApproxPre+col*matrixInRow,matrixInRow,
			      matrixOutApprox+col*matrixInRow,matrixInRow,
				  matrixOutColDetail+col*matrixInRow,matrixInRow,
			      filterLowM,filterHiM,filterLen);
	 }
     free (pMatrixOutApproxPre);

	 for (col = 0; col < matrixInCol; col++)
     {
		swt_conv(pMatrixOutDetailPre+col*matrixInRow,matrixInRow,
			      matrixOutRowDetail+col*matrixInRow,matrixInRow,
				  matrixOutDetail+col*matrixInRow,matrixInRow,
			      filterLowM,filterHiM,filterLen);
	 }
     free (pMatrixOutDetailPre);
	 return;
 }

 void swt2_output4_step(double *matrixIn, int matrixInRow, int matrixInCol,
	               double *matrixOutApprox, double *matrixOutColDetail,
				   double *matrixOutRowDetail, double *matrixOutDetail,
				   int matrixOutRow, int matrixOutCol,
				   double *filterLow, double *filterHi,
				   int filterLength, int step)
 {
     int count;
	 double *matrixApproxTemp,*matrixApproxOut;

	 matrixApproxTemp = malloc(matrixInRow*matrixInCol*sizeof(double));
     matrixApproxOut = malloc(matrixInRow*matrixInCol*sizeof(double));
	 matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixApproxTemp, matrixInRow, matrixInCol);

	 for(count=0;count<step;count++)
	 {
		 swt2_output4(matrixApproxTemp, matrixInRow, matrixInCol,
	               matrixApproxOut, matrixOutColDetail+count*matrixInRow*matrixInCol,
				   matrixOutRowDetail+count*matrixInRow*matrixInCol,
				   matrixOutDetail+count*matrixInRow*matrixInCol,
				   matrixOutRow, matrixOutCol, filterLow, filterHi,
				   filterLength, count);
		 verbatim_copy(matrixApproxOut,matrixInRow*matrixInCol,
			           matrixOutApprox+count*matrixInRow*matrixInCol,
					   matrixInRow*matrixInCol);
         matrix_tran (matrixApproxOut, matrixInCol, matrixInRow,
	                  matrixApproxTemp, matrixInRow, matrixInCol);
	 }
	 free(matrixApproxOut);
	 free(matrixApproxTemp);

	 return;
 }

 void swt2_output1_step(double *matrixIn, int matrixInRow,
	                    int matrixInCol,  double *matrixOut,
				        int matrixOutRow, int matrixOutCol,
				        double *filterLow, double *filterHi,
				        int filterLength, int step)
 {
     int count;
	 double *matrixApproxTemp,*matrixApproxOut;

	 matrixApproxTemp = malloc(matrixInRow*matrixInCol*sizeof(double));
     matrixApproxOut = malloc(matrixInRow*matrixInCol*sizeof(double));
	 matrix_tran (matrixIn, matrixInCol, matrixInRow,
	       matrixApproxTemp, matrixInRow, matrixInCol);

	 for(count=0;count<step;count++)
	 {
		 swt2_output4(matrixApproxTemp, matrixInRow, matrixInCol,
	               matrixApproxOut, matrixOut+count*matrixInRow*matrixInCol,
				   matrixOut+step*matrixInRow*matrixInCol+count*matrixInRow*matrixInCol,
				   matrixOut+2*step*matrixInRow*matrixInCol+count*matrixInRow*matrixInCol,
				   matrixInRow, matrixOutCol, filterLow, filterHi,
				   filterLength, count);

         matrix_tran (matrixApproxOut, matrixInCol, matrixInRow,
	                  matrixApproxTemp, matrixInRow, matrixInCol);
	 }
	 verbatim_copy(matrixApproxOut,matrixInRow*matrixInCol,
			       matrixOut+3*step*matrixInRow*matrixInCol,
				   matrixInRow*matrixInCol);
	 free(matrixApproxOut);
	 free(matrixApproxTemp);

	 return;
 }

 void iswt_conv (double *approx, double *detail, int sigInLength,
	             double *sigOut, int sigOutLength, double *filterLow,
				 double *filterHi, int filterLength)
 {
     int length1, length2, count;
	 double *approxTempOdd, *detailTempOdd;
	 double *approxTempEven, *detailTempEven;
	 double *approxOutOdd, *detailOutOdd;
	 double *approxOutEven, *detailOutEven;
	 double *approxT, *detailT;
	 double *sigOdd, *sigEven;

	 length1 = (int)floor(sigInLength/2.0);
	 approxTempOdd = malloc(length1*sizeof(double));
     detailTempOdd = malloc(length1*sizeof(double));
     approxTempEven = malloc(length1*sizeof(double));
     detailTempEven = malloc(length1*sizeof(double));

     dyaddown_1D_keep_odd (approx, sigInLength, approxTempOdd, length1);
     dyaddown_1D_keep_even (approx, sigInLength, approxTempEven, length1);
     dyaddown_1D_keep_odd (detail, sigInLength, detailTempOdd, length1);
     dyaddown_1D_keep_even (detail, sigInLength, detailTempEven, length1);

	 length2 = 2 * length1;
	 approxOutOdd = malloc(length2*sizeof(double));
	 approxOutEven = malloc(length2*sizeof(double));
	 detailOutOdd = malloc(length2*sizeof(double));
	 detailOutEven = malloc(length2*sizeof(double));

	 for (count=0;count<length1;count++)
	 {
          approxOutOdd[2*count] = approxTempOdd[count];
		  detailOutOdd[2*count] = detailTempOdd[count];
		  approxOutOdd[2*count+1] = 0;
		  detailOutOdd[2*count+1] = 0;
		  approxOutEven[2*count] = 0;
		  detailOutEven[2*count] = 0;
		  approxOutEven[2*count+1] = approxTempEven[count];
		  detailOutEven[2*count+1] = detailTempEven[count];
	 }

     free(approxTempOdd);
	 free(detailTempOdd);
	 free(approxTempEven);
	 free(detailTempEven);

	 approxT = malloc(length2*sizeof(double));
	 detailT = malloc(length2*sizeof(double));
	 sigOdd = malloc(length2*sizeof(double));
	 sigEven = malloc(length2*sizeof(double));

	 i_conv(approxOutOdd,length2,approxT,length2,filterLow,filterLength);
     i_conv(detailOutOdd,length2,detailT,length2,filterHi,filterLength);
     for(count=0;count<length2;count++)
		 sigOdd[count] = approxT[count] + detailT[count];
	 free(approxOutOdd);
	 free(detailOutOdd);

	 i_conv(approxOutEven,length2,approxT,length2,filterLow,filterLength);
     i_conv(detailOutEven,length2,detailT,length2,filterHi,filterLength);
     for(count=0;count<length2;count++)
		 sigEven[count] = approxT[count] + detailT[count];
	 free(approxOutEven);
	 free(detailOutEven);
	 free(approxT);
	 free(detailT);

	 for(count=(sigOutLength-filterLength-1);count<sigOutLength;count++)
		 sigOut[count+1+filterLength-sigOutLength] = (sigOdd[count] + sigEven[count]) / 2.0;
     for(count=0;count<(sigOutLength - filterLength - 1);count++)
		 sigOut[count+filterLength+1] = (sigOdd[count]+sigEven[count])/2.0;
	 free(sigOdd);
	 free(sigEven);

	 return;
 }

 void iswt_conv_step (double *approx, double *detail, int sigInLength,
	                  double *sigOut, int sigOutLength, double *filterLow,
				      double *filterHi, int filterLength, int level)
 {
     int length1, length2, count, twt, mu, su, startpoint;
	 double *approxTempOdd, *detailTempOdd;
	 double *approxTempEven, *detailTempEven;
	 double *approxOutOdd, *detailOutOdd;
	 double *approxOutEven, *detailOutEven;
	 double *approxT, *detailT;
	 double *sigOdd, *sigEven;
	 double *filterLowTemp, *filterHiTemp;

	 length1 = (int)floor(sigInLength/2.0);
	 approxTempOdd = malloc(length1*sizeof(double));
     detailTempOdd = malloc(length1*sizeof(double));
     approxTempEven = malloc(length1*sizeof(double));
     detailTempEven = malloc(length1*sizeof(double));

     dyaddown_1D_keep_odd (approx, sigInLength, approxTempOdd, length1);
     dyaddown_1D_keep_even (approx, sigInLength, approxTempEven, length1);
     dyaddown_1D_keep_odd (detail, sigInLength, detailTempOdd, length1);
     dyaddown_1D_keep_even (detail, sigInLength, detailTempEven, length1);

	 length2 = 2 * length1;
	 approxOutOdd = malloc(length2*sizeof(double));
	 approxOutEven = malloc(length2*sizeof(double));
	 detailOutOdd = malloc(length2*sizeof(double));
	 detailOutEven = malloc(length2*sizeof(double));

	 for (count=0;count<length1;count++)
	 {
          approxOutOdd[2*count] = approxTempOdd[count];
		  detailOutOdd[2*count] = detailTempOdd[count];
		  approxOutOdd[2*count+1] = 0;
		  detailOutOdd[2*count+1] = 0;
		  approxOutEven[2*count] = 0;
		  detailOutEven[2*count] = 0;
		  approxOutEven[2*count+1] = approxTempEven[count];
		  detailOutEven[2*count+1] = detailTempEven[count];
	 }

     free(approxTempOdd);
	 free(detailTempOdd);
	 free(approxTempEven);
	 free(detailTempEven);

	 swt_exp2(level-1,&twt);
	 if (level==1)
	 {
		 mu = 1;
		 su = 0;
	 }
	 else
	 {
		 mu = twt;
		 su = twt - 1;
	 }

	 filterLowTemp = malloc(mu*filterLength*sizeof(double));
	 filterHiTemp = malloc(mu*filterLength*sizeof(double));

     for(count=0;count<mu*filterLength;count++)
	 {
		 filterLowTemp[count] = 0;
		 filterHiTemp[count] = 0;
	 }

     for(count=0;count<filterLength;count++)
	 {
		 filterLowTemp[mu*count] = filterLow[count];
		 filterHiTemp[mu*count] = filterHi[count];
	 }

	 approxT = malloc(length2*sizeof(double));
	 detailT = malloc(length2*sizeof(double));
	 sigOdd = malloc(length2*sizeof(double));
	 sigEven = malloc(length2*sizeof(double));

	 i_conv(approxOutOdd,length2,approxT,length2,filterLowTemp,mu*filterLength);
     i_conv(detailOutOdd,length2,detailT,length2,filterHiTemp,mu*filterLength);
     for(count=0;count<length2;count++)
		 sigOdd[count] = approxT[count] + detailT[count];
	 free(approxOutOdd);
	 free(detailOutOdd);

	 i_conv(approxOutEven,length2,approxT,length2,filterLowTemp,mu*filterLength);
     i_conv(detailOutEven,length2,detailT,length2,filterHiTemp,mu*filterLength);
     for(count=0;count<length2;count++)
		 sigEven[count] = approxT[count] + detailT[count];
	 free(approxOutEven);
	 free(detailOutEven);
	 free(approxT);
	 free(detailT);
	 free(filterLowTemp);
	 free(filterHiTemp);


	 startpoint = sigInLength - filterLength * mu - su - 1;

	 for(count=(startpoint);count<sigOutLength;count++)
		 sigOut[count-startpoint] = (sigOdd[count] + sigEven[count]) / 2.0;
     for(count=0;count<startpoint;count++)
		 sigOut[count+1+filterLength*mu+su] = (sigOdd[count]+sigEven[count])/2.0;
	 free(sigOdd);
	 free(sigEven);

	 return;
 }

 void iswt_input1 (double *matrixIn, int rowLength, int colLength,
	              double *sigOut, int sigOutLength, double *filterLow,
				  double *filterHi, int filterLength)
 {
     int step, count;
	 double *approxTemp, *approxOut, *matrixInTemp;

	 matrixInTemp = malloc(rowLength*colLength*sizeof(double));
	 matrix_tran(matrixIn, colLength, rowLength,
		 matrixInTemp, colLength, rowLength);

	 step  = rowLength - 1;
	 approxTemp = malloc(colLength*sizeof(double));
	 approxOut = malloc(colLength*sizeof(double));


     verbatim_copy(matrixInTemp+step*colLength,colLength,approxTemp,colLength);

	 for(count=0;count<step;count++)
	 {
		 iswt_conv_step(approxTemp,matrixInTemp+(step-count-1)*colLength,colLength,
			       approxOut,colLength,filterLow,filterHi,filterLength,step-count);
		 verbatim_copy(approxOut,colLength,approxTemp,colLength);
	 }
	 verbatim_copy(approxTemp,colLength,sigOut, colLength);
	 free(matrixInTemp);
	 free(approxTemp);
	 free(approxOut);

	 return;
 }



 void iswt_input2 (double *matrixApproxIn, double *matrixDetailIn,
	               int rowLength, int colLength,
	              double *sigOut, int sigOutLength, double *filterLow,
				  double *filterHi, int filterLength)
 {
	 int count, step;
	 double *approxTemp, *approxOut;
	 double *matrixApproxInTemp, *matrixDetailInTemp;

     step = rowLength;
	 approxTemp = malloc(colLength*sizeof(double));
	 approxOut = malloc(colLength*sizeof(double));
	 matrixApproxInTemp = malloc(rowLength*colLength*sizeof(double));
	 matrixDetailInTemp = malloc(rowLength*colLength*sizeof(double));

	 matrix_tran(matrixApproxIn, colLength, rowLength,
		         matrixApproxInTemp,colLength, rowLength);
     matrix_tran(matrixDetailIn, colLength, rowLength,
		         matrixDetailInTemp,colLength, rowLength);

     verbatim_copy(matrixApproxInTemp+(step-1)*colLength,colLength,approxTemp,colLength);
     for(count=0;count<step;count++)
	 {
		 iswt_conv_step(approxTemp,matrixDetailInTemp+(step-count-1)*colLength,colLength,
			       approxOut,colLength,filterLow,filterHi,filterLength,step-count);
		 verbatim_copy(approxOut,colLength,approxTemp,colLength);
	 }
	 verbatim_copy(approxTemp,colLength,sigOut, colLength);
	 free(approxTemp);
	 free(approxOut);
	 free(matrixApproxInTemp);
	 free(matrixDetailInTemp);

	 return;
 }


void iswt2(double *matrixInApprox, double *matrixInColDetail,
		   double *matrixInRowDetail, double *matrixInDetail,
		   int matrixInRow, int matrixInCol,
		   double *matrixOut, int matrixOutRow, int matrixOutCol,
		   double *filterLow, double *filterHi,
		   int filterLength, int step)
 {
	 int row, col;//, count;
	 //int filterLen, twt;
	 double *matrixOutApproxPre, *matrixOutDetailPre;
	 double *matrixOutPre;
	 double *pMatrixOutApproxPre, *pMatrixOutDetailPre;


	 matrixOutApproxPre = malloc (matrixInRow *
			       matrixInCol * sizeof (double));
     matrixOutDetailPre = malloc (matrixInRow *
			       matrixInCol * sizeof (double));


	 for (col=0;col<matrixInCol;col++)
	 {
		 iswt_conv_step (matrixInApprox+col*matrixInRow, matrixInColDetail+col*matrixInRow,
			             matrixInRow, matrixOutApproxPre+col*matrixInRow, matrixInRow,
						 filterLow, filterHi, filterLength, step);
		 iswt_conv_step (matrixInRowDetail+col*matrixInRow, matrixInDetail+col*matrixInRow,
			             matrixInRow, matrixOutDetailPre+col*matrixInRow, matrixInRow,
						 filterLow, filterHi, filterLength, step);
	 }


	 pMatrixOutApproxPre = malloc (matrixInRow *
				matrixInCol * sizeof (double));
     matrix_tran (matrixOutApproxPre, matrixInCol, matrixInRow,
	       pMatrixOutApproxPre, matrixInRow, matrixInCol);
     free (matrixOutApproxPre);

     pMatrixOutDetailPre = malloc (matrixInRow *
				matrixInCol * sizeof (double));
     matrix_tran (matrixOutDetailPre, matrixInCol, matrixInRow,
	       pMatrixOutDetailPre, matrixInRow, matrixInCol);

     free (matrixOutDetailPre);

	 matrixOutPre = malloc(matrixInRow*matrixInCol*sizeof(double));
	 for (row = 0; row < matrixInRow; row++)
     {
		iswt_conv_step (pMatrixOutApproxPre+row*matrixInCol, pMatrixOutDetailPre+row*matrixInCol,
			             matrixInCol, matrixOutPre+row*matrixInCol, matrixInCol,
						 filterLow, filterHi, filterLength, step);
	 }
     free (pMatrixOutApproxPre);
     free (pMatrixOutDetailPre);

     matrix_tran (matrixOutPre, matrixInRow, matrixInCol,
	              matrixOut, matrixInRow, matrixInCol);
	 free(matrixOutPre);

	 return;
 }


void iswt2_input4_step(double *matrixInApprox, double *matrixInColDetail,
		               double *matrixInRowDetail, double *matrixInDetail,
		               int matrixInRow, int matrixInCol,
		               double *matrixOut, int matrixOutRow, int matrixOutCol,
		               double *filterLow, double *filterHi,
		               int filterLength, int step)
{
    int count;
	double *matrixApproxTemp, *matrixApproxOut;

	matrixApproxTemp = malloc(matrixInRow*matrixInCol*sizeof(double));
 	matrixApproxOut = malloc(matrixInRow*matrixInCol*sizeof(double));

    verbatim_copy(matrixInApprox+(step-1)*matrixInRow*matrixInCol,
		          matrixInRow*matrixInCol,
				  matrixApproxTemp, matrixInRow*matrixInCol);

	for(count=0;count<step;count++)
	{
		iswt2(matrixApproxTemp, matrixInColDetail+(step-count-1)*matrixInRow*matrixInCol,
		   matrixInRowDetail+(step-count-1)*matrixInRow*matrixInCol,
		   matrixInDetail+(step-count-1)*matrixInRow*matrixInCol,
		   matrixInRow, matrixInCol, matrixApproxOut, matrixOutRow, matrixOutCol,
		   filterLow, filterHi, filterLength, step-count);
		if (count!=(step-1))
		{
			verbatim_copy(matrixApproxOut, matrixInRow*matrixInCol,
				  matrixApproxTemp, matrixInRow*matrixInCol);
		}
	}
	free(matrixApproxTemp);

    verbatim_copy(matrixApproxOut, matrixInRow*matrixInCol,
				  matrixOut, matrixInRow*matrixInCol);
    free(matrixApproxOut);

	return;
}

void iswt2_input1_step(double *matrixIn,  int matrixInRow, int matrixInCol,
		               double *matrixOut, int matrixOutRow, int matrixOutCol,
		               double *filterLow, double *filterHi,
		               int filterLength, int step)
{
    int count;
	double *matrixApproxTemp, *matrixApproxOut;

	matrixApproxTemp = malloc(matrixInRow*matrixInCol*sizeof(double));
 	matrixApproxOut = malloc(matrixInRow*matrixInCol*sizeof(double));

    verbatim_copy(matrixIn+3*step*matrixInRow*matrixInCol,
		          matrixInRow*matrixInCol,
				  matrixApproxTemp, matrixInRow*matrixInCol);

	for(count=0;count<step;count++)
	{
		iswt2(matrixApproxTemp, matrixIn+(step-count-1)*matrixInRow*matrixInCol,
		   matrixIn+(step-count-1)*matrixInRow*matrixInCol+step*matrixInRow*matrixInCol,
		   matrixIn+(step-count-1)*matrixInRow*matrixInCol+2*step*matrixInRow*matrixInCol,
		   matrixInRow, matrixInCol, matrixApproxOut, matrixOutRow, matrixOutCol,
		   filterLow, filterHi, filterLength, step-count);
		if (count!=(step-1))
		{
			verbatim_copy(matrixApproxOut, matrixInRow*matrixInCol,
				  matrixApproxTemp, matrixInRow*matrixInCol);
		}
	}
	free(matrixApproxTemp);

    verbatim_copy(matrixApproxOut, matrixInRow*matrixInCol,
				  matrixOut, matrixInRow*matrixInCol);
    free(matrixApproxOut);

	return;
}
