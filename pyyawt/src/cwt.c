/*
 * -------------------------------------------------------------------------
 * cwt.c -- continuous wavelet transform
 * PYYAWT - yet another python wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
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
  #include "kiss_fft.h"


/*void haar(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;

	for(count=0;count<sigInLength;count++)
	{
		if (x[count]<0)
			psi[count] = -1/ys;
		else
			psi[count] = 1/ys;
	}
	return;
}*/

double
powof (double x, double alpha)
{
  double         resu;


  if (x >= 0) /* in this case, no problem */
    {
      if (x == 0)
	resu = 0.0;
      else
	resu = exp (alpha * log (x));
    }
  else /* there may be problems */
    {
      if (alpha == (int)(alpha)) /* if alpha is an integer */
	{
	  /* then x^alpha is real-valued */
	  resu = powof ( -x, alpha);
	  /* and the sign of the results depends on
	     wether alpha is ODD or EVEN */
	  if (ISODD(alpha) == 1)
	    {
	      /* alpha is ODD */
	      resu = -resu;
	    }
	}
      else
	{
	  //Scierror (999,"Attempt to compute x^alpha with x<0 : complex valued result\n");
	  return 0;
	}
    }
  return resu;
}

/*-------------------------------------------
 * Sinus Scale Filter Generation
 *-----------------------------------------*/

void sinus(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;

	for(count=0;count<sigInLength;count++)
		psi[count] = -1*sin(2*PI*x[count])/ys;
	return;
}

/*-------------------------------------------
 * Poisson Scale Filter Generation
 *-----------------------------------------*/

void poisson(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count] * 4.0 / 9.0;
		psi[count] = (1-x2)/(PI*(1+x2)*ys);
	}
	return;
}

/*-------------------------------------------
 * Mexican Hat Scale Filter Generation
 *-----------------------------------------*/

void mexihat(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2, con;

    //con = 2*sqrt(sqrt(PI))/sqrt(3);
    con = 2/(sqrt(3)*sqrt(sqrt(PI)));
	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = (1-x2)*exp(-x2/2)*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Morlet Scale Filter Generation
 *-----------------------------------------*/

void morlet(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = cos(5*x[count])*exp(-x2/2)/ys;
	}
	return;
}

/*-------------------------------------------
 * DoG Scale Filter Generation
 *-----------------------------------------*/

void DOGauss(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = exp(-x2/2)/ys-exp(-x2/8)/(2*ys);
	}
	return;
}

/*-------------------------------------------
 * Gauss Scale Filter Generation
 *-----------------------------------------*/

void Gauss(double *x, int sigInLength, double *psi, int sigOutLength, int n, double ys)
{
	switch (n) {
		case 1:
			Gaus1(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 2:
			Gaus2(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 3:
			Gaus3(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 4:
			Gaus4(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 5:
			Gaus5(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 6:
			Gaus6(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 7:
			Gaus7(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 8:
			Gaus8(x, sigInLength, psi, sigOutLength, ys);
			break;
		default:
			break;
	}


	return;
}

void Gaus1(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = -2*x[count]*exp(-x2)/sqrt(sqrt(PI/2));
	}
	return;
}


void Gaus2(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = 2*(2*x2-1)*exp(-x2)/sqrt(3*sqrt(PI/2));
	}
	return;
}


void Gaus3(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3;

	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		x3 = x2 * x[count];
		psi[count] = -4*(2*x3-3*x[count])*exp(-x2)/sqrt(15*sqrt(PI/2));
	}
	return;
}


void Gaus4(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4;

	for(count=0;count<sigInLength;count++)
	{
	 	x2 = x[count] * x[count];
		x4 = x2 * x2;
		psi[count] = 4*(-12*x2+4*x4+3)*exp(-x2)/sqrt(105*sqrt(PI/2));
	}
	return;
}


void Gaus5(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3,x5;

	for(count=0;count<sigInLength;count++)
	{
         x2 = x[count] * x[count];
	 	 x3 = x2* x[count];
		 x5 = x3 * x2;
		psi[count] = 8*(-4*x5+20*x3-15*x[count])*exp(-x2)/sqrt(105*9*sqrt(PI/2));
	}
	return;
}


void Gaus6(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4,x6;

	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
		x4 = x2 * x2;
		x6 = x4 * x2;
		psi[count] = 8*(8*x6-60*x4+90*x2-15)*exp(-x2)/sqrt(105*9*11*sqrt(PI/2));
	}
	return;
}


void Gaus7(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3,x5,x7;

	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
	    x3 = x2 * x[count];
		x5 = x3 * x2;
		x7 = x5 * x2;
		psi[count] = 16*(-8*x7+84*x5-210*x3+105*x[count])*exp(-x2)/sqrt(105*9*11*13*sqrt(PI/2));
	}
	return;
}


void Gaus8(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4,x6,x8;

	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
    	x4 = x2 * x2;
	    x6 = x4 * x2;
		x8 = x6 * x2;
		psi[count] = 16*(16*x8-224*x6+840*x4-840*x2+105)*exp(-x2)/sqrt(105*9*11*13*15*sqrt(PI/2));
	}
	return;
}

/*-------------------------------------------
 * Complex Gauss Scale Filter Generation
 *-----------------------------------------*/
void cgauss(double *x, int sigInLength, int p, double *psir, double *psii, int sigOutLength, double ys)
{

    switch (p) {
		case 1:
			cgau1(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 2:
			cgau2(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 3:
			cgau3(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 4:
			cgau4(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 5:
			cgau5(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 6:
			cgau6(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 7:
			cgau7(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 8:
			cgau8(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		default:
			break;
	}
	return;
}

void cgau1(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-2*x[count]*cosx-sinx)*exp(-x2)/sqrt(2*sqrt(PI/2));
	   psii[count]=(2*x[count]*sinx-cosx)*exp(-x2)/sqrt(2*sqrt(PI/2));
	}
	return;
}

void cgau1_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau1(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau2(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(4*x2*cosx+4*x[count]*sinx-3*cosx)*exp(-x2)/sqrt(10*sqrt(PI/2));
	   psii[count]=(-4*x2*sinx+4*x[count]*cosx+3*sinx)*exp(-x2)/sqrt(10*sqrt(PI/2));
	}
	return;
}

void cgau2_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau2(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau3(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-8*x3*cosx-12*x2*sinx+18*x[count]*cosx+7*sinx)*exp(-x2)/sqrt(76*sqrt(PI/2));
	   psii[count]=(8*x3*sinx-12*x2*cosx-18*x[count]*sinx+7*cosx)*exp(-x2)/sqrt(76*sqrt(PI/2));
	}
	return;
}

void cgau3_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau3(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau4(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
 	   psir[count]=(16*x4*cosx+32*x3*sinx-72*x2*cosx-56*x[count]*sinx+25*cosx)*exp(-x2)/sqrt(764*sqrt(PI/2));
 	   psii[count]=(-16*x4*sinx+32*x3*cosx+72*x2*sinx-56*x[count]*cosx-25*sinx)*exp(-x2)/sqrt(764*sqrt(PI/2));
	}
	return;
}

void cgau4_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau4(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau5(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
 	   psir[count]=(-32*x5*cosx-80*x4*sinx+240*x3*cosx+280*x2*sinx-250*x[count]*cosx-81*sinx)*exp(-x2)/sqrt(9496*sqrt(PI/2));
 	   psii[count]=(32*x5*sinx-80*x4*cosx-240*x3*sinx+280*x2*cosx+250*x[count]*sinx-81*cosx)*exp(-x2)/sqrt(9496*sqrt(PI/2));
	}
	return;
}

void cgau5_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau5(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau6(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(64*x6*cosx+192*x5*sinx-720*x4*cosx-1120*x3*sinx+1500*x2*cosx+972*x[count]*sinx-331*cosx)*exp(-x2)/sqrt(140152*sqrt(PI/2));
	   psii[count]=(-64*x6*sinx+192*x5*cosx+720*x4*sinx-1120*x3*cosx-1500*x2*sinx+972*x[count]*cosx+331*sinx)*exp(-x2)/sqrt(140152*sqrt(PI/2));
	}
	return;
}

void cgau6_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau6(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau7(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6, x7;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   x7 = x4 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-128*x7*cosx-448*x6*sinx+2016*x5*cosx+3920*x4*sinx-7000*x3*cosx-6804*x2*sinx+4634*x[count]*cosx+1303*sinx)*exp(-x2)/sqrt(2390480*sqrt(PI/2));
	   psii[count]=(128*x7*sinx-448*x6*cosx-2016*x5*sinx+3920*x4*cosx+7000*x3*sinx-6804*x2*cosx-4634*x[count]*sinx+1303*cosx)*exp(-x2)/sqrt(2390480*sqrt(PI/2));
	}
	return;
}

void cgau7_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau7(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau8(double *x, int sigInLength,
			double *psir, double *psii,
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6, x7, x8;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   x7 = x4 * x3;
	   x8 = x4 * x4;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(256*x8*cosx+1024*x7*sinx-5376*x6*cosx-12544*x5*sinx+28000*x4*cosx+36288*x3*sinx-37072*x2*cosx-20848*x[count]*sinx+5937*cosx)*exp(-x2)/sqrt(46206736*sqrt(PI/2));
	   psii[count]=(-256*x8*sinx+1024*x7*cosx+5376*x6*sinx-12544*x5*cosx-28000*x4*sinx+36288*x3*cosx+37072*x2*sinx-20848*x[count]*cosx-5937*sinx)*exp(-x2)/sqrt(46206736*sqrt(PI/2));
	}
	return;
}

void cgau8_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	cgau8(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}


/*-------------------------------------------
 * Complex Morlet Scale Filter Generation
 *-----------------------------------------*/

void cmorlet(double *x, int sigInLength,
			double fb, double fc, double *psir, double *psii,
			int sigOutLength, double ys)
{
	int count;
	double con, econ, x2;

	con = 1/sqrt(PI*fb);
	for(count=0;count<sigInLength;count++)
	{
		x2=x[count]*x[count];
		econ = exp(-x2/fb);
		psir[count] = cos(2*PI*fc*x[count])*econ*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*econ*con/ys;
	}
	return;
}

void cmorlet_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ, x2;

	con = 1/sqrt(PI);
	for(count=0;count<sigInLength;count++)
	{
		x2=x[count]*x[count];
		econ = exp(-x2);
		f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Shannon Scale Filter Generation
 *-----------------------------------------*/

void shanwavf(double *x, int sigInLength,
			double fb, double fc, double *psir, double *psii,
			int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = sqrt(fb);
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*fb*PI)/(x[count]*fb*PI);
		else
			econ = 1;
		psir[count] = cos(2*PI*fc*x[count])*econ*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*econ*con/ys;
	}
	return;
}

void shanwavf_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = 1;
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*PI)/(x[count]*PI);
		else
			econ = 1;
		f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Frequency B-Spline Scale Filter Generation
 *-----------------------------------------*/
void fbspwavf(double *x, int sigInLength,int m,
			double fb, double fc, double *psir, double *psii,
			int sigOutLength, double ys)
{
	int count, i;
	double con, econ;

	con = sqrt(fb);
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*fb*PI/(double)m)/(x[count]*fb*PI/(double)m);
		else
			econ = 1;
//         for(i=0;i<m;i++)
// 			econ*=econ;

		psir[count] = cos(2*PI*fc*x[count])*powof(econ,m)*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*powof(econ,m)*con/ys;
	}
	return;
}

void fbspwavf_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = 1;
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*PI)/(x[count]*PI);
		else
			econ = 1;
     	f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Cauchy Scale Filter Generation
 *-----------------------------------------*/

void cauchy(double *x, int sigInLength,
			double fb, double fc, double *psir, double *psii,
			int sigOutLength, double ys)
{
	int count;
	double con;


	for(count=0;count<sigInLength;count++)
	{
		con=2*PI*x[count]*fc/fb;
		//econ = exp(-x2/fb);
		psir[count] = (1.0-con*con)/((1.0+con*con)*(1.0+con*con)*ys);
		psii[count] = -2.0*con/((1.0+con*con)*(1.0+con*con)*ys);
	}
	return;
}

void cauchy_neo(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys)
{
    int count;
	double x2, con, econ;

	for (count=0;count<sigInLength;count++)
	{
		x2 = x[count]*x[count];
		con = 1-x2;
		econ = (1-x2)*(1-x2)+4*x2;
		psir[count]=con/(2*PI*econ*ys);
		psii[count]=2*x2/(2*PI*econ*ys);
	}
	return;
}

void cauchy_packet(double *x, int sigInLength,
			double *f, int sigOutLength, double ys)
{
	//int count;
	//double con;

	//con = 1/sqrt(PI);
	/*for(count=0;count<sigInLength;count++)
	{
		con=PI*x[count];

		f[count] = (1.0-con*con)/((1.0+con*con)*(1.0+con*con)*ys);;
		f[count+sigInLength] = -2.0*con/((1.0+con*con)*(1.0+con*con)*ys);
	}*/
    cauchy_neo(x, sigInLength, f, f+sigInLength, sigOutLength, ys);
	return;
}


/*-------------------------------------------
 * Meyer Filter Generation
 *-----------------------------------------*/


void meyer_phi(double *x, int sigInLength,
			double lb, double ub, double *phir, double *phii,
			int sigOutLength, double ys)
{
	int count;
	double con,delta,omega,xhat;
	double *xhat_r,*xhat_i;
	 xhat_r = (double *)malloc(sigInLength*sizeof(double));
	 xhat_i = (double *)malloc(sigInLength*sizeof(double));
	delta = (ub-lb)/sigInLength;

	for(count=0;count<sigInLength;count++)
	{
		xhat_r[count]=0;
		xhat_i[count]=0;
		xhat=0;
	        if (abs(x[count]) <(2*PI/3))
		  xhat=1;
		if (abs(x[count]) >=(2*PI/3) && abs(x[count]) <(4*PI/3)){
		  meyeraux(3/2/PI*abs(x[count])-1,&con);
		  xhat=cos(PI/2*con);
		}

		//tmp omega
		omega=(-sigInLength+2*count)/(2*sigInLength*delta);
		//xhat ohne fftshift
		xhat_r[count]=xhat*cos(2*PI*omega*lb)/delta;
		xhat_i[count]=xhat*sin(2*PI*omega*lb)/delta;

	}
 	fftshift(xhat_r,phir,sigInLength);
 	fftshift(xhat_i,phii,sigInLength);

 	ifft (sigInLength, sigInLength, phir, phii);
	for(count=0;count<sigInLength;count++)
	{
	  phir[count]=phir[count]*ys;
	  phii[count]=phii[count]*ys;
	}
	free(xhat_r);
	free(xhat_i);
	return;
}

void meyeraux(double x, double *y)
{
    double x4, x5, x6, x7;


	x4 = x*x*x*x;
	x5 = x4*x;
	x6 = x5*x;
	x7 = x6*x;

	*y = 35*x4-84*x5+70*x6-20*x7;
	return;
}


/*-------------------------------------------
 * CWT Utility
 *-----------------------------------------*/
void cwt_fun_parser(char *wname, int *ind)
{
    int count;

	*ind = -1;
	for(count=0;count<cwtIdentityNum;count++)
	{
		if (strcmp(wname,ci[count].wname) == 0)
		{
			*ind = count;
			break;
		}
	}
	return;
}


void full_range_scalef (char *wname, double *f, int sigOutLength)
{
    int level, ind, family, member, count, s1, s2, l;
	double one, *lowfltr, *hifltr, *trange;
	char d[2]="d";
    Func syn_fun, ana_fun;
    swt_wavelet pWaveStruct;

    level = 10;
	one = 1.0;
	wavelet_fun_parser (wname, &ind);
    if ((ind!=-1) && (wi[ind].rOrB==ORTH))
    {
      wavelet_parser(wname,&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  upcoef_len_cal (1, pWaveStruct.length, level,
	       &s1, &s2);
	  l=1;
	  //l=(int)(floor((sigOutLength-s1)/2));
	  for(count=0;count<sigOutLength;count++)
		  f[count] = 0;
	  upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, &(f[l]),
	      s1, s1, d, level);
	  if ((family==COIFLETS) || (family==SYMLETS) || (family==DMEY))
	  {
		  for(count=0;count<sigOutLength;count++)
		      f[count] = -1*f[count];
	  }
	  for(count=0;count<sigOutLength;count++)
		      f[count] = f[count]*pow(sqrt(2),level);
	  filter_clear();
	  return;
	 }

	if ((ind!=-1) && (wi[ind].rOrB==BIORTH))
    {
      wavelet_parser(wname,&family,&member);
	  ana_fun = wi[ind].analysis;
      (*ana_fun)(member, &pWaveStruct);
	  upcoef_len_cal (1, pWaveStruct.length, level,
	       &s1, &s2);
	  //l=(int)(floor((sigOutLength-s1)/2));
	  l=1;
   	  for(count=0;count<sigOutLength;count++)
	       f[count]=0;
	  lowfltr = malloc(pWaveStruct.length*sizeof(double));
	  hifltr = malloc(pWaveStruct.length*sizeof(double));
      wrev(pWaveStruct.pLowPass, pWaveStruct.length, lowfltr, pWaveStruct.length);
	  qmf_wrev(lowfltr,pWaveStruct.length,hifltr,pWaveStruct.length);
      upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, &(f[l]),
	      s1, s1, d, level);
	  free(lowfltr);
	  free(hifltr);
	  filter_clear();
      for(count=0;count<sigOutLength;count++)
	  	   f[count] = -1*f[count]*pow(sqrt(2),level);
	  return;
   }

	cwt_fun_parser(wname, &ind);
  if ((ind!=-1) && (ci[ind].realOrComplex==REAL))
  {
	  trange = malloc(sigOutLength*sizeof(double));
	  linspace(ci[ind].lb, ci[ind].ub, sigOutLength, trange, sigOutLength);
	  (*(ci[ind].scalef))(trange,sigOutLength,f,sigOutLength,ci[ind].cpsi);
	  free(trange);
	  return;
  }

  if ((ind!=-1) && (ci[ind].realOrComplex==COMPLEX))
  {
	  trange = malloc(sigOutLength*sizeof(double)/2);
	  linspace(ci[ind].lb, ci[ind].ub, sigOutLength/2, trange, sigOutLength/2);
	  (*(ci[ind].scalef))(trange,sigOutLength/2,f,sigOutLength,ci[ind].cpsi);
	  free(trange);
	  return;
  }

	return;
}

void cwt_len_cal (int sigInLength, int scale, int *sigOutLength, double *delta)
{
    *sigOutLength = scale;
	*delta = (double)(sigInLength-1)/(double)(scale-1);
	return;
}

void scale_real (double *f, int sigInLength, double delta, double *fout, int sigOutLength)
{
    int count;
	for(count=0;count<sigOutLength;count++)
	{
		fout[count] = f[(int)(floor(count*delta))];
		//sciprint("%d\n",count*delta+1);
	}
	return;
}

/*void scale_complex (double *f, int sigInLength, int delta, double *fout, int sigOutLength)
{

	return;
}*/

void cwt_conv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength)
{
	int len;
    double *fTemp, *buf;

	len = sigInLength+filterLen-1;
	buf = malloc(len*sizeof(double));
    fTemp = malloc(filterLen*sizeof(double));

	wrev(f, filterLen, fTemp, filterLen);
	conv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	free(fTemp);
    //for (i=1;i<len;i++)
	//	buf[i]=buf[i]-buf[i-1];
	wkeep_1D_center(buf,len,sigOut,sigOutLength);
	free(buf);
	return;
}

void cwt_iconv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength)
{
	int len;
    double *fTemp, *buf;

	len = sigInLength+filterLen-1;
	buf = malloc(len*sizeof(double));
    fTemp = malloc(filterLen*sizeof(double));

	wrev(f, filterLen, fTemp, filterLen);
	//iconv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	i_conv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	free(fTemp);
	wkeep_1D_center(buf,len,sigOut,sigOutLength);
	free(buf);
	return;
}

void cwt_conv_complex (double *sigIn, int sigInLength, double *fr, double *fi, int filterLen,
					   double *sigOutR, double *sigOutI, int sigOutLength)
{
    cwt_conv_real(sigIn,sigInLength,fr,filterLen,sigOutR,sigOutLength);
	cwt_conv_real(sigIn,sigInLength,fi,filterLen,sigOutI,sigOutLength);
	return;
}

void cwt_conv_complex_complex (double *a, double *b, int sigInLength,double *c, double *d,
							   int filterLen, double *sigOutR, double *sigOutI, int sigOutLength)
{
	int count;
    double *ac, *bd, *bc, *ad;
	ac = malloc(sigOutLength*sizeof(double));
	bd = malloc(sigOutLength*sizeof(double));
	bc = malloc(sigOutLength*sizeof(double));
	ad = malloc(sigOutLength*sizeof(double));

	cwt_conv_real(a,sigInLength,c,filterLen,ac,sigOutLength);
	cwt_conv_real(b,sigInLength,d,filterLen,bd,sigOutLength);
	cwt_conv_real(b,sigInLength,c,filterLen,bc,sigOutLength);
	cwt_conv_real(a,sigInLength,d,filterLen,ad,sigOutLength);

    for(count=0;count<sigOutLength;count++)
	{
		sigOutR[count]=ac[count]-bd[count];
		sigOutI[count]=bc[count]+ad[count];
	}
	free(ac);
	free(bd);
	free(bc);
	free(ad);

	return;
}


// void
// cwt_upcoef_len_cal (int sigInLength, int filterLen, int stride,
// 		int *sigOutLength, int *sigOutLengthDefault)
// {
//   int count;
//   *sigOutLength = sigInLength;
//   *sigOutLengthDefault = sigInLength;
//       for(count=0;count<stride;count++)
//       {
// 	// original version
// 	*sigOutLengthDefault = 2*(*sigOutLengthDefault) + filterLen - 1;
// 	*sigOutLength = 2*(*sigOutLength) + filterLen - 2;
//
//       }
//   return;
// }

// void
// cwt_upcoef (double *sigIn, int sigInLength, double *lowRe,double *hiRe,
// 	int filterLen, double *sigOut, int sigOutLength,
// 	int defaultLength, char *coefType, int step)
// {
//   int count, sigInLengthTemp, leng;
//   double *sigInTemp, *sigOutTemp;
//
//   // works with wavefun, cwt
//   sigInLengthTemp = 2 * sigInLength + filterLen - 2;
//
//
//
//   //sigInLengthTemp = 2 * sigInLength + filterLen - 1;
//   sigInTemp = (double *) malloc(defaultLength*sizeof(double));
//
//   if (strcmp(coefType,"a")==0)
//   {
// 	  //sciprint("recognized\n");
// // 	  printf("sigInLength %d, filterLen%d, sigInLengthTemp %d\n",sigInLength,filterLen,sigInLengthTemp);
// 	  idwt_approx_neo (sigIn, sigInLength, lowRe, filterLen,
// 		 sigInTemp, sigInLengthTemp);
// // 	  sciprint("recognized\n");
//   }
//   else
//     idwt_detail_neo (sigIn, sigInLength, hiRe, filterLen,
// 		 sigInTemp, sigInLengthTemp);
//
//   if (step > 1)
//     {
//       sigOutTemp = (double *) malloc(defaultLength*sizeof(double));
//       for(count=0;count<defaultLength;count++)
// 	sigOutTemp[count] = 0;
//       leng = sigInLengthTemp;
//       for(count=0;count<(step-1);count++) //for cwt
// 	{
// 	  //printf("leng %d, filterLen%d, leng*2-filterLen+2 %d\n",leng,filterLen,leng*2-filterLen+2);
// 	  // original version
// 	  idwt_approx_neo (sigInTemp, leng, lowRe, filterLen,
// 	               sigOutTemp, leng*2+filterLen-2);
// 	  leng = leng*2+filterLen-2;
//
// 	  verbatim_copy (sigOutTemp, leng, sigInTemp, leng);
// 	}
//       sigInLengthTemp = leng;
//       free(sigOutTemp);
//     }
//
//
//   wkeep_1D_center (sigInTemp, sigInLengthTemp, sigOut, sigOutLength);
//   free(sigInTemp);
//   return;
// }


/*====================================================================*
 * Name of the function : fftshift                                    *
 * Date of creation     : 02 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * swaps the first and second halves of a vector. Example             *
 * [1 2 3 4 5 6 7 8 9 10 ] becomes [6 7 8 9 10 1 2 3 4 5]             *
 * The parameters to pass are :                          	      *
 *   - the input vector		                            	      *
 *   - the output vector	                            	      *
 *   - its length                                                     *
 * if the length is odd, example [1 2 3 4 5] becomes [4 5 1 2 3]      *
 *====================================================================*/
int
fftshift (double *vector_in, double *vector_out, int vector_length)

{
  double inter1, inter2;
  int i, half_length;


  /* computation of the half length in case of odd or even length */
  half_length = (int) (vector_length/2.0);


  /* case where the length is odd */
  if (ISODD(vector_length)==1)
    {
      inter2=vector_in[half_length];
      for (i=0; i<half_length; i++)
	{
	  inter1 = vector_in[i];
	  vector_out[i] = vector_in[half_length+i+1];
	  vector_out[half_length + i ] = inter1;
	}
      vector_out[vector_length-1]=inter2;
    }
  /* case where the length is even */
  else
    {
      for (i=0; i<half_length; i++)
	{
	  inter1 = vector_in[half_length + i ];
	  vector_out[half_length + i] = vector_in[i];
	  vector_out[i] = inter1;
	}
    }
  /* fftshifting of the vector */


return 0;
}





int
ifft (int Signal_Length, int Nfft, double *sig_real, double *sig_imag)
{
      kiss_fft_cpx * buf;
    kiss_fft_cpx * bufout;
  kiss_fft_cfg cfg = kiss_fft_alloc(Signal_Length , 1, 0, 0);
//      fftw_complex * in=NULL;
//     fftw_complex * out=NULL;
//     fftw_plan p;

        int            i, j, k, n, n2;
      double         c, s, e, a, t1, t2;
  /*------------------------------------------------------------------*/
  /*          when the signal length is a power of two                */
  /*------------------------------------------------------------------*/
  if (Signal_Length == (int) powof (2, Nfft) + 1)
    {


      j = 0;			/* bit-reverse  */
      n2 = Signal_Length / 2;
      for (i = 1; i < Signal_Length - 1; i++)
	{
	  n = n2;
	  while (j >= n)
	    {
	      j = j - n;
	      n = n / 2;
	    }
	  j = j + n;

	  if (i < j)
	    {
	      t1 = sig_real[i];
	      sig_real[i] = sig_real[j];
	      sig_real[j] = t1;
	      t1 = sig_imag[i];
	      sig_imag[i] = sig_imag[j];
	      sig_imag[j] = t1;
	    }
	}


      n = 0;			/*IFFT  */
      n2 = 1;

      for (i = 0; i < Nfft; i++)
	{
	  n = n2;
	  n2 = n2 + n2;
	  e = 6.283185307179586 / n2;
	  a = 0.0;

	  for (j = 0; j < n; j++)
	    {
	      c = cos (a);
	      s = sin (a);
	      a = a + e;

	      for (k = j; k < Signal_Length; k = k + n2)
		{
		  t1 = c * sig_real[k + n] - s * sig_imag[k + n];
		  t2 = s * sig_real[k + n] + c * sig_imag[k + n];
		  sig_real[k + n] = sig_real[k] - t1;
		  sig_imag[k + n] = sig_imag[k] - t2;
		  sig_real[k] = sig_real[k] + t1;
		  sig_imag[k] = sig_imag[k] + t2;
		}
	    }
	}
      /* divide by Signal_Length */
      for (k = 0; k < Signal_Length; k++)
	{
	  sig_real[k] = sig_real[k] / Signal_Length;
	  sig_imag[k] = sig_imag[k] / Signal_Length;
	}
  	  free(cfg);
    }
  /*------------------------------------------------------------------*/
  /*        when the signal length is NOT a power of two              */
  /*            Calls the matlab subroutine ifft                      */
  /*------------------------------------------------------------------*/
  else
    {
//       cfg = kiss_fft_alloc(Signal_Length , 1, 0, 0);
//       mxArray       *outputArray[1];
//       mxArray       *inputArray[1];
//       mxArray       *array_ptr;


//       num_out = 1;
//       num_in = 1;

      /* recopy the real and imag parts of the signal in matrices */
//       array_ptr = mxCreateDoubleMatrix (1, Signal_Length, mxCOMPLEX);
//       memcpy (mxGetPr (array_ptr), sig_real, Signal_Length * sizeof (double));
//       memcpy (mxGetPi (array_ptr), sig_imag, Signal_Length * sizeof (double));
//       inputArray[0] = array_ptr;
       buf	= (kiss_fft_cpx*)KISS_FFT_MALLOC(Signal_Length * sizeof(kiss_fft_cpx));
       bufout	= (kiss_fft_cpx*)KISS_FFT_MALLOC(Signal_Length * sizeof(kiss_fft_cpx));
      for (i=0;i<Signal_Length;i++){
	buf[i].r=sig_real[i];
	buf[i].i=sig_imag[i];
      }
//       printf("sig in %f\n",buf[0].r);
  /*in=fftw_malloc(sizeof(fftw_complex) * Signal_Length);
    out=fftw_malloc(sizeof(fftw_complex) * Signal_Length);
    for (i=0;i<Signal_Length;++i ) {
        in[i][0] = sig_real[i];
        in[i][1] = sig_imag[i];
    }
      p = fftw_plan_dft_1d(Signal_Length, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);*/
//        for (i=0;i<Nfft;i++)
//          fftw_execute(p);
//
/*  sig_copy_real= (double*)malloc(sizeof(double) * 1*Signal_Length);
  sig_copy_imag= (double*)malloc(sizeof(double) * 1*Signal_Length);

       memcpy (sig_copy_real, sig_real, Signal_Length * sizeof (double));
      memcpy ( sig_copy_imag, sig_imag, Signal_Length * sizeof (double));*/

//        for (i=0;i<Signal_Length;i++){
// 	sig_real[i]=(double)out[i][0];
// 	sig_imag[i]=(double)out[i][1];
//       }
      /* calls the MATLAB function */
      //mexCallMATLAB (num_out, outputArray, num_in, inputArray, "ifft");
//        for (j = 0; j < Nfft; j++)
       kiss_fft(cfg, buf, bufout);
//
      for (i=0;i<Signal_Length;i++){
	sig_real[i]=(double)bufout[i].r;
	sig_imag[i]=(double)bufout[i].i;
      }
//        printf("sig out %f\n",bufout[0].r);
      /* recovers the output */

//       if (mxIsComplex (outputArray[0]))
// 	{
// 	   for (i=0;i<Signal_Length;i++){
// 	sig_copy_real[i]=(double)buf[i].r;
// 	sig_copy_imag[i]=(double)buf[i].i;
//       }
//  	  memcpy (sig_imag, sig_copy_real, Signal_Length * sizeof (double));
//       memcpy (sig_real, mxGetPr (outputArray[0]), Signal_Length * sizeof (double));
//       if (mxIsComplex (outputArray[0]))
// 	{
//    memcpy (sig_imag, sig_copy_imag, Signal_Length * sizeof (double));
// 	  memcpy (sig_imag, mxGetPi (outputArray[0]), Signal_Length * sizeof (double));
// 	}
//       else
// 	{
// 	  for (i = 0; i < Signal_Length; i++)
// 	    sig_imag[i] = 0;
// 	}

      /* free memory */
       free(cfg);free(buf);free(bufout);
//         free(sig_copy_real);free(sig_copy_imag);
//          fftw_destroy_plan(p);
//           fftw_free(in); fftw_free(out);

      //mxDestroyArray (outputArray[0]);
      //mxDestroyArray (inputArray[0]);
    }
  return 0;
}


/*void real_scale (double lb, double ub, double scale, int length, double *f, double ys, RWScaleFunc w)
{
    int ns, count;
	double *x, *ft;

    ns = (int)(floor(ub*scale)-ceil(lb*scale));
	x = malloc(ns*sizeof(double));
	ft = malloc(ns*sizeof(double));
    for(count=0;count<ns;count++)
		x[count]=ceil(lb*scale)+count;
	(*w)(x, ns, ft, ns, ys);
	free(x);
    wextend_1D_center (ft, ns, f, length, ZPD);
	free(ft);
	return;
}*/


//void complex_scale(int lb, int ub, double scale, int length, double *fr, double *fi, double ys, WScaleFunc w)
//{
//	int ns, count;
//	double *x, *ftr, *fti;
//
  //  ns = (int)(floor(ub*scale)-ceil(lb*scale));
//	x = malloc(ns*sizeof(double));
//	ftr = malloc(ns*sizeof(double));
//	fti = malloc(ns*sizeof(double));
//    for(count=0;count<ns;count++)
//		x[count]=ceil(lb*scale)+count;
//	(*w)(x, ns, ft, ns, ys);
//	free(x);
//    wextend_1D_center (ft, ns, fr, length, ZPD);
//	free(ft);
//	return;
//}


//void cwt (double *sigIn, int sigInLength,
//		  double *matrixOut, int matrixOutRow,
//		  int matrixCol, char *wname, int *scale)
//{
//	int count;
//    double *scaleFunc;

//	cwt_fun_parser(wname,&ind);

//	scaleFunc = malloc(sigInLength*sizeof(double));
//	if (ci[count].realOrComplex == 0)
//	{
//	    for(count=0:count<matrixOutRow;count++)
//	    {
//		    ci[count].scalef(ci[ind].lb, ci[ind],ub, scale[count], sigInLength, scaleFunc,sqrt(scale[count]));
//		    iconv(sigIn,sigInLength,matrixOut+count*sigInLength,sigInLength,scaleFunc,sigInLength);
//      }
//	}
//	else
//	{
//
//	}
//    free(scaleFunc);
//
//	return;
//}
