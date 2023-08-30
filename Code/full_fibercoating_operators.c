#include "full_fibercoating_operators.h"
#include<math.h>
#include<stdio.h>
#include <stdlib.h>

/********************** functions on disjoining pressure and mobility **************/

double pps_mob(double um, double u0, double alpha, double lambda,int N2)
{
	//use simpson's method for the numerical integration
	//mobility using M(u) = u^3*phi(alpha*u)/(3*phi(alpha))

	double t1, t2, t3, t4, t5, t9, t10, t11, t12,t13, t14, t22, t16,C,integ;
	if (fabs(um-u0)<TINY2)
	{
		// return m(u0), lambda = 0 assumed
		t2 = alpha*u0+1.0;
		t3 = t2*t2;
		t4 = t3*t3;
		t5 = log(t2);
		t11 = 1.0+alpha;
		t12 = t11*t11;
		t13 = t12*t12;
		t14 = log(t11);
		t22 = (t4*(4.0*t5-3.0)+4.0*t3-1.0)/(t13*(4.0*t14-3.0)+4.0*t12-1.0)/3.0;
		return t22;
	}
	else
	{
		//use simpson's method
		integ = simpson(um,u0,alpha,lambda,N2);
		return (u0-um)/integ;
	}
}

double Del_mob_operator(double um,double u0,double alpha, double lambda, int param, int N2)
{
	//derivative of m(um,u0)
	double t1,t2,t3,t4,t5,t16,t19,t18,t27,t9,t6,t7,t10,t17;
	if (fabs(um-u0)<TINY2)
	{
		//m'(u0)
		t2 = alpha*u0+1.0;
		t3 = t2*t2;
		t4 = t3*t2;
		t5 = log(t2);
		t16 = 1.0+alpha;
		t17 = t16*t16;
		t18 = t17*t17;
		t19 = log(t16);
		t27 = (4.0*t4*(4.0*t5-3.0)*alpha+4.0*t4*alpha+8.0*t2*alpha)/(t18*(4.0*t19
	-3.0)+4.0*t17-1.0)/3.0;
		return t27;
	}
	else
	{
		t1 = simpson(um,u0,alpha,lambda,N2);//int(Gdprime(s),s,um,u0);
		t4 = t1*t1;
		
		if (param==1)
		{
			
			t7 = pps_Gdprime(um, alpha, lambda);//Gdprime(um);
			t9 = -1/t1+(u0-um)/t4*t7;
			//printf("simpson = %g t4 = %g t7 = %g t9 = %g\n", t1,t3,t7,t9);
			return t9;
		}
		else if(param ==2)
		{
			//t1 = int(Gdprime(s),s,um,u0);
			//t4 = t1*t1;
			t7 = pps_Gdprime(u0, alpha, lambda);//Gdprime(u0);
			t9 = 1/t1-(u0-um)/t4*t7;
			return t9;

		}
		else
		{
			printf("param should be either 1 or 2");
			exit(1);
			return 0;
		}
	}
}

double simpson(double um, double u0, double alpha, double lambda, int N2)
{
	//calculate int_um^u0 Gdprime(s)ds using simpsons rule with N2 grids
	double a,b,fa,fb,ab,fab,integ, partial_int;
	double dx2 = (u0-um)/N2;
	integ = 0.0;
	for (int i = 0; i < N2; i++)
	{
		// for each grid int_a^b f(x) = (b-a)/6 (f(a)+4f(ab)+f(b))
		a = um+(i*dx2);
		b = um+((i+1)*dx2);
		ab = 0.5*(a+b);
		fa = pps_Gdprime(a, alpha, lambda);
		fb = pps_Gdprime(b, alpha, lambda);
		fab = pps_Gdprime(ab,alpha, lambda);
		partial_int = dx2*(fa+fb+4*fab)/6.0;
		integ+= partial_int;
		//printf("a = %12f ab = %12f b= %12f fa = %12f fab = %12f fb =%12f partial_int = %12f integ = %12f\n",a,ab,b,fa,fab,fb,partial_int,integ);
	}
	//printf("integ = %12f\n",integ);
	return integ;
}





double pps_Gdprime(double u0, double alpha, double lambda)
{
	// calculate 1/M(u0)
	double t2, t3,t4, t5,t12,t13,t14,t15,t22;

	t2 = alpha*u0+1.0;
	t3 = t2*t2;
	t4 = t3*t3;
	t5 = log(t2);
	t12 = 1.0+alpha;
	t13 = t12*t12;
	t14 = t13*t13;
	t15 = log(t12);
	t22 = 3.0/(t4*(4.0*t5-3.0)+4.0*t3-1.0)*(t14*(4.0*t15-3.0)+4.0*t13-1.0);
	return t22;
}

double Zpdiff(double um,double u0, double Hamaker,double dx)
{
	double t0;
	//Zpdiff(um,u0) = (Zp(u0)- Zp(um))/dx = (-Hamaker/u0^3 +Hamaker/um^3)/dx
	if(fabs(um-u0)<=TINY2)
	{
		return 0.0;
	}
	else
	{
	  //printf("um = %g, u0 = %g, A = %g, dx = %g\n", um,u0,Hamaker,dx);
	  t0 = (-Hamaker/(u0*u0*u0)+Hamaker/(um*um*um))/dx;
		return t0;
	}
}

double Zmdiff(double um, double u0, double eta, double alpha, double dx)
{
	double t0;
	//Zmdiff(um,u0) = (alpha/eta/(1+alpha*u0) -alpha/eta/(1+alpha*u,))/dx
	// HJI: The original one eta = 0.001 was way too small
	if(fabs(um-u0)<TINY2)
	{
		return 0.0;
	}
	else
	{
		t0 = (alpha/eta/(alpha*u0+1.0)-alpha/eta/(alpha*um+1.0))/dx;
		return t0;
	}

}


double Del_Zpdiff(double um,double u0, double Hamaker,double dx, int param)
{
	double t1,t2,t3,t4,t7,t8,t14,t9;
	//Zpdiff(um,u0) = (-Hamaker/u0^3+Hamaker/um^3)/dx
	//d1zpdiff = -3*Hamaker*um^-4/dx
	//d2zpdiff = 3*Hamaker*u0^-4/dx
	if (fabs(um-u0)<TINY2)
	{
		t3 = u0*u0;
		t4 = t3*t3;
		t9 = -Hamaker*3.0/t4/dx;
		if(param ==1)
		{
			return t9;
		}
		else if(param ==2)
		{
			return -t9;
		}
		else
		{
			printf("param should be either 1 or 2");
			exit(1);
			return 0;
		}
	}
	else
	{
		if (param==1)
		{
			t1 = um*um;
      		t2 = t1*t1;
      		t8 = -3.0*Hamaker/t2/dx;
			return t8;
		}
		else if(param ==2)
		{
			t1 = u0*u0;
      		t2 = t1*t1;
      		t7 = 3.0*Hamaker/t2/dx;
			return t7;
		}
		else
		{
			printf("param should be either 1 or 2");
			exit(1);
			return 0;
		}
	}
}

double Del_Zmdiff(double um,double u0, double eta, double alpha, double dx, int param)
{
	return 0.0;
}