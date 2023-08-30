/* Adaptive time stepping algorithm for BEM*/
/* parameters for section 5.3.3 of the reference paper were used.*/

#include<stdio.h>
#include<string.h>
#include<math.h>
#include"myalloc.h"
/*--------------------------------------------------------*/
#define sqr(x) ((x)*(x))
#define sgn(x) ((x)>=0.0 ? 1.0 : -1.0)
#define pi M_PI
#define MAXLINE 80
/*------------------------------------------------------------------------*/
#define swap(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#include "borderband.h"
#include "borderband.c"
#include "full_fibercoating_operators.h"
#include "full_fibercoating_operators.c"
/*------------------------------------------------------------------------*/
#include "sysvars.h"
#include <time.h>


FILE *out;
FILE *myout;
FILE *params;
FILE *cpu_time_pps;
FILE *negative_data;
FILE *mintrack;

// FILES for unit testing
FILE *rhstest;

/*variables used for adaptive time stepping*/
double *U_pps_prev;
double *dim_less_err;
double *dim_less_errNew;
double maxLTE;
double tol1 = 0.0;
int LTE_count;
int LTE_maxct;
int dtup_used;
double step_up_percent;

clock_t start_pps, end_pps;
double cpu_pps;
double tot_time_calc_pps;
double mass;
double dt_old_pps;

int memory_setup();
int define_system_parameters();
int run_pps();
int build_rhs_jac_pps(double *U_pps, double *U0_pps);
int timestep_pps(double tt_pps);
int open_files();
int recording_minimax;
int recording_out;



int define_system_parameters()
{
    int i;
	double xx;

	// parameters set for section 5.3.3
	
	param_alpha = 5.0;
    param_eta = 0.005;
    param_A = 0.0;
    param_lambda= 0.0;
	N = 100;
	L = 1.0; 
	dx=L/(N+1);
    tt_pps=0.0;	
    dt0 = 1e-3;
	dt_old_pps = dt0;
	dt_pps = dt0;
	dtnext_pps=dt0;
	tnext_pps = dtnext_pps;
	tend=5.0; // simulation end time
    err_tol = 1e-9;
	N2 = 2;
	adaptive_time = 1;
	recording_out = 1;
	fprintf(params, "%g\t %g\t %g\t %g\t %d\t   %g\t %g\t %g\t %g\t %d\n",param_alpha,param_eta,param_A,param_lambda,N,L,dx,dt0,err_tol,N2);
	

	memory_setup();
    
    for (i=0;i<N+1;i++)
 	{	
		//initial condition for adaptive time stepping
		hN = 0.45; 
 		xx = i*dx; 		
  		U0[i] =hN*(1+0.01*sin(xx*pi/L));//$h_0(x) = 0.45+0.01\sin(\pi x)$

		fprintf(params, "%g\t %g\n",xx,U0[i]);
	
 	}
    return 0;
}

int main(int argc,char *argv[])
{
	open_files();
	define_system_parameters();
	
    for (int i = 0 ; i <N+1; ++i)
    {
        U0_pps[i] = U0[i];
		U_pps_prev[i] = U0[i];
        U_pps[i] = U0_pps[i];
    }
	
	if(recording_out)
    { 
		// recording initial data of BEM method
		for(int i=0;i<N+1;++i)
		{			
			fprintf(myout,"%g\t %g\t %g\t %g\t %g\n",
			i*dx,U_pps[i],tt_pps,U0_pps[i],U_pps_prev[i]);
		}
		fprintf(myout,"\n\n");
	}
	run_pps();
    return 0;
}



int open_files()
{
	params=fopen("sim_data/parameters_for_adaptive_time.dat","w");
 	myout=fopen("sim_data/BEM_data_adaptive_time.dat","w");
	cpu_time_pps=fopen("sim_data/cpu_time_data_for_BEM_adaptive_time.dat","w");
	newton_test_pps=fopen("sim_data/newton_test_data_for_BEM_adaptive_time.dat","w");
	return 0;
}


int run_pps()
{

    int i,b_count_pps=0;
	int test;
	int g_count_pps = 0;
	int reg_count = 0;
	int tot_count = 0;
	int flag = 0;
	BAD = 0;
	LTE_count = 0;
	LTE_maxct = 3;
	tol1 = 1e-3;
	dtup_used = 0;
	step_up_percent = 1.2;
	while(tt_pps<tend)
	{
		start_pps = clock();
		if(timestep_pps(tt_pps))
		{
			flag+=1;
			g_count_pps+=1;
			tot_count+=1;

			// Newton successful U_pps is the new time step U0_pps is the old time step
			end_pps = clock();
            cpu_pps = ((double) (end_pps - start_pps)) / CLOCKS_PER_SEC;
			tot_time_calc_pps +=cpu_pps;
			tt_pps+=dt_pps; // increment time
			fprintf(cpu_time_pps,"%g\t %g\t %g\t %g\n",tt_pps,dt_pps,cpu_pps,tot_time_calc_pps);
			maxLTE = 0.0;
			BAD=0;

			if(recording_out)
			{ 
				// recording data of BEM method
				for(int i=0;i<N+1;++i)
				{			
					fprintf(myout,"%g\t %g\t %g\t %g\t %g\n",
					i*dx,U_pps[i],tt_pps,U0_pps[i],U_pps_prev[i]);
				}
				fprintf(myout,"\n\n");
			}

			if(adaptive_time ==1)
			{
				// calculate max of local truncation error
				for(i=0; i<N+1;++i)
				{
					dim_less_err[i] = (U0_pps[i]-U_pps_prev[i])/U_pps_prev[i];
					dim_less_errNew[i] = (U_pps[i]-U0_pps[i])/U0_pps[i];
					maxLTE = max(maxLTE, fabs(dim_less_errNew[i]-((dt_pps/dt_old_pps)*dim_less_err[i])));
				}
				
				dt_old_pps = dt_pps; 
				dt_pps*= 1.01;
				reg_count++;

				// update time step and data 
				for(i=0;i<N+1;++i)
				{
					U_pps_prev[i] = U0_pps[i]; // update time at u^{k-1}
				}
				//if(maxLTE < tol1 && maxLTE > tol1*(0.1*0.1))
				if(maxLTE < tol1)
				{
					LTE_count+=1;
					if(LTE_count == LTE_maxct)
					{
						dtup_used+=1;
						dt_pps*= step_up_percent;
						LTE_count = 0;
					}
				}
				
			
			}
			
			for(i=0;i<N+1;++i)
			{
				U0_pps[i]=U_pps[i];  // update timestep at u^{k}
			}
			
		}
		else
		{
			dt_pps*=0.5;
			BAD=1;
			b_count_pps++;
			tot_count+=1;
			if(b_count_pps>4) 
			{
				printf("Sim failed stop here\n");
				printf("tt_pps = %g adaptive time used = %d times newton succeded = %d bcount = %d, reg_count = %d tot_count = %d\n",tt_pps,dtup_used, g_count_pps,b_count_pps, reg_count, tot_count);
				exit(1);
			}
		
			continue; 
		}
	}
	printf("SUCCESS! tt_pps = %g adaptive time used = %d times newton succeded = %d bcount = %d, reg_count = %d tot_count = %d,  tol1 = %g\n",tt_pps,dtup_used, g_count_pps,b_count_pps, reg_count, tot_count,tol1);
	
	return(0);

}


int timestep_pps(double tt_pps) /* return 1=good, 0=bad */
{
    
    int i,j;
	double maxerr_pps;

	for(i=0;i<N+1;++i)
    {
        U_pps[i]=U0_pps[i];
    }

	for(j=0;j<15;++j)
	{
		maxerr_pps=0.0;

		build_rhs_jac_pps(U_pps,U0_pps);
		border2bandsolve(A_pps, b1_pps,b2_pps, c1_pps,c2_pps, d_pps, N-1, 2,2, A1_pps,indx_pps,f_pps, g_pps);
		
		for(i=0;i<N-1;++i)
		{
			maxerr_pps=max(maxerr_pps,fabs(f_pps[i]));
			U_pps[i]+=f_pps[i];
		}
		maxerr_pps=max(maxerr_pps,fabs(g_pps[0]));
		maxerr_pps=max(maxerr_pps,fabs(g_pps[1]));
		U_pps[N-1]+= g_pps[0];
		U_pps[N]+= g_pps[1];
		
        if(maxerr_pps<err_tol/10.0)
		{
			fprintf(newton_test_pps, "%g\t %d\t %g\n", tt_pps,j,maxerr_pps);
			break;
		}
			
	
	}
	
	
	if(maxerr_pps<err_tol) 
	{	
		return(1);
	}
	else
	{
		return(0);
	}
}

int build_rhs_jac_pps(double *U_pps, double *U0_pps)
{
    #include"tempvars.h"
	int i,j,ni;
    double umm,um,u0,up,upp,umold,u0old,upold;
    double mobumu0,mobu0up;
    double zpdumu0,zpdu0up, zmdumu0,zmdu0up;
    double d1zpumu0,d2zpumu0,d1zpu0up,d2zpu0up;
    double d1mobumu0,d2mobumu0,d1mobu0up,d2mobu0up;

	for(i=0;i<N+1;++i) 
	{	
		rhs_pps[i]=0.0;
		for(j=0;j<5;++j)
		{
			J_pps[i][j]=0.0;
		}
	}
	
    for(i=0;i<N+1;++i)
    {
		ni=N-i;

		u0=U_pps[i];
		up=U_pps[i+1];
		upp=U_pps[i+2];
		um=U_pps[i-1];
		umm=U_pps[i-2];
		u0old=U0_pps[i];
		upold = U0_pps[i+1];
		umold = U0_pps[i-1];

		switch(i) // ghost points for left periodic BCs
		{
			case 0:
					um = U_pps[N]; umold = U0_pps[N];
					umm = U_pps[N-1];

					break;
			case 1:
					umm = U_pps[N];

					break;
		}
		switch(ni) // ghost points for right periodic BCs 
		{
			case 0:
					up=U_pps[0]; upold = U0_pps[0];
					upp=U_pps[1];
					break;
			case 1:
					upp=U_pps[0];
					break;
		}

		//mobility definition
		mobumu0 = pps_mob(um,u0,param_alpha,param_lambda,N2); mobu0up = pps_mob(u0,up,param_alpha,param_lambda,N2);
		
		//definition on zpdiff and zmdiff//
		zpdumu0 = Zpdiff(um,u0,param_A,dx); zpdu0up = Zpdiff(u0,up,param_A,dx);
		zmdumu0 = Zmdiff(umold,u0old,param_eta,param_alpha,dx);zmdu0up = Zmdiff(u0old,upold,param_eta,param_alpha,dx);

		d1zpumu0 = Del_Zpdiff(um,u0, param_A,dx, 1);d2zpumu0 = Del_Zpdiff(um,u0, param_A,dx, 2);
		d1zpu0up = Del_Zpdiff(u0,up, param_A,dx, 1);d2zpu0up = Del_Zpdiff(u0,up, param_A,dx, 2);

		d1mobumu0 = Del_mob_operator(um,u0,param_alpha,param_lambda,1,N2);d2mobumu0 = Del_mob_operator(um,u0,param_alpha,param_lambda,2,N2);
		d1mobu0up = Del_mob_operator(u0,up,param_alpha,param_lambda,1,N2);d2mobu0up = Del_mob_operator(u0,up,param_alpha,param_lambda,2,N2);

		t9 = mobu0up;
		t11 = 3.0*u0;
		t13 = dx*dx;
		t15 = 1/t13/dx;
		t17 = zpdu0up;
		t18 = zmdu0up;
		t21 = mobumu0;
		t25 = zpdumu0;
		t26 = zmdumu0;
		t32 = (1.0+0.5*param_alpha*(u0+u0old))*(u0-u0old)/dt_pps+(t9*(0.1E1+(upp-3.0*up+t11-
		um)*t15-t17-t18)-t21*(0.1E1+(up-t11+3.0*um-umm)*t15-t25-t26))/dx;

		rhs_pps[i] = t32;
	
		// djdumm	
		t1 = mobumu0;//m(um,u0);
		t2 = dx*dx;
		t3 = t2*t2;
		t5 = t1/t3;

		J_pps[i][0] = t5;

		//djdum
		t1 = mobu0up;
		t2 = dx*dx;
		t4 = 1/t2/dx;
		t6 = d1mobumu0;
		t11 = zpdumu0;
		t12 = zmdumu0;
		t15 = mobumu0;
		t17 = d1zpumu0;
		t22 = (-t1*t4-t6*(0.1E1+(up-3.0*u0+3.0*um-umm)*t4-t11-t12)-t15*(3.0*t4-
		t17))/dx;

		J_pps[i][1] = t22;

		//dJdu0
		t3 = 1/dt_pps;
		t11 = d1mobu0up;
		t13 = 3.0*u0;
		t15 = dx*dx;
		t17 = 1/t15/dx;
		t19 = zpdu0up;
		t20 = zmdu0up;
		t23 = mobu0up;
		t24 = 3.0*t17;
		t25 = d1zpu0up;
		t28 = d2mobumu0;
		t32 = zpdumu0;
		t33 = zmdumu0;
		t36 = mobumu0;
		t37 = d2zpumu0;
		t43 = 0.5*param_alpha*(u0-u0old)*t3+(1.0+0.5*param_alpha*(u0+u0old))*t3+(t11*(0.1E1+(
		upp-3.0*up+t13-um)*t17-t19-t20)+t23*(t24-t25)-t28*(0.1E1+(up-t13+3.0*um-umm)*
		t17-t32-t33)-t36*(-t24-t37))/dx;
		J_pps[i][2] = t43;

		//djdup
		t1 = d2mobu0up;
		t5 = dx*dx;
		t7 = 1/t5/dx;
		t9 = zpdu0up;
		t10 = zmdu0up;
		t13 = mobu0up;
		t15 = d2zpu0up;
		t18 = mobumu0;
		t22 = (t1*(0.1E1+(upp-3.0*up+3.0*u0-um)*t7-t9-t10)+t13*(-3.0*t7-t15)-t18*
		t7)/dx;

		J_pps[i][3] = t22;

		//djdupp
		t1 = mobu0up;
		t2 = dx*dx;
		t3 = t2*t2;
		t5 = t1/t3;

		J_pps[i][4] = t5;

	}

	
	for(int i=0;i<N-1;++i)
	{
		for(j=0;j<5;++j)
		{
			A_pps[i][j] = J_pps[i][j];
		}

		b1_pps[i] = 0.0;
		b2_pps[i] = 0.0;
		c1_pps[i] = 0.0;
		c2_pps[i] = 0.0;
	} 

	b1_pps[0] = J_pps[0][0];//b1[0] = da0dumm
	b1_pps[N-3] = J_pps[N-3][4]; //daN-3dpp
	b1_pps[N-2] = J_pps[N-2][3]; // dan-2dup

	b2_pps[0] = J_pps[0][1];//b2[0] = da0dum
	b2_pps[1] = J_pps[1][0];//b2[1] = da1dumm
	b2_pps[N-2] = J_pps[N-2][4]; //daNdupp

	c1_pps[0] = J_pps[N-1][4];//daN-1/upp
	c1_pps[N-3] = J_pps[N-1][0];//daN-1dumm
	c1_pps[N-2] = J_pps[N-1][1];//daN-1dum

	c2_pps[0] = J_pps[N][3];//daNdup
	c2_pps[1] = J_pps[N][4];//daNdupp
	c2_pps[N-2] = J_pps[N][0];//dNdumm

	d_pps[0][0] = J_pps[N-1][2];//daN-1du0
	d_pps[0][1] = J_pps[N-1][3];//daN-1dup
	d_pps[1][0] = J_pps[N][1];//daNdum
	d_pps[1][1] = J_pps[N][2];//daNdu0

		
	for(i=0;i<N+1;++i)
	{
		rhs_pps[i]*= -1.0;
	}	

	for(i=0;i<N-1;++i)
	{
		f_pps[i] = rhs_pps[i];
	} 
		
	////last two values of RHS 
	g_pps[0] = rhs_pps[N-1];
	g_pps[1] = rhs_pps[N];
		
	return 0;
}


int memory_setup()
{
	// only allocate memory for BEM method
    U_pps=myalloc(N+1);
    U0_pps=myalloc(N+1);
    U0=myalloc(N+1);
	U_pps_prev=myalloc(N+1);
	dim_less_err=myalloc(N+1);
	dim_less_errNew =myalloc(N+1);

    rhs_pps=myalloc(N+1);
    J_pps=matrix(N+1,5);
    J1_pps=matrix(N+1,5);
    
    indx_pps=myalloci(N-1);
    A1_pps = matrix(N-1,5);
    A_pps = matrix(N-1,5);
    f_pps = myalloc(N-1);
    g_pps = myalloc(2);
    b1_pps = myalloc(N-1);
    b2_pps = myalloc(N-1);
    c1_pps = myalloc(N-1);
    c2_pps = myalloc(N-1);
    d_pps = matrix(2,2);
    
    return(0);
}