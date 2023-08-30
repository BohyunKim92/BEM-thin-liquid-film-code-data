// fixed time step method for BEM and GM//
//parameteres for section 5.3.3 in the reference paper were used.
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
FILE *cpu_time_gen;
FILE *cpu_time_pps;
FILE *negative_data;
FILE *mintrack;

// FILES for unit testing
FILE *rhstest;

clock_t start_gen, end_gen;
clock_t start_pps, end_pps;
double cpu_gen,cpu_pps;
double tot_time_calc_gen,tot_time_calc_pps;
double mass;
int BAD_gen;
int memory_setup();
int define_system_parameters();
int run_gen();
int run_pps();
int build_rhs_jac_gen(double *U_gen, double *U0_gen);
int build_rhs_jac_pps(double *U_pps, double *U0_pps);
int timestep_gen(double tt_gen);
int timestep_pps(double tt_pps);
int open_files();
int record_minmax_values_gen();
double umin_gen,umax_gen;
int imin_gen, imax_gen;
int recording_minimax;
int recording_out;



int define_system_parameters()
{
    int i;
	double xx;

	//parameters for table 1
	
	param_alpha = 5.0;
    param_eta = 0.005;
    param_A = 0;
    param_lambda= 0.0;
	N = 100;
	L = 1.0; 
	dx=L/(N+1);
    tt_gen=0.0;
    tt_pps=0.0;	
    dt0 = 1e-3;
    dt_gen = dt0;
    dt_pps = dt0;
 	dtnext_gen=dt0;
    dtnext_pps=dt0;
    tnext_gen = dtnext_gen;
	tnext_pps = dtnext_pps;
	tend= 3.4765;
    err_tol = 1e-9;
	N2 = 2;
	adaptive_time = 0;
	recording_out = 1;
	recording_minimax = 1;
	fprintf(params, "%g\t %g\t %g\t %g\t %d\t   %g\t %g\t %g\t %g\t %d\n",param_alpha,param_eta,param_A,param_lambda,N,L,dx,dt0,err_tol,N2);
	

	memory_setup();
    
    for (i=0;i<N+1;i++)
 	{	
		//initial data for computational efficiency simulation
		hN = 0.45; 
 		xx = i*dx; 		
  		U0[i] =hN*(1+0.01*sin(xx*pi/L));
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
        U0_gen[i] = U0[i];
        U0_pps[i] = U0[i];
        U_gen[i] = U0_gen[i];
        U_pps[i] = U0_pps[i];
	
    }
	
	if(recording_out)
    { 
		for(int i=0;i<N+1;++i)
		{			
			fprintf(out,"%g\t %g\t %g\n",
			i*dx,U_gen[i],tt_gen);
		}
		fprintf(out,"\n\n");
		for(int i=0;i<N+1;++i)
		{			
			fprintf(myout,"%g\t %g\t %g\n",
			i*dx,U_pps[i],tt_pps);
		}
		fprintf(myout,"\n\n");
	}
   
    
    record_minmax_values_gen();
	run_gen();
	//run_pps();// pps can be ran as well with fixed time step.
    return 0;
}


int open_files()
{
	params=fopen("sim_data/parameters_for_fixed_time.dat","w");
 	out=fopen("sim_data/GM_data_for_fixed_time.dat","w");
    myout=fopen("sim_data/BEM_data_for_fixed_time.dat","w");
	cpu_time_gen=fopen("sim_data/cpu_time_GM_data.dat","w");
	cpu_time_pps=fopen("sim_data/cpu_time_BEM_data.dat","w");
	newton_test_pps=fopen("sim_data/newton_test_data_BEM.dat","w");
	newton_test_gen=fopen("sim_data/newton_test_data_GM.dat","w");
	negative_data = fopen("sim_data/negative_value_data_GM.dat","w");
	mintrack=fopen("sim_data/BEM_min_max_data.dat","w");
    

	return 0;
}

int record_minmax_values_gen()
{
    umin_gen = INFINITY;
    umax_gen = -INFINITY;
    imin_gen = 0;
    imax_gen = 0;
    int i;
    for(i=0;i<N+1;++i)
    {			
        if(U0_gen[i]>umax_gen)
        {
            imax_gen=i;
            umax_gen=U0_gen[i];
        }
    }
    for(i=1;i<N+1;++i)
    {
        if(U0_gen[i]<umin_gen)
        {
            imin_gen=i; 
            umin_gen=U0_gen[i];
        }
    }

    if (umin_gen <0)
    {
        fprintf(negative_data,"negative values of u found!\n");
        fprintf(negative_data,"printing all the parameters ... \n");
        fprintf(negative_data,"alpha=%g\t eta=%g\t A=%g\t lambda=%g\t N=%d\t L=%g\t mass=%g\n",param_alpha,param_eta,param_A,param_lambda,N,L, mass);
        fprintf(negative_data,"occuring at tt = %g, umin_gen = %g, imin_gen = %d\n\n", tt_gen, umin_gen,imin_gen);
    }	
    if (recording_minimax)
    {
        fprintf(mintrack,"%g\t %g\t %d\t %g\t %d\t %g\n", tt_gen,dt_gen, imin_gen,umin_gen, imax_gen,umax_gen);
    }
    
    return 0;
}


int run_pps()
{

    int i,b_count_pps=0;
	int test;
	int g_count_pps = 0;
	int reg_count = 0;
	int tot_count_pps = 0;
	int flag = 0;
	BAD = 0;
    
	while(tt_pps<tend)
	{
		start_pps = clock();
		if(timestep_pps(tt_pps))
		{
			flag+=1;
			g_count_pps+=1;
			tot_count_pps+=1;
			// Newton successful U_pps is the new time step U0_pps is the old time step
			end_pps = clock();
            cpu_pps = ((double) (end_pps - start_pps)) / CLOCKS_PER_SEC;
			tot_time_calc_pps +=cpu_pps;
			tt_pps+=dt_pps; // increment time
			fprintf(cpu_time_pps,"%g\t %g\t %g\t %g\n",tt_pps,dt_pps,cpu_pps,tot_time_calc_pps);
			BAD=0;

			if(recording_out)
			{ 
				// recording data of BEM method
				for(int i=0;i<N+1;++i)
				{			
					fprintf(myout,"%g\t %g\t %g\n",
					i*dx,U_pps[i],tt_pps);
				}
				fprintf(myout,"\n\n");
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
			tot_count_pps+=1;
			printf("BAD=%d\t changed dt=%g\t ttpps = %g\t b_count = %d\n",BAD,dt_pps,tt_pps, b_count_pps);
			if(b_count_pps>4) 
			{
				printf("stop here\n");
				printf("tt_pps = %g times newton succeded = %d bcount = %d, tot_count_pps = %d\n",tt_pps, g_count_pps,b_count_pps, tot_count_pps);
				exit(1);
			}
		
			continue; 
		}
	}
	printf("tt_pps = %g times newton succeded = %d bcount = %d,tot_count = %d\n",tt_pps, g_count_pps,b_count_pps, tot_count_pps);
	
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
		printf("maxerr of pps scheme=%g\n",maxerr_pps);
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
		t1 = mobumu0;
		t2 = dx*dx;
		t3 = t2*t2;
		t5 = t1/t3;

		J_pps[i][0] = t5;

		//djdum
		t1 = mobu0up;
		t2 = dx*dx;
		t4 = 1/t2/dx;
		t6 = d1mobumu0;//diff(m(um,u0),um);
		t11 = zpdumu0;//Zpdiff(um,u0);
		t12 = zmdumu0;//Zmdiff(umold,u0old);
		t15 = mobumu0;//m(um,u0);
		t17 = d1zpumu0;//diff(Zpdiff(um,u0),um);
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

int run_gen()
{

    int i,b_count_gen=0;
	int test;
	int g_count_gen = 0;
	int reg_count = 0;
	int tot_count_gen = 0;
	int flag = 0;
	BAD_gen = 0;
    while(tt_gen<tend)
	{
		start_gen = clock();
		if(timestep_gen(tt_pps))
		{
            record_minmax_values_gen();
			flag+=1;
			g_count_gen+=1;
			tot_count_gen+=1;
			// Newton successful U_pps is the new time step U0_pps is the old time step
			end_gen = clock();
            cpu_gen = ((double) (end_gen - start_gen)) / CLOCKS_PER_SEC;
			tot_time_calc_gen +=cpu_gen;
			tt_gen+=dt_gen; // increment time
			fprintf(cpu_time_gen,"%g\t %g\t %g\t %g\n",tt_gen,dt_gen,cpu_gen,tot_time_calc_gen);
			BAD_gen=0;

			if(recording_out)
			{ 
				// recording data of GM method
				for(int i=0;i<N+1;++i)
				{			
					fprintf(out,"%g\t %g\t %g\n",
					i*dx,U_gen[i],tt_gen);
				}
				fprintf(out,"\n\n");
			}
			
			for(i=0;i<N+1;++i)
			{
				U0_gen[i]=U_gen[i];  // update timestep at u^{k}
			}
            
			
		}
		else
		{
			dt_gen*=0.5;
			BAD_gen=1;
			b_count_gen++;
			tot_count_gen+=1;
			printf("BAD_gen=%d\t changed dt_gen=%g\t ttgen = %g\t b_count_gen = %d\n",BAD_gen,dt_gen,tt_gen, b_count_gen);
			if(b_count_gen>4) 
			{
				printf("Sim failed. stop here\n");
				printf("tt_gen = %g times newton succeded = %d bcount_gen = %d, tot_count_gen = %d\n",tt_gen, g_count_gen,b_count_gen, tot_count_gen);
				exit(1);
			}
		
			continue; 
		}
	}
	printf("SUCCESS! tt_gen = %g times newton succeded = %d bcount_gen = %d, tot_count_gen = %d\n",tt_gen, g_count_gen,b_count_gen, tot_count_gen);
				
	
	return(0);

}

int timestep_gen(double tt_gen) /* return 1=good, 0=bad */
{
    
    int i,j;
	double maxerr_gen;

	for(i=0;i<N+1;++i)
    {
        U_gen[i]=U0_gen[i];
    }

	for(j=0;j<15;++j)
	{
		maxerr_gen=0.0;

		build_rhs_jac_gen(U_gen,U0_gen);

		border2bandsolve(A_gen, b1_gen,b2_gen, c1_gen,c2_gen, d_gen, N-1, 2,2, A1_gen,indx_gen,f_gen, g_gen);
		
		for(i=0;i<N-1;++i)
		{
			maxerr_gen=max(maxerr_gen,fabs(f_gen[i]));
			U_gen[i]+=f_gen[i];
		}
		maxerr_gen=max(maxerr_gen,fabs(g_gen[0]));
		maxerr_gen=max(maxerr_gen,fabs(g_gen[1]));
		U_gen[N-1]+= g_gen[0];
		U_gen[N]+= g_gen[1];
		
		if(maxerr_gen<err_tol/10.0)
		{
			fprintf(newton_test_gen, "%g\t %d\t %g\n", tt_gen,j,maxerr_gen);
			break;
		}
	
	}
	
	
	if(maxerr_gen<err_tol) 
	{

		return(1);
	}
	else
	{
		printf("maxerr of generic scheme=%g\n",maxerr_gen);
		return(0);
	}
}

int build_rhs_jac_gen(double *U_gen, double *U0_gen)
{
	int i,j,ni;

	double u0,up,upp,um,umm,uold;

	#include"tempvars.h"

	for(i=0;i<N+1;++i) /* zero to start with */
	{	
		rhs_gen[i]=0.0;
		for(j=0;j<5;++j)
        {
            J_gen[i][j]=0.0;
        }
	}

    for(i=0;i<N+1;++i)
    {
    /*-------------------------*/
    /*-------------------------*/
        ni=N-i;

        u0=U_gen[i];
        up=U_gen[i+1];
        upp=U_gen[i+2];
        um=U_gen[i-1];
        umm=U_gen[i-2];
        uold=U0_gen[i];

        switch(i) /* ghost points for left symmetry BCs */
        {
                case 0:

                        um = U_gen[N];
                        umm = U_gen[N-1];

                        break;
                case 1:
                        umm = U_gen[N];

                        break;
        }
        switch(ni) /* ghost points for right symmetry BCs */
        {
                case 0:
            
                        up=U_gen[0];
                        upp=U_gen[1];
            break;
                case 1:
                        upp=U_gen[0];
            break;
        }

        /*-------------------------*/
        t2 = 1/dt_gen;
          t4 = u0*u0;
          t5 = uold*uold;
          t10 = up+u0;
          t16 = 1.0+param_alpha*t10/2.0;
          t17 = t16*t16;
          t18 = t17*t17;
          t19 = log(t16);
          t26 = t10*t10;
          t29 = param_alpha*param_alpha;
          t31 = t29/param_eta;
          t35 = 1/dx;
          t38 = up*up;
          t44 = param_A/t4/u0;
          t48 = 3.0*u0;
          t50 = dx*dx;
          t52 = 1/t50/dx;
          t56 = 1.0+param_alpha;
          t57 = t56*t56;
          t58 = t57*t57;
          t59 = log(t56);
          t65 = 1/(t58*(4.0*t59-3.0)+4.0*t57-1.0);
          t68 = u0+um;
          t74 = 1.0+param_alpha*t68/2.0;
          t75 = t74*t74;
          t76 = t75*t75;
          t77 = log(t74);
          t84 = t68*t68;
          t92 = um*um;
          t107 = (u0-uold)*t2+param_alpha*(t4-t5)*t2/2.0+8.0/3.0*(t10*t10*(t18*(4.0
    *t19-3.0)+4.0*t17-1.0)/t26*(1.0+t31/t17*(up-u0)*t35-(-param_A/t38/up+t44)*t35+(
    upp-3.0*up+t48-um)*t52)*t65/8.0-t68*t68*(t76*(4.0*t77-3.0)+4.0*t75-1.0)/t84*(
    1.0+t31/t75*(u0-um)*t35-(-t44+param_A/t92/um)*t35+(up-t48+3.0*um-umm)*t52)*t65/
    8.0)*t35;


        rhs_gen[i]= t107;
        //printf("%g\t %g\n",u0,rhs_gen[i]);

        /*-------------------------*/


        /*--------------*/
          t1 = u0+um;
          t7 = 1.0+param_alpha*t1/2.0;
          t8 = t7*t7;
          t9 = t8*t8;
          t10 = log(t7);
          t17 = t1*t1;
          t20 = dx*dx;
          t21 = t20*t20;
          t24 = 1.0+param_alpha;
          t25 = t24*t24;
          t26 = t25*t25;
          t27 = log(t24);
          t36 = t1*t1*(t9*(4.0*t10-3.0)+4.0*t8-1.0)/t17/t21/(t26*(4.0*t27-3.0)+4.0*
    t25-1.0)/3.0;

        J_gen[i][0]=t36;

        /*--------------*/
        t1 = up+u0;
          t7 = 1.0+param_alpha*t1/2.0;
          t8 = t7*t7;
          t9 = t8*t8;
          t10 = log(t7);
          t17 = t1*t1;
          t20 = dx*dx;
          t22 = 1/t20/dx;
          t24 = 1.0+param_alpha;
          t25 = t24*t24;
          t26 = t25*t25;
          t27 = log(t24);
          t33 = 1/(t26*(4.0*t27-3.0)+4.0*t25-1.0);
          t37 = u0+um;
          t38 = t37*t37/4.0;
          t41 = 1.0+param_alpha*t37/2.0;
          t42 = t41*t41;
          t43 = t42*t42;
          t44 = log(t41);
          t46 = 4.0*t44-3.0;
          t49 = t43*t46+4.0*t42-1.0;
          t51 = t37*t37;
          t53 = 1/t51/t37;
          t54 = param_alpha*param_alpha;
          t55 = 1/param_eta;
          t56 = t54*t55;
          t57 = 1/t42;
          t58 = u0-um;
          t60 = 1/dx;
          t63 = u0*u0;
          t67 = um*um;
          t77 = 1.0+t56*t57*t58*t60-(-param_A/t63/u0+param_A/t67/um)*t60+(up-3.0*u0
    +3.0*um-umm)*t22;
          t79 = t53*t77*t33;
          t83 = t38*t37/2.0;
          t84 = t42*t41;
          t96 = t83*t49;
          t97 = t51*t51;
          t111 = t67*t67;
          t123 = (-t1*t1*(t9*(4.0*t10-3.0)+4.0*t8-1.0)/t17*t22*t33/3.0-4.0*t38*t49*
    t79-8.0/3.0*t83*(2.0*t84*t46*param_alpha+4.0*t41*param_alpha+2.0*t84*
    param_alpha)*t79+8.0*t96/t97*t77*t33-8.0/3.0*t96*t53*(-t54*param_alpha*t55/t84*
    t58*t60-t56*t57*t60+3.0*param_A/t111*t60+3.0*t22)*t33)*t60;

        J_gen[i][1]=t123;

        /*--------------*/
        t1 = 1/dt_gen;
          t4 = up+u0;
          t5 = t4*t4/4.0;
          t8 = 1.0+param_alpha*t4/2.0;
          t9 = t8*t8;
          t10 = t9*t9;
          t11 = log(t8);
          t13 = 4.0*t11-3.0;
          t16 = t10*t13+4.0*t9-1.0;
          t18 = t4*t4;
          t20 = 1/t18/t4;
          t21 = param_alpha*param_alpha;
          t22 = 1/param_eta;
          t23 = t21*t22;
          t24 = 1/t9;
          t25 = up-u0;
          t27 = 1/dx;
          t30 = up*up;
          t34 = u0*u0;
          t37 = param_A/t34/u0;
          t41 = 3.0*u0;
          t43 = dx*dx;
          t45 = 1/t43/dx;
          t47 = 1.0+t23*t24*t25*t27-(-param_A/t30/up+t37)*t27+(upp-3.0*up+t41-um)*
    t45;
          t49 = 1.0+param_alpha;
          t50 = t49*t49;
          t51 = t50*t50;
          t52 = log(t49);
          t58 = 1/(t51*(4.0*t52-3.0)+4.0*t50-1.0);
          t59 = t20*t47*t58;
          t63 = t5*t4/2.0;
          t64 = t9*t8;
          t76 = t63*t16;
          t77 = t18*t18;
          t84 = t21*param_alpha*t22;
          t91 = t34*t34;
          t95 = 3.0*param_A/t91*t27;
          t96 = 3.0*t45;
          t102 = u0+um;
          t103 = t102*t102/4.0;
          t106 = 1.0+param_alpha*t102/2.0;
          t107 = t106*t106;
          t108 = t107*t107;
          t109 = log(t106);
          t111 = 4.0*t109-3.0;
          t114 = t108*t111+4.0*t107-1.0;
          t116 = t102*t102;
          t118 = 1/t116/t102;
          t119 = 1/t107;
          t120 = u0-um;
          t124 = um*um;
          t133 = 1.0+t23*t119*t120*t27-(-t37+param_A/t124/um)*t27+(up-t41+3.0*um-
    umm)*t45;
          t135 = t118*t133*t58;
          t139 = t103*t102/2.0;
          t140 = t107*t106;
          t152 = t139*t114;
          t153 = t116*t116;
          t172 = t1+t1*param_alpha*u0+(4.0*t5*t16*t59+8.0/3.0*t63*(2.0*t64*t13*
    param_alpha+2.0*t64*param_alpha+4.0*t8*param_alpha)*t59-8.0*t76/t77*t47*t58+8.0
    /3.0*t76*t20*(-t84/t64*t25*t27-t23*t24*t27+t95+t96)*t58-4.0*t103*t114*t135-8.0/
    3.0*t139*(2.0*t140*t111*param_alpha+4.0*t106*param_alpha+2.0*t140*param_alpha)*
    t135+8.0*t152/t153*t133*t58-8.0/3.0*t152*t118*(-t84/t140*t120*t27+t23*t119*t27-
    t95-t96)*t58)*t27;


        J_gen[i][2]=t172;

        /*--------------*/
         t1 = up+u0;
          t2 = t1*t1/4.0;
          t5 = 1.0+param_alpha*t1/2.0;
          t6 = t5*t5;
          t7 = t6*t6;
          t8 = log(t5);
          t10 = 4.0*t8-3.0;
          t13 = t7*t10+4.0*t6-1.0;
          t15 = t1*t1;
          t17 = 1/t15/t1;
          t18 = param_alpha*param_alpha;
          t19 = 1/param_eta;
          t20 = t18*t19;
          t21 = 1/t6;
          t22 = up-u0;
          t24 = 1/dx;
          t27 = up*up;
          t31 = u0*u0;
          t40 = dx*dx;
          t42 = 1/t40/dx;
          t44 = 1.0+t20*t21*t22*t24-(-param_A/t27/up+param_A/t31/u0)*t24+(upp-3.0*
    up+3.0*u0-um)*t42;
          t46 = 1.0+param_alpha;
          t47 = t46*t46;
          t48 = t47*t47;
          t49 = log(t46);
          t55 = 1/(t48*(4.0*t49-3.0)+4.0*t47-1.0);
          t56 = t17*t44*t55;
          t60 = t1*t2/2.0;
          t61 = t6*t5;
          t73 = t60*t13;
          t74 = t15*t15;
          t88 = t27*t27;
          t99 = u0+um;
          t105 = 1.0+param_alpha*t99/2.0;
          t106 = t105*t105;
          t107 = t106*t106;
          t108 = log(t105);
          t115 = t99*t99;
          t123 = (4.0*t2*t13*t56+8.0/3.0*t60*(2.0*t61*t10*param_alpha+4.0*t5*
    param_alpha+2.0*t61*param_alpha)*t56-8.0*t73/t74*t44*t55+8.0/3.0*t73*t17*(-t18*
    param_alpha*t19/t61*t22*t24+t20*t21*t24-3.0*param_A/t88*t24-3.0*t42)*t55-t99*
    t99*(t107*(4.0*t108-3.0)+4.0*t106-1.0)/t115*t42*t55/3.0)*t24;


        J_gen[i][3]=t123;

        /*--------------*/
         t1 = up+u0;
          t7 = 1.0+param_alpha*t1/2.0;
          t8 = t7*t7;
          t9 = t8*t8;
          t10 = log(t7);
          t17 = t1*t1;
          t20 = dx*dx;
          t21 = t20*t20;
          t24 = 1.0+param_alpha;
          t25 = t24*t24;
          t26 = t25*t25;
          t27 = log(t24);
          t36 = t1*t1*(t9*(4.0*t10-3.0)+4.0*t8-1.0)/t17/t21/(t26*(4.0*t27-3.0)+4.0*
    t25-1.0)/3.0;



        J_gen[i][4]=t36;

        /*-------------------------*/
        /*-------------------------*/
    }

    /* set-up periodic boundary conditions */

	/*-------------------------*/
	for(i=0;i<N-1;++i)
    {
        for(j=0;j<5;++j)
        {
			A_gen[i][j] = J_gen[i][j];
        }
        b1_gen[i] = 0.0;
        b2_gen[i] = 0.0;
        c1_gen[i] = 0.0;
        c2_gen[i] = 0.0;
    }

	b1_gen[0] = J_gen[0][0];
	b1_gen[N-3] = J_gen[N-3][4];
	b1_gen[N-2] = J_gen[N-2][3];
		
	b2_gen[0] = J_gen[0][1];
	b2_gen[1] = J_gen[1][0];
	b2_gen[N-2] = J_gen[N-2][4];

	c1_gen[0] = J_gen[N-1][4];
	c1_gen[N-3] = J_gen[N-1][0];
	c1_gen[N-2] = J_gen[N-1][1];

	
	c2_gen[0] = J_gen[N][3];
	c2_gen[1] = J_gen[N][4];
	c2_gen[N-2] = J_gen[N][0];

	d_gen[0][0] = J_gen[N-1][2];
	d_gen[0][1] = J_gen[N-1][3];
	d_gen[1][0] = J_gen[N][1];
	d_gen[1][1] = J_gen[N][2];

	for(i=0;i<N+1;++i)
    {
        rhs_gen[i]*= -1.0;
    }
		

	for(i=0;i<N-1;++i)
    {
        f_gen[i] = rhs_gen[i];
    } 
	g_gen[0] = rhs_gen[N-1];
	g_gen[1] = rhs_gen[N];

	return(0);

}



int memory_setup()
{
    U_gen=myalloc(N+1);U_pps=myalloc(N+1);
    U0_gen=myalloc(N+1);U0_pps=myalloc(N+1);
    U0=myalloc(N+1);        
    rhs_gen=myalloc(N+1);rhs_pps=myalloc(N+1);
    J_gen=matrix(N+1,5);J_pps=matrix(N+1,5);
    J1_gen=matrix(N+1,5); J1_pps=matrix(N+1,5);
    
    indx_gen=myalloci(N-1);indx_pps=myalloci(N-1);
    A1_gen = matrix(N-1,5);A1_pps = matrix(N-1,5);
    A_gen = matrix(N-1,5);A_pps = matrix(N-1,5);
    f_gen = myalloc(N-1);f_pps = myalloc(N-1);
    g_gen = myalloc(2);g_pps = myalloc(2);
    b1_gen = myalloc(N-1);b1_pps = myalloc(N-1);
    b2_gen = myalloc(N-1); b2_pps = myalloc(N-1);
    c1_gen = myalloc(N-1);c1_pps = myalloc(N-1);
    c2_gen = myalloc(N-1);c2_pps = myalloc(N-1);
    d_gen = matrix(2,2);d_pps = matrix(2,2);
    
    return(0);
}