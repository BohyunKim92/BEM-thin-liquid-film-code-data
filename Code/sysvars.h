/*=============================================================*/
int N,N2,flag_plot;
double L,dx;
double err_tol;
double tend;
/*---------------------------*/

double *x;
double *U_gen; double *U_pps;
double *U0; double *U0_gen; double *U0_pps;

double *rhs_gen;double *rhs_pps;
double **J_gen;double **J_pps;
double **J1_gen;double **J1_pps;
int *indx_gen;int *indx_pps;

// for bordered band matrix with periodic bcs
double **A_gen; double **A_pps;
double **A1_gen;double **A1_pps;
double *b1_gen;double *b1_pps;
double *b2_gen;double *b2_pps;
double *c1_gen;double *c1_pps;
double *c2_gen;double *c2_pps;
double **d_gen;double **d_pps;
double *f_gen;double *f_pps;
double *g_gen;double *g_pps;

/*-------------------------*/

double param_alpha;
double param_eta; 
double param_lambda;
double param_A;
double hmax0,hmin0;
double hN; // characteristic uniform film thickness
double mass0;

double dt0,tt_gen,tt_pps,dtmax,dt_gen,dt_pps;
double tnext_gen,tnext_pps,dtnext_gen,dtnext_pps;
int BAD;


/*=====================================================*/
//variables used for testing purpose 

// FILES for l2 err testing

// FILES for unit testing
FILE *rhstest;
FILE *jactest;
FILE *mobtest;
FILE *newton_test_pps;
FILE *newton_test_gen;
FILE *Ztest;
FILE *uxxxtest;
FILE *pdifftest;
FILE *mpdifftest;
int initial_data_test;
int adaptive_time;
