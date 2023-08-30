double Zpdiff(double um,double u0, double Hamaker,double dx);
double Zmdiff(double um, double u0, double eta, double alpha, double dx);
double Del_Zpdiff(double um,double u0, double Hamaker,double dx, int param);
double Del_Zmdiff(double um,double u0, double eta, double alpha, double dx, int param);
double phi(double alpha);
double pps_mob(double um, double u0, double alpha, double lambda, int N2);
double pps_Gdprime(double u0, double alpha, double lambda);
double Del_mob_operator(double um,double u0,double alpha, double lambda, int param,int N2);
double simpson(double um, double u0, double alpha, double lambda, int N2);
#define TINY2 1.0e-9