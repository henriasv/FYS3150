

#include "Helium_Solver.h"
#include "omp.h"
#include <cmath>
#include "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/lib/cppLibrary/lib.h"

#define PI 3.14150265358979323


Helium_Solver::Helium_Solver(int inN) {
	alpha = 2;
	N = inN;
}

Gauss_Legendre::Gauss_Legendre(int inN, double a, double b) : Helium_Solver(inN){
	N = inN;
	x = new double[N];
	w = new double[N];
}

Gauss_Laguerre::Gauss_Laguerre(int inN) : Helium_Solver(inN) {
	alph = 2;
	N = inN;
	x_r = new double[N+1];
	w_r = new double[N+1];
	x_phi = new double[N];
	w_phi = new double[N];
	x_theta = new double[N];
	w_theta = new double[N];
	gauleg(0, 2*PI, x_phi, w_phi, N);
	gauleg(0, PI, x_theta, w_theta, N);
	gaulag(x_r, w_r, N, alph);
;}

double Gauss_Legendre::helium_function(double x1, double x2, double y1, double y2, double z1, double z2) 
{
	alpha = 2;
	double r1r2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
	double exp1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	double exp2 = sqrt(x2*x2 + y2*y2 + z2*z2);
	if (r1r2>1e-6) {
	return exp(-2*alpha*(exp1+exp2))/r1r2;
	} else return 0;
}


double Gauss_Laguerre::helium_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2) {
	alpha = 2;
}

double Gauss_Legendre::solve(int N_threads) {
	;
}

double Gauss_Laguerre::solve(int N_threads) {
	;
}

#define MAXIT 10
#define EPS 3.0e-14
void Gauss_Laguerre::gaulag(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
double Gauss_Laguerre::gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#undef MAXIT
#undef EPS
