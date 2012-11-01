/* 
 * File:   Gauss_Laguerre.cpp
 * Author: henriasv
 * 
 * Created on 30. oktober 2012, 21:10
 */
#define PI 3.1415926536
#include "Gauss_Laguerre.h"
#include <cmath>
#include "omp.h"

using namespace std;

Gauss_Laguerre::Gauss_Laguerre(int _N) : Helium_Solver(_N) {
	alph = 2;
	N = _N;
	x_r = new double[N+1];
	w_r = new double[N+1];
	x_phi = new double[N];
	w_phi = new double[N];
	x_theta = new double[N];
	w_theta = new double[N];
	gauleg(0, 2*PI, x_phi, w_phi, N);
	gauleg(0, PI, x_theta, w_theta, N);
	gaulag(x_r, w_r, N, alph);
        
        /*
        for (int i = 0; i<N; i++) {
            cout << x_phi[i] << " " << x_theta[i] << endl;
            double x_tmp = x_theta[i];
            x_theta[i] = cos(x_tmp);
            w_theta[i] = sin(x_tmp)*w_theta[i];//w_theta[i] *PI/2*sin(x_tmp);//;sin(acos(x_theta[i])) *w_theta[i];
        }
        //*/
}

Gauss_Laguerre::~Gauss_Laguerre() {
}

double Gauss_Laguerre::helium_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2) {
    
        double denom = (r1*r1 + r2*r2 -2*(r1*r2*(cos(theta1)*cos(theta2) + sqrt(1-cos(theta1)*cos(theta1))*sqrt(1-cos(theta2)*sin(theta2))*cos(phi1-phi2))));
        double nom = sin(theta1)*sin(theta2);
        if(denom > 1e-15) {
        return nom/(1024.0*sqrt(denom));
        } else return 0;
}

double Gauss_Laguerre::solve(int N_threads) {
        int num_threads = N_threads;
	double *integral = new double[num_threads];
        for (int i = 0; i<num_threads; i++)
            integral[i] = 0;
	
#pragma omp parallel num_threads(num_threads)
{
#pragma omp for
	// integration loop for six dimensional integral
	for (int i = 1; i<N+1; i++) {
	cout << "outer loop step " << i << " / " << N << endl;
	for (int j = 1; j<N+1; j++) 
	for (int k = 0; k<N; k++) 
	for (int l = 0; l<N; l++) 
	for (int m = 0; m<N; m++) 
	for (int n = 0; n<N; n++) {
		integral[omp_get_thread_num()] += w_r[i]*w_r[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[n] * helium_function(x_r[i], x_r[j], x_theta[k], x_theta[l], x_phi[m], x_phi[n]);
		if (omp_get_thread_num() > num_threads) {
			cout << "error in thread numbering!, quitting" << endl;
			exit(1);
		}
	}
        }
#pragma omp barrier
}

	double sum_integral = 0;
	for (int i = 0; i<num_threads; i++) {
		sum_integral += integral[i];
	}
	
	return sum_integral;
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
