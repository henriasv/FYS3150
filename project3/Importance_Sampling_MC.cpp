/* 
 * File:   Importance_Sampling_MC.cpp
 * Author: henriasv
 * 
 * Created on 31. oktober 2012, 02:41
 */

#include "Importance_Sampling_MC.h"
#define PI 3.1415926536
#define EPS 1e-14
#include "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/lib/cppLibrary/lib.h"'
#include "omp.h"


Importance_Sampling_MC::Importance_Sampling_MC(int N) : Helium_Solver(N) {
    this->N = N;
}


Importance_Sampling_MC::~Importance_Sampling_MC() {
}


double Importance_Sampling_MC::helium_function(double r1, double r2, double cos1, double cos2, double phi1, double phi2) {
    
        double denom = 4*4*sqrt((r1*r1 + r2*r2 -2*(r1*r2*(cos1*cos2 + sqrt(1-cos1*cos1)*sqrt(1-cos2*cos2)*cos(phi1-phi2)))));
        double nom = r1*r1*r2*r2;
        if(denom > 1e-10) {
        return nom/denom;
        } else return 0;
}

double Importance_Sampling_MC::solve(int N_threads) {
    integrals = new double[N_threads];
    double integral = 0;
    int j;
    long idum;
    double x1;
    double x2;
    double tmpSum = 0;
#pragma omp parallel num_threads(N_threads) private(j) private(idum) private(x1) private(x2) private(tmpSum)
{
    j = omp_get_thread_num();
    idum = -j;
#pragma omp for 
    for (int i = 0; i<N; i++) {
        // drawing numbers
        x1 = ran0(&idum);
        x2 = ran0(&idum);
        tmpSum += helium_function(exp_trans(x1), exp_trans(x2), 2*ran0(&idum)-1, 2*ran0(&idum)-1, ran0(&idum)*2*PI, ran0(&idum)*2*PI);
        
        if (tmpSum/(integrals[j]+EPS) > 0.001) {
            integrals[j] += tmpSum;
            tmpSum = 0;
        }
    }
#pragma omp critical
    {
        integral += integrals[j];
        integral += tmpSum;
    }
    integral /= N;
    integral *= 2*2*2*PI*2*PI;
}
    return integral;
}

double Importance_Sampling_MC::exp_trans(double rannum) {
    return -log(1-rannum)/4;
}

double Importance_Sampling_MC::exp_func(double x) {
    return 4*exp(-4*x);
}