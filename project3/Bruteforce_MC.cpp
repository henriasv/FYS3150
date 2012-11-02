/* 
 * File:   Bruteforce_MC.cpp
 * Author: henriasv
 * 
 * Created on 31. oktober 2012, 00:13
 */

#include "Bruteforce_MC.h"
#include "Gauss_Legendre.h"
#include "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/lib/cppLibrary/lib.h"
#include "omp.h"
#define EPS 1e-14

Bruteforce_MC::Bruteforce_MC(double range, int _N) : Helium_Solver (_N){
    l = range;
    N = _N;
}


Bruteforce_MC::~Bruteforce_MC() {
}

double Bruteforce_MC::solve(int N_threads) {
    integrals = new double[N_threads];
    for (int i = 0; i<N_threads; i++) {
        integrals[i] = 0;
    }
    double integral = 0;
    int j;
    long idum;
#pragma omp parallel num_threads(N_threads) private(j) private(idum)
{
    double tmpSum = 0;
    j = omp_get_thread_num();
    idum = -j;
#pragma omp for 
    for (int i = 0; i<N; i++) {
        tmpSum += Helium_Solver::helium_function((2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l);
        if (tmpSum/(integrals[j]+EPS) > 0.001) {
                integrals[j] += tmpSum;
                tmpSum = 0;
        }
    }
#pragma omp critical
    {
        integral += integrals[j];
        integral += tmpSum;
        cout << integral << endl;
    }
}
    integral /= N;
    cout << integral << endl;
    integral *= pow(2*l, 6);
    cout << integral << endl;

    return integral;
}
