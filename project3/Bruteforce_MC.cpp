/* 
 * File:   Bruteforce_MC.cpp
 * Author: henriasv
 * 
 * Created on 31. oktober 2012, 00:13
 */

#include "Helium_Solver.h"
#include "Bruteforce_MC.h"
#include "Gauss_Legendre.h"
#include "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/lib/cppLibrary/lib.h"'
#include "omp.h"

Bruteforce_MC::Bruteforce_MC(double range, int N) : Helium_Solver (N){
    l = range;
}


Bruteforce_MC::~Bruteforce_MC() {
}

double Bruteforce_MC::solve(int N_threads) {
    integrals = new double[N_threads];
    double integral = 0;
    int j;
    long idum;
#pragma omp parallel num_threads(N_threads) private(j) private(idum)
{
    j = omp_get_thread_num();
    idum = -j;
#pragma omp for 
    for (int i = 0; i<N; i++) {
        integrals[j] += Helium_Solver::helium_function((2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l, (2*ran0(&idum)-1)*l);
    }
#pragma omp barrier
}
    for (int i = 0; i<N_threads; i++) {
        integral += integrals[i];
    }
    integral /= N;
    integral *= pow(2*l, 6);
    return integral;
}
