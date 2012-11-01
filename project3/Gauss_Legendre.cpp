/* 
 * File:   Gauss_Legendre.cpp
 * Author: henriasv
 * 
 * Created on 30. oktober 2012, 21:10
 */

#include "Gauss_Legendre.h"
#include <cmath>
#include "omp.h"
#include <iostream>

using namespace std;

Gauss_Legendre::Gauss_Legendre(int inN, double a, double b) : Helium_Solver(inN){
	N = inN;
	x = new double[N];
	w = new double[N];
        gauleg(a, b, x, w, N);
}


Gauss_Legendre::~Gauss_Legendre() {
}

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

double Gauss_Legendre::solve(int N_threads) {
        int num_threads = N_threads;
	integrals = new double[num_threads];
	cout << "In Gauss-Legendre solver" << endl;
#pragma omp parallel num_threads(num_threads)
{
#pragma omp for
	// integration loop for six dimensional integral
	for (int i = 0; i<N; i++) {
	cout << "outer loop step " << i+1 << " / " << N << endl;
	for (int j = 0; j<N; j++) {
	for (int k = 0; k<N; k++) {
	for (int l = 0; l<N; l++) {
	for (int m = 0; m<N; m++) {
	for (int n = 0; n<N; n++) {
		integrals[omp_get_thread_num()] += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] * helium_function(x[i], x[j], x[k], x[l], x[m], x[n]);
		if (omp_get_thread_num() > num_threads) {
			cout << "error in thread numbering!, quitting" << endl;
			exit(1);
		}
	}}}}}}
#pragma omp barrier
}

	double sum_integral = 0;
	for (int i = 0; i<num_threads; i++) {
		sum_integral += integrals[i];
	}
	
	return sum_integral;
}
