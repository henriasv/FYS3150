#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

/*
 * Program for solving a tridiagonal matrix problem Ax = v => x=A^-1v
 */

// forward declarations:
void forwardSubstitution(double *, double *, double *, double *, int);
void backwardSubstitution(double *, double *, double *, double *, int);
void outputMatrix(double *, double *, double *, double *, int);
void rescaleDiagonal(double *, double*, int);
double f(double);

using namespace std;

int main() {
	// defining arrays
	int N = 35;
	double h = 1.0/(N+1);
	double* a = new double[N];
	double* b = new double[N];
	double* c = new double[N];
	double* v = new double[N];
	
	for (int i = 0; i<N; i++) {
		a[i] = -1;
		b[i] = 2;
		c[i] = -1;
		v[i] = h*h*f((i+1)*h);
	}
	outputMatrix(a, b, c, v, N);
	forwardSubstitution(a, b, c, v, N);
	outputMatrix(a, b, c, v, N);
	backwardSubstitution(a, b, c, v, N);
	outputMatrix(a, b, c, v, N);
	rescaleDiagonal(b, v, N);
	outputMatrix(a, b, c, v, N);
	
}

void forwardSubstitution(double *a, double *b, double *c, double *v, int N) {
	for (int i = 1; i<N; i++) {
		b[i] = b[i]-c[i-1]*a[i]/b[i-1];
		v[i] = v[i]-v[i-1]*a[i]/b[i-1];
		a[i] = a[i]-b[i-1]*a[i]/b[i-1]; // = 0 To be removed, but kept for testing!
	}
}

void backwardSubstitution(double *a, double *b, double *c, double *v, int N) {
	for (int i = N-1; i>=1; i--) {
		v[i-1] = v[i-1] - v[i]*c[i-1]/b[i];
		b[i-1] = b[i-1] - a[i]*c[i-1]/b[i];
		c[i-1] = c[i-1] - b[i]*c[i-1]/b[i];
	}
}
void rescaleDiagonal(double *b, double *v, int N) {
	for (int i = 0; i<N; i++) {
		v[i] = v[i]/b[i];
		b[i] = b[i]/b[i];
	}
}

double f(double x) {
 return 100.0*exp(-10*x);
 }

void outputMatrix(double *a, double *b, double *c, double *v, int N) {
	cout << endl;
	for (int i = 0; i<N; i++) {
		int j;
		for (j = 0; j<i-1; j++) {
			printf("%6.2f", 0.0);
		}
		if(i>=1 && i<N-1) {
			printf("%6.2f%6.2f%6.2f", a[i], b[i], c[i]);
		}
		else if (i<1){
			printf("%6.2f%6.2f", b[i], c[i]);
			j--;
		}
		else if (i>=N-1) {
			printf("%6.2f%6.2f", a[i], b[i]);
		}
		for (j+=3; j<N; j++) {
			printf("%6.2f", 0.0);
		}
		printf("%6.2f\n", v[i]);
	}
}