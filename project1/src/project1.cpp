#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

/*
 * Program for solving a tridiagonal matrix problem Ax = v => x=A^-1v
 */
using namespace std;

// forward declarations:
void forwardSubstitution(double *, double *, double *, double *, int);
void backwardSubstitution(double *, double *, double *, double *, int);
void outputMatrix(double *, double *, double *, double *, int);
void rescaleDiagonal(double *, double*, int);
double f(double);
void vectorToFile(double*, int, string);
double solve(double *, double *, double *, double *, int);



int main() {
	string filename = "/scratch/henriasv/FYS3150/project1/";
	string timername = "time.txt";
	string tmp_string = filename+timername;
	cout << filename << endl;
	// timer outfile
	ofstream timer_out(tmp_string.c_str(), ios::out);
	// defining arrays
	int num = 7;
	int N_array[] = {10, 100, 1000, 10000, 100000, 1000000, 10000000};
	for (int i = 0; i<num; i++) {
		int N = N_array[i];
		double h = 1.0/(N+1);
		double* a = new double[N];
		double* b = new double[N];
		double* c = new double[N];
		double* v = new double[N];
		double cpu_time = solve(a, b, c, v, N);
		cout << cpu_time << endl;

		string endString = static_cast<ostringstream*>( &(ostringstream() << N) )->str();
		string filenameOut = filename+endString;
		vectorToFile(v, N, filenameOut);

		string cpu_time_string = static_cast<ostringstream*>( &(ostringstream() << cpu_time) )->str();
		timer_out << endString << " " << cpu_time_string << "\n";

	}
	return 0;
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

void vectorToFile(double *a, int N, string filename) {
	ofstream * outfile = new ofstream(filename.c_str(), ios::out | ios::binary);
	double tmpnum = 0;
	outfile->write((char*) &tmpnum, sizeof(double));
	for (int i = 0; i<N; i++) {
		outfile->write((char*) &a[i], sizeof(double));
	}
	outfile->write((char*) &tmpnum, sizeof(double));
	outfile->close();
}

double f(double x) {
 return 100.0*exp(-10*x);
 }

double solve(double *a, double *b, double *c, double *v, int N) {
	double h = 1.0/(N+1);
	for (int i = 0; i<N; i++) {
		a[i] = -1;
		b[i] = 2;
		c[i] = -1;
		v[i] = h*h*f((i+1)*h);
	}
	time_t start, finish;
	start = clock();

	forwardSubstitution(a, b, c, v, N);
	backwardSubstitution(a, b, c, v, N);
	rescaleDiagonal(b, v, N);
	
	finish = clock();
	return (double)(finish-start)/CLOCKS_PER_SEC;
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
