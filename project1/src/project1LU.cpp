#include <stdio.h>
#include <stdlib.h>
#include "lib.h"
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

void vectorToFile(double *, int, string);
double f(double);

int main() {
	string filename = "/scratch/henriasv/FYS3150/project1/";
	string timername = "timeLU.txt";
	string tmp_string = filename+timername;
	cout << filename << endl;
	// timer outfile
	ofstream timer_out(tmp_string.c_str(), ios::out);
	// defining arrays
	int num = 3;
	int N_array[] = {3000};//, 1000};//, 1000};//, 10000, 100000};
	for (int r = 0; r<num; r++) {
		int N = N_array[r];
		double h = 1.0/(N+1);
		double * f_val = new double[N];
		double ** A = (double**) matrix(N, N, sizeof(double));
		int * indx = new int[N];
		double d;
		cout << endl;

		A[0][0] = 2.0;
		A[0][1] = -1.0;
		for (int i = 1; i<N-1; i++) {
			A[i][i-1] = -1.0;
			A[i][i] = 2.0;
			A[i][i+1] = -1.0;
		}
		A[N-1][N-1] = 2.0;
		A[N-1][N-2] = -1.0;

		time_t start, finish;
		start = clock();
		ludcmp(A, N, indx, &d);
		lubksb(A, N, indx, f_val);
		finish = clock();

		cout << "N = " << N << "; " << "Time spent: " << (double)(finish-start)/CLOCKS_PER_SEC << endl;

		ostringstream outname;
		outname << filename << "LU" << N;
		vectorToFile(f_val, N, outname.str());

		//string cpu_time_string = static_cast<ostringstream*>( &(ostringstream() << cpu_time) )->str();
		//timer_out << endString << " " << cpu_time_string << "\n";
		free_matrix((void**) A);
		delete [] indx;
		delete [] f_val;
		
	}
	return 0;
}

double f(double x) {
 return 100.0*exp(-10*x);
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
