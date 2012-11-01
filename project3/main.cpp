/* 
 * File:   main.cpp
 * Author: henriasv
 *
 * Created on 30. oktober 2012, 21:08
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "Helium_Solver.h"
#include "Gauss_Laguerre.h"
#include "Gauss_Legendre.h"
#include "Bruteforce_MC.h"
#include "Importance_Sampling_MC.h"
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

// function predefinitions
double run_gauleg(int, double, double, int);
double run_gaulag(int, int);
double run_bruteforce_mc(double, int, int);
double run_importance_sampling_mc(int N, int);

/*
 * 
 */
int main(int argc, char** argv) {
    cout << cos(1e-30) << endl;
    cout << "Project 3 - FYS3150" << endl << "-----------------------" << endl;
    int N = atoi(argv[1]);
    double range = atof(argv[2]);
    string method = argv[3];
    int N_threads = atoi(argv[5]);
    double a = -range;
    double b = range;
    double integral;
    cout << "Running " << N_threads << " threads" << endl;
    
    ostringstream op; // output_path
    op << argv[4] << "/" << "#method%" << method << "#N%" << N;
    
    
    //running simulation
    if (method.compare("gauleg") == 0) {
        integral = run_gaulag(N, N_threads);
        op << "#range%" << range;
    }
    else if (method.compare("gaulag") == 0)
        integral = run_gauleg(N, a, b, N_threads);
    else if (method.compare("bf_mc") == 0) {
        integral = run_bruteforce_mc(range, N, N_threads);
        op << "#range%" << range;
    }
    else if (method.compare("is_mc") == 0)
        integral = run_importance_sampling_mc(N, N_threads);
    else
        cout << "Provide integration method as argv[3]" << endl;
    
    ostringstream mkdircmd;
    mkdircmd << "mkdir -p -v " << op.str();
    cout <<mkdircmd.str().c_str() <<endl;
    system(mkdircmd.str().c_str());
    
    op << "/energy.bin";
    ofstream outfile = ofstream(op.str().c_str(), ios::out | ios::binary);
    outfile.write((char*) &integral, sizeof(double));
    
    return 0;
}

double run_gauleg(int N, double a, double b, int N_threads) {
    Gauss_Legendre solver = Gauss_Legendre(N, a, b);
    double integral = solver.solve(N_threads);
    cout << "Solved integral with Gauss-Legendre" << endl << "N = " << N << endl
            << "range [" << a << ", " << b << "]" << endl << "Calculated integral: "
            << integral << endl;
    return integral;
}

double run_gaulag(int N, int N_threads) {
    Gauss_Laguerre solver = Gauss_Laguerre(N);
    double integral = solver.solve(N_threads);
      cout << "Solved integral with Gauss-Laguerre" << endl << "N = " << N << endl
            << "Calculated integral: " << integral << endl;
    return integral;
}

double run_bruteforce_mc(double range, int N, int N_threads) {
    Bruteforce_MC solver = Bruteforce_MC(range, N);
    double integral = solver.solve(N_threads);
    cout << "Solved integral with Brute Force Monte Carlo" << endl << "N_cycles = " << N << endl
            << "Calculated integral: " << integral << endl;
    return integral;
}

double run_importance_sampling_mc(int N, int N_threads) {
    Importance_Sampling_MC solver = Importance_Sampling_MC(N);
    double integral = solver.solve(N_threads);
    cout << "Solved integral with Importance Sampling Monte Carlo" << endl << "N_cycles = " << N << endl
            << "Calculated integral: " << integral << endl;
    return integral;
}