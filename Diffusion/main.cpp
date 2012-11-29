/* 
 * File:   main.cpp
 * Author: henriasv
 *
 * Created on 24. november 2012, 17:15
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>

#include "Diffusion.h"
#include "ExplicitEuler.h"
#include <string>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    cout << "Running program" << endl;
    string outfile = "/home/scratch/henriasv/FYS3150/Diffusion/test.bin";
    double L = 1.0;
    int N = 10;
    double h = L/(N+1);
    double dt = 0.005;
    Diffusion* solver = new ExplicitEuler(N, L, dt, T);
    solver->set_outfile(outfile);
    solver->set_initial_condition();
    solver->output();
    solver->step();
    solver->output();
    solver->close();
    return 0;
}

