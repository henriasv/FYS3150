/* 
 * File:   main.cpp
 * Author: henriasv
 *
 * Created on 8. november 2012, 16:37
 */

#include <cstdlib>
#include "IsingLattice2D.h"
#include <iostream>
#include <sstream>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    ostringstream basepath;
    basepath << "/scratch/henriasv/FYS3150/IsingModel_once/";
    int N = atoi(argv[1]);
    int n_mcs = atoi(argv[3]);
    int n_threads = 1;
    long int idum = -2;
    double T = atof(argv[2]);
    
    // Set output folder ;
    basepath << "N" << N << "/T" << T;
    ostringstream mkdircmd;
    mkdircmd << "mkdir -p -v " << basepath.str();
    system(mkdircmd.str().c_str());
    
    IsingLattice2D* lattice = new IsingLattice2D(N, T, n_threads, idum);
    lattice->set_output_folder(basepath.str());
    lattice->initializeUp();
    lattice->print();
    for (int i = 0; i<n_mcs; i++) {
        lattice->Metropolis();
        lattice->output();
    }
    lattice->close();
    return 0;
}

