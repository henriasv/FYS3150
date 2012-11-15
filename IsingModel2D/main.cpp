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
    basepath << "/scratch/henriasv/FYS3150/IsingModel/";
    int N = 20;
    int n_threads = 1;
    long int idum = -2;
    double T = 2.0;
    
    // Set output folder ;
    basepath << "N" << N << "/T" << T;
    
    IsingLattice2D* lattice = new IsingLattice2D(N, T, n_threads, idum);
    lattice->initializeRandom();
    lattice->print();
    for (int i = 0; i<100; i++) {
        lattice->Metropolis();
        lattice->print();
    }
    return 0;
}

