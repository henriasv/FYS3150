/* 
 * File:   IsingLattice2D.h
 * Author: henriasv
 *
 * Created on 14. november 2012, 01:13
 */

#ifndef ISINGLATTICE2D_H
#define	ISINGLATTICE2D_H

#include "../lib/cppLibrary/lib.h"
#include <armadillo>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

using namespace arma;


/** The 2 dimensional Ising model
 * Units: 
 * 
 * @param N
 * @param n_threads
 * @param rng_seed
 */
class IsingLattice2D {
public:
    /**
     * There are parameters passed to this function that are not in use: num_cycles.
     * The cycling is taken from outside, by calling th Metropolis function from global scope.
     * @param _N
     * @param T
     * @param num_cycles
     * @param rng_seed
     */
    IsingLattice2D(int _N, double T, int num_cycles, long int rng_seed);
    virtual ~IsingLattice2D();
    void Metropolis();
    
    // Initializing functions
    void initializeUp();
    void initializeDown();
    void initializeRandom();
    
    
    // Measuring functions
    void calculateEnergy();
    void calculateMagnetization();
    
    // MISC
    void print();
    int randomSpin();
    int periodic(int);
    void flipSpin(int, int, bool);
    void set_output_folder(string);
    void close();
    
    /**
     * Output data to files. Folder has to be specified;
     */
    void output();
    
    
protected:
    bool isInitialized;
    mat spin_matrix;
    int num_cycles;
    long int idum;
    int N;
    double E;
    double M;
    double T;
    int accepted;
    int rejected;
    
    // Hardcoded acceptance probabilities for positive energy changes
    double p4;
    double p8;
    
    string output_folder;
    ofstream out_E;
    ofstream out_M;
    ofstream out_Esq;
    ofstream out_Msq;
    ofstream out_absM;
};

#endif	/* ISINGLATTICE2D_H */

