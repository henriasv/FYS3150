/* 
 * File:   IsingLattice2D.cpp
 * Author: henriasv
 * 
 * Created on 14. november 2012, 01:13
 */

#include "IsingLattice2D.h"


IsingLattice2D::IsingLattice2D(int _N, double _T, int _num_cycles, long int rng_seed) {
    E = 0;
    M = 0;
    accepted = 0;
    rejected = 0;
    isInitialized = false;
    N = _N;
    T = _T;
    num_cycles = _num_cycles;
    idum = rng_seed;
    spin_matrix = mat(N, N);
    p4 = exp(-4.0/T);
    p8 = exp(-8.0/T);
    output_folder = "False";
}


IsingLattice2D::~IsingLattice2D() {
}

void IsingLattice2D::Metropolis() {
    for (int n = 0; n<N*N; n++) {
        int i = (int) (ran1(&idum)*N);
        int j = (int) (ran1(&idum)*N);
        flipSpin(i, j);
    }
}

void IsingLattice2D::initializeDown() {
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            spin_matrix(i, j) = -1;
    calculateEnergy();
    isInitialized = true;
}

void IsingLattice2D::initializeUp() {
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            spin_matrix(i, j) = 1;
    calculateEnergy();
    isInitialized = true;
}

void IsingLattice2D::initializeRandom() {
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++) {
            spin_matrix(i, j) = this->randomSpin();
        }
    calculateEnergy();
    isInitialized = true;
}

void IsingLattice2D::print() {
    //spin_matrix.print();
    cout << "Energy = " << E << endl;
    cout << "Magnetization = " << M << endl;
    cout << "Accepted moves = " << accepted << endl;
    cout << "Rejected moves = " << rejected << endl;
}

int IsingLattice2D::randomSpin() {
    bool tmp = (ran1(&idum)> 0.5);
    return (int) (tmp - !tmp);
}

void IsingLattice2D::calculateEnergy() {
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++) {
            E -=    spin_matrix(i,j)*(
                    spin_matrix(i, periodic(j+1)) +
                    spin_matrix(i ,periodic(j-1)) +
                    spin_matrix(periodic(i+1), j) +
                    spin_matrix(periodic(i-1), j));
        }
}

void IsingLattice2D::calculateMagnetization() {
    M = sum(sum(spin_matrix));
}

/**
 * Transforming indices outside the lattice corresponding to periodic boundary conditions.
 * NOTE: Works for nearest neighbor calls only! Crashes the program for requests more than one step out of the grid.
 * @param i Index of spin to access
 * @param N Size of lattice
 * @return Index to be used
 */
int IsingLattice2D::periodic(int i) {
    if (i<N && i>=0)
        return i;
    else if (i == N)
        return 0;
    else if (i == -1)
        return N-1;
    else {
        cout << "Index larger than N, quitting!" << endl;
        exit(1);
    }
}


/**
 * Tries to flip spin i, j. If accepted, the energy and magnetization is updated
 * 
 * @param i
 * @param j
 */
void IsingLattice2D::flipSpin(int i, int j) {
    double dE =     2*spin_matrix(i,j)*(
                    spin_matrix(i, periodic(j+1)) +
                    spin_matrix(i ,periodic(j-1)) +
                    spin_matrix(periodic(i+1), j) +
                    spin_matrix(periodic(i-1), j));
    if (dE <= 0) {
        M += -2*spin_matrix(i, j);
        E += dE;
        spin_matrix(i,j) *= -1;
        accepted ++;
    }
    else {
        // Check for acceptance;
        if (dE == 4) {
            if (ran1(&idum) <= p4) {
                // accept move
                M += -2*spin_matrix(i, j);
                E += dE;
                spin_matrix(i, j) *= -1;
                accepted ++;
            } else rejected ++;
        }
        if (dE == 8) {
            if (ran1(&idum) <= p8) {
                // Accept move
                M += -2*spin_matrix(i, j);
                E += dE;
                spin_matrix(i, j) *= -1;
                accepted ++;
            } else rejected ++;
        }
    }
}

void IsingLattice2D::output() {
    if (!(output_folder.compare("False") == 0)) {
        
    }
}

void IsingLattice2D::set_output_folder(string folder) {
    ostringstream path_M, path_E, path_Esq, path_Msq, path_absM;
    path_E << folder << "/E.bin";
    path_M << folder << "/M.bin";
    path_Esq << folder << "/Esq.bin";
    path_Msq << folder << "/Msq.bin";
    path_absM << folder << "/absM.bin";
    
    out_E.open(path_E);
    out_M.open(path_M);
    out_Esq.open(path_Esq);
    out_Msq.open(path_Msq);
    out_absM.open(path_absM);
}
