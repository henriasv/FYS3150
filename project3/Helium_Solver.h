/* 
 * File:   Helium_Solver.h
 * Author: henriasv
 *
 * Created on 30. oktober 2012, 21:09
 */

#ifndef HELIUM_SOLVER_H
#define	HELIUM_SOLVER_H

#include <cmath>
#include <iostream>

using namespace std;

class Helium_Solver {
    public:
	// Constructor
	Helium_Solver(int);
	// Destructor
        virtual ~Helium_Solver();
	virtual double solve(int N_threads){cout << "ua" <<endl;};
        // methods
	static double helium_function(double, double, double, double, double, double);
        void gauleg(double, double, double*, double*, int);
	// attributes:
        double *integrals; //pointer to integral array for parallellization 
	double alpha;
	int N;


};

#endif	/* HELIUM_SOLVER_H */

