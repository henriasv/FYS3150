/* 
 * File:   Gauss_Legendre.h
 * Author: henriasv
 *
 * Created on 30. oktober 2012, 21:10
 */

#ifndef GAUSS_LEGENDRE_H
#define	GAUSS_LEGENDRE_H

#include "Helium_Solver.h"

class Gauss_Legendre : public Helium_Solver {
public:
        Gauss_Legendre(int N, double a, double b);
        virtual ~Gauss_Legendre();
        double solve(int);
	double helium_function(double, double, double, double, double, double);
        
        double *w;
	double *x;
        

};

#endif	/* GAUSS_LEGENDRE_H */

