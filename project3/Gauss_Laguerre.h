/* 
 * File:   Gauss_Laguerre.h
 * Author: henriasv
 *
 * Created on 30. oktober 2012, 21:10
 */

#ifndef GAUSS_LAGUERRE_H
#define	GAUSS_LAGUERRE_H

#include "Helium_Solver.h"

class Gauss_Laguerre : public Helium_Solver {
public:
	/**
	* N is the grid resolution, limit is the upper limit of integration on the r-axis.	
	* Angles run from 0 to pi and 0 to 2pi
	*/
	Gauss_Laguerre(int N);
        double solve(int);
        virtual ~Gauss_Laguerre();
protected:
	void gaulag(double *, double *, int, double);
	double gammln(double);
	
	double helium_function(double, double, double, double, double, double);
	double alph;
        double *x_r;
	double *w_r;
	double *x_phi;
	double *w_phi;
	double *x_theta;
	double *w_theta;
};

#endif	/* GAUSS_LAGUERRE_H */

