/* 
 * File:   Importance_Sampling_MC.h
 * Author: henriasv
 *
 * Created on 31. oktober 2012, 02:41
 */

#ifndef IMPORTANCE_SAMPLING_MC_H
#define	IMPORTANCE_SAMPLING_MC_H

#include "Helium_Solver.h"

using namespace std;

class Importance_Sampling_MC : public Helium_Solver{
public:
    Importance_Sampling_MC(int);

    virtual ~Importance_Sampling_MC();
    double solve(int);
    double helium_function(double r1, double r2, double cos1, double cos2, double phi1, double phi2);
    double exp_trans(double);
    double exp_func(double);
    
private:

};

#endif	/* IMPORTANCE_SAMPLING_MC_H */

