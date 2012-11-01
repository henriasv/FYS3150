/* 
 * File:   Bruteforce_MC.h
 * Author: henriasv
 *
 * Created on 31. oktober 2012, 00:13
 */

#ifndef BRUTEFORCE_MC_H
#define	BRUTEFORCE_MC_H

#include "Helium_Solver.h"

class Bruteforce_MC  : public Helium_Solver  {
public:
    Bruteforce_MC(double, int);
    Bruteforce_MC(const Bruteforce_MC& orig);
    double solve(int N_cycles);
    virtual ~Bruteforce_MC();
    
    double l;
    double *integrals;
private:

};

#endif	/* BRUTEFORCE_MC_H */

