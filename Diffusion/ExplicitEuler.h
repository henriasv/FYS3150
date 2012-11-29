/* 
 * File:   ExplicitEuler.h
 * Author: henriasv
 *
 * Created on 28. november 2012, 20:28
 */

#ifndef EXPLICITEULER_H
#define	EXPLICITEULER_H

#include "Diffusion.h"

class ExplicitEuler : public Diffusion{
public:
    ExplicitEuler(int, double, double);
    virtual ~ExplicitEuler();
    
    void step();
    virtual void solve();
private:

};

#endif	/* EXPLICITEULER_H */

