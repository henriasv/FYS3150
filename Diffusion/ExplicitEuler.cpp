/* 
 * File:   ExplicitEuler.cpp
 * Author: henriasv
 * 
 * Created on 28. november 2012, 20:28
 */

#include "ExplicitEuler.h"

ExplicitEuler::ExplicitEuler(int N, double dt, double L) : Diffusion(N, dt, L){
    
}

ExplicitEuler::~ExplicitEuler() {
}

void ExplicitEuler::solve() {
    ;
}

void ExplicitEuler::step() {
    new_solution[0] = gamma*tmp_solution[1] + (1-2*gamma)*tmp_solution[0];
    for (int i = 1; i<N-1; i++) {
        new_solution[i] = gamma*(tmp_solution[i+1]+tmp_solution[i-1]) + (1-2*gamma)*tmp_solution[i];
    }
    new_solution[N-1] = gamma*tmp_solution[N-2] + (1-2*gamma)*tmp_solution[N-1];
    double * tmp = tmp_solution;
    tmp_solution = new_solution;
    new_solution = tmp_solution;
}
