/* 
 * File:   Diffusion.h
 * Author: henriasv
 *
 * Created on 28. november 2012, 20:16
 */

#ifndef DIFFUSION_H
#define	DIFFUSION_H

#include <fstream>
#include <sstream>
#include <string>

class Diffusion {
public:
    Diffusion(int, double, double);
    virtual ~Diffusion();
    virtual void solve() = 0;
    virtual void set_initial_condition();
    double initfunc(double);
    void output();
    void set_outfile(std::string);
    void close();
    void solve(double);
    virtual void step() = 0;
    
    // Attributes
    int N;
    double L;
    double h;
    double dt;
    double gamma;
    double t;
    double* tmp_solution;
    double* new_solution;
    std::ofstream* outfiledata;
    std::ofstream* outfiletime;
    std::ofstream* outfilegrid;
    bool isOutfile;
    
private:

};

#endif	/* DIFFUSION_H */

