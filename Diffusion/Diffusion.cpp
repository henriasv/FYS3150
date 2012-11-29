/* 
 * File:   Diffusion.cpp
 * Author: henriasv
 * 
 * Created on 28. november 2012, 20:16
 */

#include <iostream>
#include "Diffusion.h"

using namespace std;

Diffusion::Diffusion(int _N, double _L, double _dt) {
    L = _L;
    N = _N-2;
    h = L/(N+1);
    dt = _dt;
    gamma = dt/h/h;
    
    t = 0;
    tmp_solution = new double[N];
    new_solution = new double[N];
    
    isOutfile = false;
}

void Diffusion::set_initial_condition() {
    for (int i = 0; i<N; i++) {
        tmp_solution[i] = initfunc((i+1)*h);
        new_solution[i] = 0;
    }
}

double Diffusion::initfunc(double x) {
    return 1./x;
}

void Diffusion::output() {
    if (isOutfile) {
    for (int i = 0; i<N; i++) {
        outfiledata->write((char*) &tmp_solution[i], sizeof(double));
    }
    outfiletime->write((char*) &t, sizeof(double));
    }
    else {
        cout << "Outfile not initialized before trying to output. Quitting!" << endl;
        exit(1);
    }
}

void Diffusion::set_outfile(std::string filename) {
    ostringstream outpathdata;
    ostringstream outpathtime;
    ostringstream outpathgrid;
    outpathdata << filename << "/data.bin"; 
    outpathtime << filename << "/time.bin"; 
    outpathgrid << filename << "/grid.bin";
    outfiledata = new std::ofstream(outpathdata.str().c_str(), ios::out | ios::binary);
    outfiletime = new std::ofstream(outpathtime.str().c_str(), ios::out | ios::binary);
    outfilegrid = new std::ofstream(outpathgrid.str().c_str(), ios::out | ios::binary);
    isOutfile = true;
}

void Diffusion::close() {
    outfiledata->close();
    outfiletime->close();
    outfilegrid->close();
}

void Diffusion::solve(double T) {
    while (t<T) {
        step();
        output();
    }
}

Diffusion::~Diffusion() {
    ;
}

