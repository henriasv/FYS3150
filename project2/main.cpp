/* 
 * File:   main.cpp
 * Author: henriasv
 *
 * Created on 13. oktober 2012, 20:04
 */

//#define DEBUG

#include <cstdlib>
#include <armadillo>
#include <cmath>
#include <sys/stat.h>

using namespace std;
using namespace arma;

// Function predefinitions
void jacobiRotation(mat &A, mat &B, vec indices, vec values);
double potential(double rho, double omega);
double diagonalElement(double rho, double omega, double h);
double offDiagonalElement(double h);
vec maxOffdiagonalElement(mat A);
vec findAngle(mat A, vec indices);
void vectorToFile(vec a, int N, string filename);

/*
 * 
 */
int main(int argc, char** argv) {
    cout << "Program  started";
    vec indices, values, rho;
    ostringstream outpath;
    outpath << "/scratch/henriasv/FYS3150/project2/with_interaction/";
    
    double omega = 5;
    int N = 100;
    double epsilon = 2.3e-6;
    double rho_max = 2;
    double rho_min = rho_max/(N+1);
    double h = (rho_max-rho_min)/N;
 
    // Create directory for data output
    outpath << "N=" << N << "/" << "omega=" << omega << "/";
    const char *dir = outpath.str().c_str();
    bool isCreatedDirectory = (mkdir(dir, 0777) == 0);
    if (isCreatedDirectory) {
        cout << "Successfully created directory " << outpath.str() << endl;
    }
    else {
        cout << "could not create directory " << outpath.str() <<
                   "\nDo you wish to continue? (old data may be overwritten)" << endl;
        string evalIsContinue;
        cin >> evalIsContinue;
        bool isContinue = (evalIsContinue.compare("yes")==0);
        if (!isContinue) {
            cout << "quitting" << endl;
            exit(0);
        }
    }
    
    // Prepare for output of sum of squared off-diagonal elements
    ofstream sqsum_output;
    sqsum_output.open((outpath.str() + "sqsum.dat").c_str());
    
    // Preparing matrix
    mat A = mat(N, N);
    cout << endl;
    A.zeros();
    A.print();
    
    // Matric for collecting the transformation to find eigenvectors.
    mat B = eye(N, N);
    rho = vec(N);
    
    for (int i = 1; i<N-1; i++) {
        rho(i) = rho_min + i*h;
        A(i, i) = diagonalElement(rho_min+i*h, omega, h);
        A(i, i+1) = offDiagonalElement(h);
        A(i, i-1) = offDiagonalElement(h);
    }
    rho(0) = rho_min;
    rho(N-1) = rho_min + (N-1)*h;
    A(0, 0) = diagonalElement(rho_min, omega, h);
    A(N-1, N-1) = diagonalElement(rho_min+(N+1)*h,omega, h);
    A(0, 1) = offDiagonalElement(h);
    A(N-1, N-2) = offDiagonalElement(h);
    // Tridiagonal matrix created
    
    // Solution loop
    int n_iterations = 0;
    double sqsum = epsilon+1; // dummy value
    while (sqsum > epsilon) {
        indices = maxOffdiagonalElement(A);
        values = findAngle(A, indices);
        jacobiRotation(A, B, indices, values);
        sqsum = indices(2);
        cout << "sqsum " << sqsum << endl;
        n_iterations ++;
        sqsum_output << sqsum << endl;
    }
    
    // Output
    vec diagonal = A.diag();
    uvec indices_out = sort_index(diagonal);
    indices_out.print();
    
    vec eigenvalues = vec(3);
    eigenvalues.zeros();
    int n_eig = 3;
    for (int i = 0; i<n_eig; i++) {
        eigenvalues(i) = diagonal(indices_out(i));
    }
    
    // Dump eigenvalues
    ostringstream outeig_path;
    outeig_path << outpath.str() << "eigenvalues.bin";
    vectorToFile(eigenvalues, eigenvalues.size(), outeig_path.str());
    
    
    // Dump grid
    ostringstream outgrid_path;
    outgrid_path << outpath.str() << "grid.bin";
    vectorToFile(rho, N, outgrid_path.str());
    
    // Dump data from specified states:
    int N_states_output = 3;
    for (int i = 0; i<N_states_output; i++) {
        ostringstream outstate_path;
        outstate_path << outpath.str();
        outstate_path << "state" << i << ".bin";
        cout << "writing state to file " << outstate_path.str() << endl;
        vectorToFile(B.col(indices_out(i)), N, outstate_path.str());
    }
    
    // Dump info to txt file
    ofstream txtFile;
    string resultsFileName("results.txt");
    txtFile.open((outpath.str() + resultsFileName).c_str());
    txtFile << "epsilon = " << epsilon << "\n";
    txtFile << "n_iterations = " << n_iterations <<"\n";
    txtFile << "N = " << N << "\n";
    txtFile << "rho_min = " << rho_min << "\n";
    txtFile << "rho_max = " << rho_max << "\n";
    txtFile << "omega_r = " << omega << "\n";
    txtFile << "eigenvalues:\n";
    for (int i= 0; i<N_states_output; i++) {
        txtFile << "lambda"<< i << " = " << diagonal(indices_out(i)) << "\n";
    }
    return 0;
    
}



void jacobiRotation(mat &A, mat &B, vec indices, vec values) {
    uword p = indices(0);
    uword q = indices(1);
    mat S = eye(A.n_rows, A.n_rows);
    double sin_theta = values(0);
    double cos_theta = values(1);
    S(p, p) = cos_theta;
    S(q, q) = cos_theta;
    S(p, q) = sin_theta;
    S(q, p) = -sin_theta;
    A = S.st()*A*S;
    B = B*S;
}

vec findAngle(mat A, vec indices) {
    uword p, q;
    double tan_theta;
    vec values = vec(2);
    p = indices(0);
    q = indices(1);
    
    double tau = (A(q,q)-A(p, p))/(2*A(p, q));
    
    if (tau>=0) {
        tan_theta = -tau + sqrt(1+tau*tau);
    }
    else {
        tan_theta = -tau - sqrt(1+tau*tau);
    }
    
    double cos_theta = 1.0/sqrt(1+tan_theta*tan_theta);
    double sin_theta = tan_theta*cos_theta;
    values(0) = sin_theta;
    values(1) = cos_theta;
    return values;
}

vec maxOffdiagonalElement(mat A) {
    uword m, n;
    vec v = vec(3);
    mat C = A%A;
    
    // Removing diagonal
    for (int i = 0; i<C.n_rows; i++) {
        C(i, i) = 0;
    }
    C.max(m, n);
    v(0) = m;
    v(1) = n;
    v(2) = accu(C); // Sum the elements
#ifdef DEBUG
    cout << "Found maximum off diagonal element at" << endl << "["<< endl;
    v.print();
    cout << "]" << endl;
#endif
    return v;
}


double potential(double rho, double omega) {
    return rho*rho*omega*omega + 1.0/rho;
}

double diagonalElement(double rho, double omega, double h) {
    return potential(rho, omega) + 2.0/(h*h);
}

double offDiagonalElement(double h) {
    return -1.0/(h*h);
}

void vectorToFile(vec a, int N, string filename) {
        ofstream * outfile = new ofstream(filename.c_str(), ios::out | ios::binary);
        double tmpnum = 0;
        //outfile->write((char*) &tmpnum, sizeof(double));
        for (int i = 0; i<N; i++) {
                outfile->write((char*) &a(i), sizeof(double));
        }
        //outfile->write((char*) &tmpnum, sizeof(double));
        outfile->close();
}