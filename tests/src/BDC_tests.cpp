#ifndef _COLORS_
#define _COLORS_

/* FOREGROUND */
#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

#endif  /* _COLORS_ */
#define PASS 1e-7

// Include standard libraries
#include <iostream>
#include <string>

// Include eigen
#include <Eigen/Dense>

// Include custom libraries
#include "Freccia/All.hpp"

#define F 20
#define P 100

Eigen::MatrixXd generateTBDMatrix(unsigned int N, unsigned int f){
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
    
    unsigned int p = N/f;
    // Populate diagonal blocks
    for(int i = 0; i < p; ++i){
        Eigen::MatrixXd B = Eigen::MatrixXd::Random(f, f) * 10;
        B *= B.transpose();
        A.block(i*f, i*f, f, f) = B;
    }

    // Populate off-diagonal blocks
    for(int i = 0; i < p-1; ++i){
        // Eigen::MatrixXd B = Eigen::MatrixXd::Random(f, f) * 10;
        // B *= B.transpose();
        A.block((i+1)*f, i*f, f, f) = Eigen::MatrixXd::Identity(f, f);
        A.block(i*f, (i+1)*f, f, f) = Eigen::MatrixXd::Identity(f, f);
    }

    return A;
}

int test_input(unsigned int t, Eigen::MatrixXd& A, unsigned int f){ 
    
    Freccia::Banded::BDCEigenSolver solver(A, f);
    unsigned int N = A.cols();
    double err = (A - solver.eigenvectors() * solver.eigenvalues().asDiagonal() * solver.eigenvectors().transpose()).norm();
    double orth = (solver.eigenvectors().transpose() * solver.eigenvectors() - Eigen::MatrixXd::Identity(N, N)).norm();

    if(err < PASS && orth < PASS){
        std::cout << BOLD(FGRN("PASSED: "));
    } else {
        std::cout << BOLD(FRED("FAILED: "));
    }
    
    std::cout << std::endl;
    std::cout << "||A - QDQ^T||_F = " << err << std::endl;
    std::cout << "||Q^TQ - I||_F = " << orth << std::endl;

    return err > PASS;
}

// Random matrix
int test1(){
    unsigned int N = 10;
    unsigned int f = 2;

    Eigen::MatrixXd A = generateTBDMatrix(N, f);

    return test_input(1, A, f);
}

int test2(){
    unsigned int N = 2000;
    unsigned int f = 20;

    Eigen::MatrixXd A = generateTBDMatrix(N, f);

    return test_input(2, A, f);
}

int main(int argc, char** argv){
    int rax = 0;
    
    rax |= test1();
    rax |= test2();

    return rax;
}