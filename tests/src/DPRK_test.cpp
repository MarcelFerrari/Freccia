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

#define PASS 1e-8

// Include standard libraries
#include <iostream>
#include <string>

// Include eigen
#include <Eigen/Dense>

// Include custom libraries
#include "Freccia/All.hpp"

int test_input(unsigned int t, Eigen::VectorXd& D, Eigen::MatrixXd& W, double rho = 1.0){

    std::cout << KYEL << "Running DPRK test " << t << "..." << RST << std::endl;

    Eigen::MatrixXd A = D.asDiagonal();
    A += rho * W * W.transpose();

    Freccia::DPRK::DPRKEigenSolver solver(D, W, rho);

    double err = (A - solver.eigenvectors() * solver.eigenvalues().asDiagonal() * solver.eigenvectors().transpose()).norm();
    double orth = (solver.eigenvectors().transpose() * solver.eigenvectors() - Eigen::MatrixXd::Identity(D.size(), D.size())).norm();
    
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

int test1(){
    unsigned int n = 200;
    unsigned int f = 15;
    Eigen::VectorXd D =  Eigen::VectorXd::Random(n);
    Eigen::MatrixXd W = Eigen::MatrixXd::Random(n, f);

    return test_input(1, D, W);
}


int test2(){
    unsigned int n = 10;
    unsigned int f = 5;
    Eigen::VectorXd D = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd W = Eigen::MatrixXd::Random(n, f);

    // Making some elements in the columns of W matrix the same
    W.col(0).head(5) = W(0, 0) * Eigen::VectorXd::Ones(5);
    W.col(1).head(5) = W(1, 1) * Eigen::VectorXd::Ones(5);
    W.col(2).head(5) = W(2, 2) * Eigen::VectorXd::Ones(5);
    W.col(3).head(5) = W(3, 3) * Eigen::VectorXd::Ones(5);

    return test_input(2, D, W);
}


int main(int argc, char** argv){
    int rax = 0;
    
    rax |= test1();
    rax |= test2();
    
    return rax;
}