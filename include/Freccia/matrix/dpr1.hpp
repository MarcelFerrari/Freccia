#ifndef FRECCIA_DPR1_MATRIX_H
#define FRECCIA_DPR1_MATRIX_H

// Include Eigen
#include<Eigen/Dense>

// Struct to store DPR1 matrices
template<typename T>
struct DPR1Matrix {
    Eigen::Array<T, Eigen::Dynamic, 1> D;
    Eigen::Array<T, Eigen::Dynamic, 1> z;
    T rho;

    // Default constructor
    DPR1Matrix() : rho(T(0)) {}

    // Empty constructor
    DPR1Matrix(unsigned int N) : rho(T(0)) {
        D = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(N);
        z = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(N);
    }

    // Parameterized constructor
    DPR1Matrix(const Eigen::Array<T, Eigen::Dynamic, 1>& D_in, 
               const Eigen::Array<T, Eigen::Dynamic, 1>& z_in, 
               T rho_in) 
        : D(D_in), z(z_in), rho(rho_in) {}

    // Assignment operator for deep copy
    DPR1Matrix& operator=(const DPR1Matrix& other) {
        if (this != &other) {
            D = other.D;
            z = other.z;
            rho = other.rho;
        }
        return *this;
    }
};


#endif