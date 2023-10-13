#ifndef FRECCIA_ARROWHEAD_MATRIX_H
#define FRECCIA_ARROWHEAD_MATRIX_H

// Include Eigen
#include <Eigen/Dense>

// Struct to store Arrowhead matrices
// Templated for double and quadruple precision support
template<typename T>
struct ArrowheadMatrix {
    Eigen::Array<T, Eigen::Dynamic, 1> D;
    Eigen::Array<T, Eigen::Dynamic, 1> w;
    T b;
    unsigned int i; // Index of the arrowhead

    // Default constructor
    ArrowheadMatrix() : b(T(0)), i(0) {}

    // Empty constructor
    ArrowheadMatrix(unsigned int N) : i(0), b(T(0)) {
        D = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(N);
        w = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(N);
    }

    // Parameterized constructor
    ArrowheadMatrix(const Eigen::Array<T, Eigen::Dynamic, 1>& D_in, 
                    const Eigen::Array<T, Eigen::Dynamic, 1>& w_in, 
                    T b_in, 
                    unsigned int i_in)
        : D(D_in), w(w_in), b(b_in), i(i_in) {}
    

    // Assignment operator for deep copy
    ArrowheadMatrix& operator=(const ArrowheadMatrix& other) {
        if (this != &other) {
            D = other.D;
            w = other.w;
            b = other.b;
            i = other.i;
        }
        return *this;
    }
};

#endif