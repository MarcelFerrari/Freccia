#include <Eigen/Dense>
#include "Freccia/banded/BandedEigenSolver.hpp"
#include <iostream>

// Declaration of the LAPACK dsbevd function
// extern "C" {
// void dsbevd_(char* jobz, char* uplo, int* n, int* kd, double* ab, int* ldab, 
//              double* w, double* z, int* ldz, double* work, int* lwork, 
//              int* iwork, int* liwork, int* info);
// }

// Include relevant headers based on which library is used
#ifdef FRECCIA_USE_MKL
    #include "mkl.h"
    #define LAPACK_INT MKL_INT
#elif defined(FRECCIA_USE_OPENBLAS)
    #include "lapacke.h"
    #define LAPACK_INT lapack_int  
#else
    #error "Freccia relies on either MKL or OpenBLAS for the ?SBEVD routine. Otherwise, set FRECCIA_NO_LAPACK to true."
#endif

void Freccia::Banded::BandedEigenSolver::dsbevd_wrapper() {
    char jobz = 'V'; // Request eigenvalues and eigenvectors
    LAPACK_INT n = A_band.cols(); // Order of the matrix
    LAPACK_INT kd = f-1; // Number of superdiagonals of the band matrix
    LAPACK_INT ldab = A_band.rows(); // Leading dimension of the array AB
    double* ab = A_band.data();
    double* w = ew.data();
    double* z = ev.data();
    LAPACK_INT ldz = n;
    LAPACK_INT lwork = -1; // Set to -1 initially to query optimal workspace size
    double wkopt;
    LAPACK_INT iwork_query;
    LAPACK_INT liwork = -1; // Set to -1 initially to query optimal workspace size
    LAPACK_INT info = 0; // Info about the computation process

    // Query and allocate the optimal workspace
    dsbevd(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, &wkopt, &lwork, &iwork_query, &liwork, &info);
    
    // Allocate the optimal workspace
    lwork = wkopt;
    std::vector<double> work(lwork); // using vector instead of raw pointer
    
    liwork = iwork_query;
    std::vector<LAPACK_INT> iwork(liwork); // using vector instead of raw pointer

    // Reset info
    info = -1;
    
    // Calculate eigenvalues and eigenvectors
    dsbevd(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work.data(), &lwork, iwork.data(), &liwork, &info);

    // Check for computation errors
    if (info != 0) {
        // Handle errors
        std::cerr << "The algorithm failed to compute eigenvalues and eigenvectors." << std::endl;
    }
}