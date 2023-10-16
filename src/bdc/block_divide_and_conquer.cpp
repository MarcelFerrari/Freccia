// Standard libraries
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <tuple>
#include <limits>

// Include eigen
#include <Eigen/Dense>

// Custom headers
#include "Freccia/bdc/BDCEigenSolver.hpp"

// Convenience aliases
using SVDTuple = std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>;
using EigenTuple = std::tuple<unsigned int, unsigned int>;

void Freccia::Banded::BDCEigenSolver::block_divide_and_conquer(const Eigen::MatrixXd &A, const unsigned int f_in) {    
    // Initialize Eigen parallel environment
    Eigen::initParallel();

    // Set up problem variables
    N = A.cols();
    f = f_in;
    p = N/f;

    // Allocate memory
    D = Eigen::VectorXd::Zero(N);
    Q = Eigen::MatrixXd::Zero(N, N);

    // Compute SVD decompositions of off-diagonal blocks
    svdBlocks = std::vector<SVDTuple>(p-1);
    #pragma omp parallel for
    for(int i = 0; i < p-1; ++i) {
        // Only need thin V and U
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A.block((i+1)*f, i*f, f, f), Eigen::ComputeThinU | Eigen::ComputeThinV);
        svdBlocks[i] = std::make_tuple(std::move(svd.matrixU()), std::move(svd.singularValues()), std::move(svd.matrixV()));
    }

    // Diagonalize diagonal blocks in parallel
    // Results are stored in Q and D
    std::vector<EigenTuple> diag_blocks(p);
    #pragma omp parallel for
    for(int i = 0; i < p; ++i){
        const Eigen::MatrixXd& Ai = A.block(i*f, i*f, f, f);
        if(i == 0){ // First block
            // Get required matrices
            const Eigen::VectorXd& S2 = std::get<1>(svdBlocks.front());
            const Eigen::MatrixXd& V2 = std::get<2>(svdBlocks.front());

            // Compute A tilde and diagonalize
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(Ai - V2 * S2.asDiagonal() * V2.transpose());
            
            // Store results
            D.head(f) = std::move(eigh.eigenvalues());
            Q.block(0, 0, f, f) = std::move(eigh.eigenvectors());

            // Encode diagonal block
            diag_blocks[i] = std::make_tuple(0, 1);

        } else if (i == p-1){ // Las block 
            // Get required matrices
            const Eigen::MatrixXd& U1 = std::get<0>(svdBlocks.back());
            const Eigen::VectorXd& S1 = std::get<1>(svdBlocks.back());

            // Compute A tilde and diagonalize
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(Ai - U1 * S1.asDiagonal() * U1.transpose());

            // Store results
            D.tail(f) = std::move(eigh.eigenvalues());
            Q.block((p-1)*f, (p-1)*f, f, f) = std::move(eigh.eigenvectors());

            // Encode diagonal block
            diag_blocks[i] = std::make_tuple((p-1), 1);

        } else { // Middle blocks
            // Get required matrices
            // Need U1 and S1 from previous block
            const Eigen::MatrixXd& U1 = std::get<0>(svdBlocks[i-1]);
            const Eigen::VectorXd& S1 = std::get<1>(svdBlocks[i-1]);
            
            // Need V2 and S2 from current block
            const Eigen::VectorXd& S2 = std::get<1>(svdBlocks[i]);
            const Eigen::MatrixXd& V2 = std::get<2>(svdBlocks[i]);

            // Compute A tilde and diagonalize
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(Ai - U1 * S1.asDiagonal() * U1.transpose() - V2 * S2.asDiagonal() * V2.transpose());

            // Store results
            D.segment(i*f, f) = std::move(eigh.eigenvalues());
            Q.block(i*f, i*f, f, f) = std::move(eigh.eigenvectors());

            // Encode diagonal block
            diag_blocks[i] = std::make_tuple(i, 1);
        }

    }

    // Merge solution using binary tree algorithm
    unsigned int O;
    while(diag_blocks.size() > 1){
        std::vector<EigenTuple> tmp(diag_blocks.size()/2);

        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < diag_blocks.size()/2; ++i){ // Integer division
            EigenTuple A1 = diag_blocks[2*i];
            EigenTuple A2 = diag_blocks[2*i+1];

            // Merge blocks
            EigenTuple A3 = mergeBlocks(A1, A2);

            // Encode result
            tmp[i] = A3;
        }

        if(diag_blocks.size() & 0x1){ // Odd number of blocks
            tmp.push_back(diag_blocks.back());
        }

        diag_blocks = std::move(tmp);
    }

    // All done!
}