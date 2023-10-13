#ifndef FRECCIA_HOUSEHOLDER_DEFLATION_H
#define FRECCIA_HOUSEHOLDER_DEFLATION_H

// Standard library
#include <vector>
#include <utility>  // for std::move

// Eigen
#include <Eigen/Dense>

// Struct to store HH deflation matrices
struct HHDeflationMatrix{
    // Encode HH deflation matrices as tuples of index, multiplicit and HH vector
    using HHBlock = std::tuple<unsigned int, unsigned int, Eigen::VectorXd>;

    // Push a block onto the stack
    void push(unsigned int i, unsigned int j, Eigen::VectorXd&& u){ // Pass u as rvalue reference
        blocks.push_back(std::make_tuple(i, j, std::move(u)));
    }

    // These templated classes are used to handle view objects resulting from matrix slicing
    // Apply the HH deflation matrix to the left of a matrix
    // HH * M = (I - 2uuT / ||u||^2) * M = M - 2/||u||^2 uuT * M);
    void applyToTheLeft(Eigen::Ref<Eigen::MatrixXd> M) const {
        #pragma omp parallel for
        for(int i = 0; i < blocks.size(); i++){ 
            const HHBlock& HH = blocks[i];
            unsigned int i = std::get<0>(HH);
            unsigned int mult = std::get<1>(HH);
            const Eigen::VectorXd& u = std::get<2>(HH);
            
            // Compute indices for convenience
            unsigned int m = i;
            unsigned int n = mult;
            unsigned int k = M.cols() - n - m;
            
            // Apply reflection in O(n^2) 
            // Check if this block exists       
            if(n > 0 && m > 0){
                M.block(m, 0, n, m) = (M.block(m, 0, n, m) - 2. * u * (u.transpose() * M.block(m, 0, n, m))).eval();
            }

            // This block always exists
            M.block(m, m, n, n) = (M.block(m, m, n, n) - 2. * u * (u.transpose() * M.block(m, m, n, n))).eval();

            // Check if block exists
            if(n > 0 && k > 0){
                M.block(m, m+n, n, k) = (M.block(m, m+n, n, k) - 2. * u * (u.transpose() * M.block(m, m+n, n, k))).eval();
            }
        }
    }

    // Apply the HH deflation matrix to the right of a matrix
    // M * HH = M * (I - 2uuT / ||u||^2) = M - 2MuuT / ||u||^2;
    void applyToTheRight(Eigen::Ref<Eigen::MatrixXd> M) const {
        #pragma omp parallel for
        for(int i = 0; i < blocks.size(); i++){  
            const HHBlock& HH = blocks[i];
            unsigned int i = std::get<0>(HH);
            unsigned int mult = std::get<1>(HH);
            const Eigen::VectorXd& u = std::get<2>(HH);
            
            // Compute indices for convenience
            unsigned int m = i;
            unsigned int n = mult;
            unsigned int k = M.rows() - n - m;
            
            // Apply reflection in O(n^2) 
            // Check if this block exists       
            if(n > 0 && m > 0){
                M.block(0, m, m, n) = (M.block(0, m, m, n) - 2. * (M.block(0, m, m, n) * u) * u.transpose()).eval();
            }

            // This block always exists
            M.block(m, m, n, n) = (M.block(m, m, n, n) - 2. * (M.block(m, m, n, n) * u) * u.transpose()).eval();

            // Check if block exists
            if(n > 0 && k > 0){
                M.block(m+n, m, k, n) = (M.block(m+n, m, k, n) - 2. * (M.block(m+n, m, k, n) * u) * u.transpose()).eval();
            }
        }
    }
    
    // These templated classes are used to handle view objects resulting from matrix slicing
    // Apply the HH deflation matrix to the left of a matrix
    // HH * M = (I - 2uuT / ||u||^2) * M = M - 2/||u||^2 uuT * M);
    template <typename View>
    void applyToTheLeft(View M) const {
        // Make sure the code is compiling correctly. This function should be templated only for special Eigen objects.
        static_assert(!std::is_same<View, Eigen::Ref<Eigen::MatrixXd>>::value, "Error: generating template for Eigen::Ref object!");
        static_assert(!std::is_same<View, Eigen::MatrixXd>::value, "Error: generating template for Eigen::MatrixXd object!");

        #pragma omp parallel for
        for(int i = 0; i < blocks.size(); i++) { 
            const HHBlock& HH = blocks[i];
            unsigned int i = std::get<0>(HH);
            unsigned int mult = std::get<1>(HH);
            const Eigen::VectorXd& u = std::get<2>(HH);
            
            // Compute indices for convenience
            unsigned int m = i;
            unsigned int n = mult;
            unsigned int k = M.cols() - n - m;
            
            // Apply reflection in O(n^2) 
            // Check if this block exists       
            if(n > 0 && m > 0){
                M.block(m, 0, n, m) = (M.block(m, 0, n, m) - 2. * u * (u.transpose() * M.block(m, 0, n, m))).eval();
            }

            // This block always exists
            M.block(m, m, n, n) = (M.block(m, m, n, n) - 2. * u * (u.transpose() * M.block(m, m, n, n))).eval();

            // Check if block exists
            if(n > 0 && k > 0){
                M.block(m, m+n, n, k) = (M.block(m, m+n, n, k) - 2. * u * (u.transpose() * M.block(m, m+n, n, k))).eval();
            }
        }
    }

    // Apply the HH deflation matrix to the right of a matrix
    // M * HH = M * (I - 2uuT / ||u||^2) = M - 2MuuT / ||u||^2;
    template <typename View>
    void applyToTheRight(View M) const{
        // Make sure the code is compiling correctly. This function should be templated only for special Eigen objects.
        static_assert(!std::is_same<View, Eigen::Ref<Eigen::MatrixXd>>::value, "Error: generating template for Eigen::Ref object!");
        static_assert(!std::is_same<View, Eigen::MatrixXd>::value, "Error: generating template for Eigen::MatrixXd object!");
        
        #pragma omp parallel for
        for(int i = 0; i < blocks.size(); i++){  
            const HHBlock& HH = blocks[i];
            unsigned int i = std::get<0>(HH);
            unsigned int mult = std::get<1>(HH);
            const Eigen::VectorXd& u = std::get<2>(HH);
            
            // Compute indices for convenience
            unsigned int m = i;
            unsigned int n = mult;
            unsigned int k = M.rows() - n - m;
            
            // Apply reflection in O(n^2) 
            // Check if this block exists       
            if(n > 0 && m > 0){
                M.block(0, m, m, n) = (M.block(0, m, m, n) - 2. * (M.block(0, m, m, n) * u) * u.transpose()).eval();
            }

            // This block always exists
            M.block(m, m, n, n) = (M.block(m, m, n, n) - 2. * (M.block(m, m, n, n) * u) * u.transpose()).eval();

            // Check if block exists
            if(n > 0 && k > 0){
                M.block(m+n, m, k, n) = (M.block(m+n, m, k, n) - 2. * (M.block(m+n, m, k, n) * u) * u.transpose()).eval();
            }
        }
    }

    private:
        // Stores the indices of the HH blocks and the corresponding HH vectors
        std::vector<HHBlock> blocks;
};

#endif