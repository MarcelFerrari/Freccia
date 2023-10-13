#ifndef FRECCIA_ARRANGEMENT_MATRIX_H
#define FRECCIA_ARRANGEMENT_MATRIX_H

// Standard library
#include <vector>
#include <algorithm>  // for std::sort and std::copy_n
#include <numeric>    // for std::iota
#include <cassert>    // for assert
#include <cmath>      // for std::abs

// Eigen
#include <Eigen/Dense>  // For Eigen::ArrayXd


// Base class for implicit permutation and partition matrices
class ArrangementBase {
    protected: // Need to be accessed by children so not private
    std::vector<unsigned int> _perm;
    std::vector<unsigned int> _invperm;
    bool invert = true; // Flag to indicate if inverse permutation needs to be recomputed

    public:

    ArrangementBase(unsigned int N) : _perm(N), _invperm(N) {
        std::iota(_perm.begin(), _perm.end(), 0);
    }

    // Get size
    unsigned int size(void){
        return _perm.size();
    }

    // Get permutation vector
    std::vector<unsigned int> & perm(void){
        return _perm;
    }
};

class Permutation : public ArrangementBase {
public:
    Permutation(unsigned int N) : ArrangementBase(N) {}

    // Argsort for eigen arrays. Default is decreasing order.
    void argsort(Eigen::ArrayXd& vec, bool reverse = true){
        unsigned int N = vec.size();
        
        // Check if vector is of correct size
        assert(N == _perm.size());

        // Sort perm
        auto comp = [&](unsigned int i, unsigned int j) {
            return reverse ? vec(i) > vec(j) : vec(i) < vec(j);
        };

        std::sort(_perm.begin(), _perm.end(), comp);

        // Compute and store inverse permutation directly after sorting
        computeInvperm();
    }

    // Get inverse permutation vector
    const std::vector<unsigned int>& invperm() const {        
        return _invperm;
    }

private:
    // Compute the inverse permutation and store it in _invperm
    void computeInvperm() {
        for(unsigned int i = 0; i < _perm.size(); ++i){
            _invperm[_perm[i]] = i;
        }
    }
};


class Partition: public ArrangementBase {
private:
    int nnz = 0; // Number of non-zero elements
    std::vector<unsigned int> _zero;
    std::vector<unsigned int> _nnzero;

public:
    Partition(unsigned int N) : ArrangementBase(N) {}

    // Partition zero and non-zero elements based on threshold value val
    // This maintains the order of the nnzero elements
    void partition(Eigen::ArrayXd& vec, double val = 1e-15) {
        // Ensure correct vector size
        assert(vec.size() == _perm.size());
        
        // Start and end of vector
        unsigned int head = 0;
        unsigned int tail = _perm.size() - 1;

        // Create space to store temporary partition result
        // Need this to preserve sorting order
        std::vector<unsigned int> tmp(_perm.size());
        
        for(unsigned int i = 0; i < _perm.size(); ++i) {
            if (std::abs(vec(i)) < val) {
                tmp[tail--] = _perm[i];
            } else {
                tmp[head++] = _perm[i];
            }
        }

        std::swap(_perm, tmp);
        
        nnz = head;

        // Store nnzero and zero indices directly after partitioning
        _nnzero.assign(_perm.begin(), _perm.begin() + nnz);
        _zero.assign(_perm.begin() + nnz, _perm.end());

        invert = true; // Need to recompute inverse
    }
    
    // Get non-zero indices
    const std::vector<unsigned int>& nnzero() const {
        return _nnzero;
    }

    // Get zero indices
    const std::vector<unsigned int>& zero() const {
        return _zero;
    }
};


#endif