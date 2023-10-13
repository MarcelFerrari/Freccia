#include <algorithm> // For std::max
#include <cmath>     // For std::abs
#include <cstdlib>  // for std::getenv
#include <string>   // for std::stoi
#include <iostream> // for std::cerr
#include <exception> // for std::exception

#include "Freccia/utils/utils.hpp"

void read_env(const char* env_var, double & var) {
    char* env = std::getenv(env_var);
    if(env == NULL) {
        return;
    } else {
        try {
            var = std::stod(env);
        } catch(const std::exception& e) {
            std::cerr << "Warn: invalid value for variable \"" << env_var << "\"."<< std::endl;
            return;
        }
    }
}

void read_env(const char* env_var, int & var) {
    char* env = std::getenv(env_var);
    if(env == NULL) {
        return;
    } else {
        try {
            var = std::stoi(env);
        } catch(const std::exception& e) {
            std::cerr << "Warn: invalid value for variable \"" << env_var << "\"." << std::endl;
            return;
        }
    }
}

void read_env(const char* env_var, unsigned int & var) {
    char* env = std::getenv(env_var);
    if(env == NULL) {
        return;
    } else {
        try {
            // Convert string to unsigned long and then to unsigned int
            var = static_cast<unsigned int>(std::stoul(env));
        } catch(const std::exception& e) {
            std::cerr << "Warn: invalid value for variable \"" << env_var << "\"." << std::endl;
        }
    }
}

bool isclose(double x, double y, double atol, double rtol) {
    double diff = std::abs(x - y);
    if(diff < atol){
        return true;
    }
    double max = std::max(std::abs(x), std::abs(y));
    if(max >= 1.){ 
        return diff <= rtol * max;
    } else if (max > 0.){ 
        return diff / max <= rtol;
    } else { 
        return true;
    }
}
