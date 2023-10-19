#ifndef FRECCIA_ROOT_FINDING_H
#define FRECCIA_ROOT_FINDING_H
// Fast root finding routines to speed up the computation of the secular equation
// This relies on the boost c++ library
// These wrappers are required to deal with the poles of the secular equation

#include <boost/math/tools/roots.hpp>
#include "Freccia/options/options.hpp"

namespace Freccia::RootFinding {
    
    template<typename FUNC, typename OPTIONS>
    double toms748(double left, double right, const FUNC& F, const OPTIONS& opt){
        // Adjust the left and/or right bounds
        const double rel_shift_tol = 1e-7;
        const double abs_shift_tol = 1e-9;

        if(!std::isfinite(F(left))){ // Left bound is on a pole
            const double shift = left >= 1. ? std::max(rel_shift_tol * std::abs(left), abs_shift_tol) : abs_shift_tol;
            left += shift;
        }

        if(!std::isfinite(F(right))){ // Right bound is on a pole
            const double shift = right >= 1. ? std::max(rel_shift_tol * std::abs(right), abs_shift_tol) : abs_shift_tol;
            right -= shift;
        }

        // Create variable for max iterations
        std::uintmax_t max_iter = opt.BISECT_MAX_ITER;

        // Solve using Toms748 from boost library
        // This method throws so it needs to be wrapped in a try-catch block
        try {
            std::pair<double, double> result = boost::math::tools::toms748_solve(F, left, right, boost::math::tools::eps_tolerance<double>(), max_iter);
            return result.second;
        } catch (const std::exception& e) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    };
}
#endif