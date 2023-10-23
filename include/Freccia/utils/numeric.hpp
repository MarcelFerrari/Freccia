#ifndef FRECCIA_NUMERIC_H
#define FRECCIA_NUMERIC_H

namespace Freccia::Numeric {

template<typename T>
inline int sign(T x){
    return (int) (x > 0.0);
}

template<typename T>
bool isclose(T x, T y, T atol, T rtol) {
    T diff = std::abs(x - y);
    if(diff < atol){
        return true;
    }
    T max = std::max(std::abs(x), std::abs(y));
    if(max >= 1.){ 
        return diff <= rtol * max;
    } else if (max > 0.){ 
        return diff / max <= rtol;
    } else { 
        return true;
    }
}

}

#endif //FRECCIA_NUMERIC_H