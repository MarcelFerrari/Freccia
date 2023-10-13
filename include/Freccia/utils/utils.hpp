#ifndef FRECCIA_UTILS_H
#define FRECCIA_UTILS_H

void read_env(const char* env_var, double & var);

void read_env(const char* env_var, int & var);

void read_env(const char* env_var, unsigned int & var);

bool isclose(double x, double y, double atol, double rtol);

inline int sign(double x){
    return (int) (x > 0.0);
}

#endif //FRECCIA_UTILS_H
