#ifndef FRECCIA_ENV_H
#define FRECCIA_ENV_H

namespace Freccia::Env {
    void read_env(const char* env_var, double & var);

    void read_env(const char* env_var, int & var);

    void read_env(const char* env_var, unsigned int & var);
}

#endif