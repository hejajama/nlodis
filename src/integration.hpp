#pragma once

#include <string>
#include <cuba.h>

struct CubaConfig
{
    const int verbose=0;
    int maxeval=5e5;
    const double epsrel=1e-4;
    const double epsabs=0;
    std::string method="suave";
};

// Wrapper that allows user to specify the Cuba method to use
void Cuba(std::string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob, CubaConfig config = CubaConfig()); 
