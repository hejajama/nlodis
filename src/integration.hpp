#pragma once

#include <string>
#include <cuba.h>

struct CubaConfig
{
    const int verbose=0;
    int maxeval=1e6;
    double epsrel=0.001;
    const double epsabs=0;
    double warning_relaccuracy_threshold = 1e-1; // If the relative accuracy of the integral is worse than this, print a warning. This is not a hard failure, just a warning that the integration may not be accurate enough for some purposes.
    std::string method="vegas"; // Recommended integration method
};

// Wrapper that allows user to specify the Cuba method to use
void Cuba(std::string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob, CubaConfig config = CubaConfig()); 
