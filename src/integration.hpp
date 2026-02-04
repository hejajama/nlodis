#ifndef _INTEGRATION_CPP_
#define _INTEGRATION_CPP_

#include <string>
#include <cuba.h>

namespace cuba_config{
    static const int verbose=0;
    static const int maxeval=1e5;
    static const double epsrel=1e-3;
    static const double epsabs=0;
}

// Wrapper that allows user to specify the Cuba method to use
void Cuba(std::string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob);


#endif 
