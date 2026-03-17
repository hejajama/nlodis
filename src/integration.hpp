#pragma once

#include <string>
#include <cuba.h>

struct CubaConfig
{
    const int verbose=0;
    int maxeval=2e6;
    double epsrel=0.001;
    const double epsabs=0;
    double warning_relaccuracy_threshold = 1e-1; // If the relative accuracy of the integral is worse than this, print a warning. This is not a hard failure, just a warning that the integration may not be accurate enough for some purposes.
    std::string method="vegas"; // Recommended integration method

};

/**
 * @brief Performs numerical integration using the Cuba library with the specified method.
 *
 * This function dispatches to one of the Cuba integration routines (Vegas, Suave, Divonne, or Cuhre)
 * based on the provided method string. It configures method-specific parameters and handles integration
 * results, including warning if the achieved accuracy is below a user-defined threshold.
 * 
 * @param method The integration method to use ("vegas", "suave", "divonne", or "cuhre").
 * @param ndim The number of dimensions for the integration.
 * @param integrand The integrand function to be evaluated.
 * @param userdata Pointer to user data passed to the integrand.
 * @param integral Pointer to the variable where the computed integral will be stored.
 * @param error Pointer to the variable where the estimated error will be stored.
 * @param prob Pointer to the variable to store the probability that the error estimate is reliable.
 * @param config Configuration parameters for the integration (maxeval, epsrel, epsabs, verbose, etc.).
 *
 * @note If the integration fails and the relative accuracy is below the warning threshold,
 *       a warning message is printed suggesting to increase the number of Monte Carlo points.
 */
void Cuba(std::string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob, CubaConfig config = CubaConfig()); 
