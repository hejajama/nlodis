#include "integration.hpp" 
#include <string>
#include <iostream>
using namespace std;    


int integrand_dip_massive(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata);

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
void Cuba(string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob,
    CubaConfig config) {

    if (method != "vegas" and method != "suave" and method != "divonne" and method != "cuhre")
    {
        throw std::runtime_error("Invalid Cuba integration method: " + method);
    }

    int maxeval = config.maxeval;
    if (integrand == integrand_dip_massive)
    {
        maxeval *= 10; // Increase maxeval for the dipole integrand, which is more challenging to integrate accurately and faster to evaluate than the qqg integrand. 
        // This is a heuristic choice based on testing, and can be adjusted if needed.
        // Note that this is the same integrand for both polarizations
    }

    // common arguments
    int ncomp=1, nvec=1, seed=0, mineval=maxeval/10;
    int nregions, neval, fail;
    void *spin=NULL;
    char *statefile=NULL;
    if(method=="vegas"){
    // Vegas-specific arguments
    //int nstart=mineval, nincrease=config.maxeval/20, nbatch=100, gridno=0;
    int nstart=mineval/10, nincrease=nstart/2, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose,seed,mineval,
        maxeval,nstart,nincrease,nbatch,gridno,statefile,
        spin,&neval,&fail,integral,error,prob);
    }
    else if(method=="suave"){
    // Suave-specific arguments
    int nnew=1e3, nmin=10; // nnew=10e3
    double flatness = 25;
    int last=4;
    Suave(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose | last,seed,mineval,
        maxeval,nnew,nmin,flatness,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="divonne"){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=1*47, key2=1, key3=1, maxpass=5, ngiven=0, nextra=0;
    double border=1e-8, maxchisq=10, mindeviation=0.25;
    Divonne(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose,seed,mineval,
        maxeval,key1,key2,key3,maxpass,border,maxchisq,
        mindeviation,ngiven,ndim,NULL,nextra,NULL,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="cuhre"){
    if(ndim==1) ndim=2;
    // Cuhre-specific arguments
    int key=9;
    int last=4;
    Cuhre(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose | last,mineval,
        config.maxeval,key,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }

    if (fail > 0 and std::abs(*error / *integral) > config.warning_relaccuracy_threshold) {
        cout << "# Warning: Cuba (" << method << ") integration (dimension " << ndim << ") accuracy " << std::abs(*error / *integral)*100 << "%,  result " << *integral << " +/- " << *error << ". It is recommended to use more MC integration points" <<  endl;
    }
}

