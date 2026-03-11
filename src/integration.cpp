#include "integration.hpp" 
#include <string>
#include <iostream>
using namespace std;    


void Cuba(string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob,
    CubaConfig config) {
    // common arguments
    int ncomp=1, nvec=1, seed=0, mineval=config.maxeval/100;
    int nregions, neval, fail;
    void *spin=NULL;
    char *statefile=NULL;
    if(method=="vegas"){
    // Vegas-specific arguments
    //int nstart=mineval, nincrease=config.maxeval/20, nbatch=100, gridno=0;
    int nstart=1000, nincrease=500, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose,seed,mineval,
        config.maxeval,nstart,nincrease,nbatch,gridno,statefile,
        spin,&neval,&fail,integral,error,prob);
    }
    else if(method=="suave"){
    // Suave-specific arguments
    int nnew=1e3, nmin=10; // nnew=10e3
    double flatness = 25;
    int last=4;
    Suave(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose | last,seed,mineval,
        config.maxeval,nnew,nmin,flatness,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="divonne"){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=1*47, key2=1, key3=1, maxpass=5, ngiven=0, nextra=0;
    double border=1e-8, maxchisq=10, mindeviation=0.25;
    Divonne(ndim,ncomp,integrand,userdata,nvec,config.epsrel,
        config.epsabs,config.verbose,seed,mineval,
        config.maxeval,key1,key2,key3,maxpass,border,maxchisq,
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

