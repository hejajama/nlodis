#include "integration.hpp" 
#include <string>
using namespace std;    


void Cuba(string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob) {
    // common arguments
    int ncomp=1, nvec=1, seed=0, mineval=0, last=4;
    int nregions, neval, fail;
    void *spin=NULL;
    char *statefile=NULL;
    if(method=="vegas"){
    // Vegas-specific arguments
    int nstart=1000, nincrease=500, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        cuba_config::maxeval,nstart,nincrease,nbatch,gridno,statefile,
        spin,&neval,&fail,integral,error,prob);
    }
    else if(method=="suave"){
    // Suave-specific arguments
    int nnew=1e3, nmin=2; // nnew=10e3
    double flatness=25; //25;
    Suave(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
        cuba_config::epsabs,cuba_config::verbose | last,seed,mineval,
        cuba_config::maxeval,nnew,nmin,flatness,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="divonne"){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=1*47, key2=1, key3=1, maxpass=5, ngiven=0, nextra=0;
    double border=1e-8, maxchisq=10, mindeviation=0.25;
    Divonne(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        cuba_config::maxeval,key1,key2,key3,maxpass,border,maxchisq,
        mindeviation,ngiven,ndim,NULL,nextra,NULL,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="cuhre"){
    if(ndim==1) ndim=2;
    // Cuhre-specific arguments
    int key=0;
    Cuhre(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
        cuba_config::epsabs,cuba_config::verbose | last,mineval,
        cuba_config::maxeval,key,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
}

