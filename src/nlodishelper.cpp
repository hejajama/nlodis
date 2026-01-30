/* 
 * NLO DIS cross section
 *
 * This file contains some internal functions needed to evaluate the NLO DIS
 * impact factors
 *  
 * Functions are defined such that they are suitable for Cuba 
*/

#include "nlodis.hpp"
#include "qcd.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>


/*
* Longitudinal photon, dipole contribution with massive quarks
* This term corresponds to second line of 2103.14549 (166)
* i.e. 4 Q^2 z^2(1-z)^2 alpha_s*Nc/pi*K_0(r*eps)*I_{v, part (c)+(d)}
* I_{v, part (c)+(d)} is given in (114)
 */

double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf, double xi, double x);

/*
 * Integrand for Cuba
 * Integration variables are
 * x[0] = z1
 * x[1] = [0,1] mapped to r = x[1]*maxr
 * x[2] = xi (integration variable xi in (114))
 * x[3] = x (integration variable x in (114))
 */
int integrand_ILdip_massive_Icd(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
     throw std::invalid_argument("integrand_ILdip_massive_Icd: ndim must equal 4");
    
    auto* p = static_cast<IntegrationParams*>(userdata);
    double Q2=p->Q2;
    double xbj=p->xbj;
    double mf=p->quark.mass;

    double z1=x[0];
    double x01=p->nlodis->GetMaxR()*x[1];
    double x01sq=SQR(x01);
    double xi=x[2]; 
    double intx=x[3]; 

    double alphabar=p->nlodis->Alphas(x01)*CF/M_PI;

    // TODO: add more user control for evolution rapidity
    double evolution_rapidity = std::log(1/xbj); //Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double dipole = p->nlodis->GetDipole().DipoleAmplitude(x01,evolution_rapidity);
    double res;

    res = dipole*(ILdip_massive_Icd(Q2,z1,x01,mf,xi,intx))*x01*alphabar;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


/*
 * https://arxiv.org/pdf/2103.14549 (114) 
*/
double ILdip_massive_Icd(double Q2, double z1, double r, double mf, double xi, double x) {
    // Note that this now uses the earlier parametrization (xi,x) instead of (chi,u). This should be better for numerics.
    double front_factor = 4.0*Q2*SQR(z1*(1.0-z1));

    double kappa_z = sqrt( z1*(1.0-z1)*Q2 + SQR(mf) );
    double bessel_inner_fun = kappa_z * r;
    double Icd_integrand = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        Icd_integrand = 0;
    }else{
        // 2103.14549 (89,97)
        double CLm1 = SQR(z1)*(1.0-xi)/(1.0-z1) * ( -SQR(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+z1*xi/(1.0-z1)) ) / (x*(1.0-xi)+xi/(1.0-z1)) ); // C^L_m(z,x,xi)
        double CLm2 = SQR(1.0-z1)*(1.0-xi)/(z1) * ( -SQR(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+(1.0-z1)*xi/(z1)) ) / (x*(1.0-xi)+xi/(z1)) );; // C^L_m(1-z,x,xi)
        double kappa1 = xi*SQR(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) ) * ( xi*(1.0-x) + x*(1.0-z1*(1.0-xi)/(1.0-z1)) ); // kappa(z,x,xi)
        double kappa2 = xi*SQR(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) ) * ( xi*(1.0-x) + x*(1.0-(1.0-z1)*(1.0-xi)/(z1)) );; // kappa(1-z,x,xi)
        Icd_integrand = gsl_sf_bessel_K0( bessel_inner_fun ) * SQR(mf)*
                     ( (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( r*sqrt( SQR(kappa_z)/(1.0-x) + kappa1 ) ) ) *
                        CLm1 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) * ( x/(1.0-x)*SQR(kappa_z) + kappa1 ) ) +
                        (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( r*sqrt( SQR(kappa_z)/(1.0-x) + kappa2 ) ) ) *
                        CLm2 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) * ( x/(1.0-x)*SQR(kappa_z) + kappa2 ) ));
    }   

    double dip_res = front_factor * Icd_integrand;
    return dip_res;
}