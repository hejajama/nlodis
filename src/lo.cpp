#include "nlodis.hpp"
#include <gsl/gsl_sf_bessel.h>

const double INTRELACC_LO = 1e-4;
using std::cout; 
using std::endl;

namespace {
struct LOIntegrationParams {
    NLODIS* nlodis;
    double Q2;
    double xbj;
    Polarization pol;
};
}

/*
 * LO Photon-proton cross section
    * Q2 [GeV^2]: photon virtuality
    * xbj: Bjorken-x
    * pol: photon polarization (T or L)
    * 
 */
double NLODIS::Photon_proton_cross_section_LO_d2b(double Q2, double xbj, Polarization pol)
{
    if (config.scheme != SubtractionScheme::UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }
    
    dipole->InitializeInterpolation( std::log(dipole->X0()/xbj) );

    LOIntegrationParams intparams = {this, Q2, xbj, pol};
    auto integrand_lo_massive = [](const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
        auto* const p = static_cast<LOIntegrationParams*>(userdata);
        *f = 0.0;

        if (*ndim != 2 or *ncomp != 1)
        {
            return 0;
        }

        // Integration variables in [0,1]
        // z is cut away from endpoints for numerical stability
        constexpr double zmin = 0; //1e-4;
        constexpr double zmax = 1.0 - zmin; //1e-4;
        constexpr double dz = zmax - zmin;

        const double z = zmin + x[0] * dz;
        const double r = p->nlodis->GetMaxR() * x[1];

        // 2*pi*r from angular integration and radial measure, and Jacobians for mappings
        const double jacobian = dz * (2.0 * M_PI * r * p->nlodis->GetMaxR());

        double res = p->nlodis->Integrand_photon_target_LO(r, z, p->xbj, p->Q2, p->pol);
        res *= jacobian;

        if (std::isfinite(res))
        {
            *f = res;
        }

        return 0;
    };

    double result = 0.0;
    double abserr = 0.0;
    double prob = 0.0;
    Cuba("cuhre", 2, integrand_lo_massive, &intparams, &result, &abserr, &prob, cuba_config);



    result *= 4.0*Constants::AlphaEM*Constants::NC/SQR(2.0*M_PI); // Include prefactors as needed
    // Note: 1/(2pi)^2 because 1708.07328 (1) includes 2pi to transverse coordiante measures
    
    // \int d^2b not included in here

    return result;


}


/*
 * Integrand for the LO cross section
 * r: dipole size [GeV^-1]
 * z: longitudinal momentum fraction of the quark
 * x: Bjorken-x
 * Q2: photon virtuality [GeV^2]
 * pol: photon polarization (T or L)
 * 
 * Note: In the UNSUB scheme, expect user to evaluate the dipole amplitude at the initial x
 * 
 * This is |\Psi|^2 2N(r), no Jacobian of r included
 * Factor 2 from the optical theorem is included in thje overall normalization factor
 * Reference: 1708.07328 (1-3), this is K_L, K_T
 * 
*/
double NLODIS::Integrand_photon_target_LO(double r, double z, double x, double Q2, Polarization pol )
{
    if (config.scheme != SubtractionScheme::UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }

    double res=0;
    
    
    for (const auto& quark : quarks) {
        double eps = sqrt( Q2*z*(1.0-z) + SQR(quark.mass) );
        if (r*eps < 1e-7 or r*eps > 5e2) { // GSL overflow/underflow
            continue;
        }

        if (pol == Polarization::T)
        {
            res += SQR(quark.Charge())*((1.0-2.0*z+2.0*SQR(z))*SQR(eps)*SQR(gsl_sf_bessel_K1(r*eps)) 
                + SQR( quark.mass*gsl_sf_bessel_K0( r*eps ) ) );
        }
        else if (pol == Polarization::L)
        { 
            res += SQR(quark.Charge()) * 4.0 * Q2 * SQR(z) * SQR(1.0 - z) * SQR(gsl_sf_bessel_K0(r*eps));
        }        
        else
        {
            throw std::runtime_error("Unknown polarization in Integrand_LO.");
        }
    }

    double evolution_rapidity = std::log(dipole->X0() / x); 
    if (evolution_rapidity < 0)
    {
        cout << "Warning: evolution rapidity " << evolution_rapidity << "< 0 in Integrand_LO. Setting to 0." << endl;
        evolution_rapidity = 0;
    }
    
    res *= dipole->DipoleAmplitude(r, evolution_rapidity); // Dipole amplitude at x

    
    return res;
}
