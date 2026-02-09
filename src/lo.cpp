#include "nlodis.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

const double INTRELACC_LO = 1e-4;
using std::cout; 
using std::endl;

/*
 * LO Photon-proton cross section
    * Q2 [GeV^2]: photon virtuality
    * xbj: Bjorken-x
    * pol: photon polarization (T or L)
    * 
 */
double NLODIS::Photon_proton_cross_section_LO_d2b(double Q2, double xbj, Polarization pol)
{
    if (scheme != UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }
    
    dipole->InitializeInterpolation( std::log(dipole->X0()/xbj) );

    // Simple approach with nested 1D integrations
    gsl_function F;
    double result = 0.0;
    double abserr = 0.0;
    size_t neval = 0;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    // For r: integrate from 0 to infinity
    // For z: integrate from 0 to 1

    struct LOParams {
        NLODIS* nlodis;
        double Q2;
        double xbj;
        double r;
        Polarization pol;
    } params_struct = {this, Q2, xbj, 0.0, pol};

    // r integrand (outer)
    F.function = [](double r, void* params) {
        auto* p = static_cast<LOParams*>(params);
        p->r = r;

        gsl_function F_z;
        // z integrand (inner)
        F_z.function = [](double z, void* z_params) {
            auto* int_params = static_cast<LOParams*>(z_params);
            return 2.0*M_PI*int_params->r*int_params->nlodis->Integrand_photon_target_LO(int_params->r, z, int_params->xbj, int_params->Q2, int_params->pol);
        };
        F_z.params = p;

        double z_result = 0.0, z_err = 0.0;
        gsl_integration_workspace *w_z = gsl_integration_workspace_alloc(1000);
        gsl_integration_qag(&F_z, 1e-4, 1.0-1e-4, 0, INTRELACC_LO, 1000, GSL_INTEG_GAUSS21,
            w_z, &z_result, &z_err);
        gsl_integration_workspace_free(w_z);

        //cout << r << " " << z_result << endl;
        return z_result; // Note: Jacobian is in the innermost function
    };
    F.params = &params_struct;

    gsl_integration_qagiu(&F, 0.0, 0, INTRELACC_LO, 1000, w, &result, &abserr);

    gsl_integration_workspace_free(w);



    result *= 4.0*ALPHA_EM*NC/SQR(2.0*M_PI); // Include prefactors as needed
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
 * This is |\Psi|^2 N(r), no Jacobian of r included
 * Reference: 1708.07328 (1-3), this is K_L, K_T
 * 
*/
double NLODIS::Integrand_photon_target_LO(double r, double z, double x, double Q2, Polarization pol )
{
    if (scheme != UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }

    double res=0;
    
    
    for (const auto& quark : quarks) {
        double eps = sqrt( Q2*z*(1.0-z) + SQR(quark.mass) );
        if (r*eps < 1e-7 or r*eps > 5e2) { // GSL overflow/underflow
            continue;
        }

        if (pol == T)
        {
            res += SQR(quark.charge)*((1.0-2.0*z+2.0*SQR(z))*SQR(eps)*SQR(gsl_sf_bessel_K1(r*eps)) 
                + SQR( quark.mass*gsl_sf_bessel_K0( r*eps ) ) );
        }
        else if (pol == L)
        { 
            res += SQR(quark.charge) * 4.0 * Q2 * SQR(z) * SQR(1.0 - z) * SQR(gsl_sf_bessel_K0(r*eps));
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
