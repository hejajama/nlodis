#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_errno.h>

#include "nlodis.hpp"
#include "qcd.hpp"

using namespace std;

inline double SQR(double x) { return x*x; }

const double INTRELACC = 1e-3;

/*
 * Structure function F2 
 *
 * Q2 [GeV^2]: photon virtuality
 * xbj: Bjorken-x
 * 
 */
double NLODIS::F2(double Q2, double xbj)
{
    double sigmaT = Photon_proton_cross_section(Q2, xbj, T);
    double sigmaL = Photon_proton_cross_section(Q2, xbj, L);

    return Q2 / (4.0 * M_PI * M_PI * ALPHA_EM) * (sigmaT + sigmaL);
}





double NLODIS::Photon_proton_cross_section(double Q2, double xbj, Polarization pol)
{
    if (scheme != UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }
    

    // Simple approach with nested 1D integrations
    gsl_function F;
    double result = 0.0;
    double abserr = 0.0;
    size_t neval = 0;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    // For r: integrate from 0 to some large cutoff
    // For z: integrate from 0 to 1
    
    // Create a struct to hold the integration parameters (must survive the integration)
    struct IntegrationParams {
        NLODIS* nlodis;
        double Q2;
        double xbj;
        double z;
        Polarization pol;
        gsl_integration_workspace* w_r;
    } params_struct = {this, Q2, xbj, 0, pol, nullptr};

    // z integrand
    F.function = [](double z, void* params) {
        auto* p = static_cast<IntegrationParams*>(params);
        p->z = z;
        gsl_function F_r;
        // r integrand
        F_r.function = [](double r, void* r_params) {
            auto* int_params = static_cast<IntegrationParams*>(r_params);
            // Below add Jacobian r and 2pi from angular integral
            return 2.0*M_PI*r*int_params->nlodis->Integrand_photon_target_LO(r, int_params->z, int_params->xbj, int_params->Q2, int_params->pol);
        };
        F_r.params = p;
        
        double r_result = 0.0, r_err = 0.0;
        gsl_integration_workspace *w_r = gsl_integration_workspace_alloc(1000);
        gsl_integration_qagiu(&F_r, 0.0, 0, INTRELACC, 1000, w_r, &r_result, &r_err);
        gsl_integration_workspace_free(w_r);
        return r_result;
    };
    F.params = &params_struct;

    gsl_integration_qag(&F, 1e-3, 1.0-1e-3, 0, INTRELACC, 1000, GSL_INTEG_GAUSS21,
        w, &result, &abserr);

    gsl_integration_workspace_free(w);

    result *= 4.0*ALPHA_EM*NC/SQR(2.0*M_PI); // Include prefactors as needed
    // Note: 1/(2pi)^2 because 1708.07328 (1) includes 2pi to transverse coordiante measures
    result /= 2.0; // Follow convention that 2\int d^2b is replaced by sigma_0
    // which is not included here
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
 * To get the cross section, this has to be integrated over d^2r and multiplied by sigma_0 [in GeV^-2]
 * i.e. we replace 2\int d^2b -> sigma_0
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

    double evolution_rapidity = std::log(dipole.X0() / x); 
    if (evolution_rapidity < 0)
    {
        cout << "Warning: evolution rapidity " << evolution_rapidity << "< 0 in Integrand_LO. Setting to 0." << endl;
        evolution_rapidity = 0;
    }
    
    res *= dipole.DipoleAmplitude(r, evolution_rapidity); // Dipole amplitude at x

    
    return res;
}


 NLODIS::NLODIS(string bkdata) : dipole(bkdata)
 {
    // Initialize quark flavors and masses
    Quark u; u.type = Quark::U; u.mass = 0.14; u.charge = 2.0/3.0;
    Quark d; d.type = Quark::D; d.mass = 0.14; d.charge = -1.0/3.0;
    Quark s; s.type = Quark::S; s.mass = 0.14;  s.charge = -1.0/3.0;
    //Quark c; c.type = Quark::C; c.mass = 1.27;   c.charge = 2.0/3.0;
    //Quark b; b.type = Quark::B; b.mass = 4.18;   b.charge = -1.0/3.0;

    quarks.push_back(u);
    quarks.push_back(d);
    quarks.push_back(s);
    //quarks.push_back(c);
    //quarks.push_back(b);

    dipole.SetOutOfRangeErrors(false);

 }