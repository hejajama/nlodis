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
#include "integration.hpp"

using namespace std;
const string cubamethod = "cuhre";



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


/*
 * qqg-target scattering amplitude
 * Ref. https://arxiv.org/pdf/2211.03504 (4)
 * In that notation, this is 1-S_{012}
 * 
 * x01, x02, x21: dipole sizes in GeV^-1
 * Y: evolution rapidity
 **/
double NLODIS::TripoleAmplitude(double x01, double x02, double x21, double Y)
{
    double S01 = 1-dipole.DipoleAmplitude(x01, Y);
    double S02 = 1-dipole.DipoleAmplitude(x02, Y);
    double S12 = 1-dipole.DipoleAmplitude(x21, Y);

    if (nc_scheme == LargeNC)
    {
        return 1.0 - S02*S12;
    }
    else if (nc_scheme == FiniteNC)
    {
        return NC/(2.0*CF)*(S02*S12 - 1./SQR(NC)*S01);
    }
    else
    {
        throw std::runtime_error("NLODIS::TripoleAmplitude: unknown NC scheme");
    }
}

/*
 * Evolution rapidity (in the qqg contribution)
 *
 * Ref https://arxiv.org/pdf/2007.01645 eq (19)
 * xbj: Bjorken-x
 * Q2: photon virtuality in GeV^2
 * z2: gluon longitudinal momentum fraction
 */
double NLODIS::EvolutionRapidity(double xbj, double Q2, double z2)
{
    double W2 = Q2 / xbj;
    return std::log(W2*z2/Q0sqr);
}

/*
 * Running coupling scale depending on the RC scheme used
 */
double NLODIS::RunningCouplinScale(double x01, double x02, double x21)
{
    if (rc_scheme == SMALLEST)
    {
        return std::min({x01, x02, x21});
    }
    else if (rc_scheme == PARENT)
    {
        return x01;
    }
    else
    {
        throw std::runtime_error("NLODIS::RunningCouplinScale: unknown running coupling scheme");
    }
}

/*
 * Photon-proton cross section [GeV^-2]
 *
 * To get the cross section, this has to be integrated over d^2r and multiplied by sigma_0 [in GeV^-2]
 * i.e. we replace 2\int d^2b -> sigma_0
 * 
 * */
double NLODIS::Photon_proton_cross_section(double Q2, double xbj, Polarization pol)
{
    if (scheme != UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }

    if (order==LO)
    {
        return Photon_proton_cross_section_LO(Q2, xbj, pol);
    }
    
    // NLO calculation

    double sigma_LO = Photon_proton_cross_section_LO(Q2, dipole.X0(), pol);
    double sigma_dip = Sigma_dip(Q2, xbj, pol);
    double  sigma_qg = 0; // TODO

    return sigma_LO + sigma_dip + sigma_qg;
}


/*
 * \sigma_dip
 * qq part of the NLO cross section
 * L polarization: https://arxiv.org/pdf/2103.14549 (166)
 * T polarization:
 */
double NLODIS::Sigma_dip(double Q2, double xbj, Polarization pol)
{
    double result=0;
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^2 below (from d^2x_{01} d^2b)
    double fac=4.0*NC*ALPHA_EM/SQR(2.0*M_PI); 
    IntegrationParams intparams;
    intparams.nlodis=this;
    intparams.Q2=Q2;
    intparams.xbj=xbj;
    intparams.pol=pol;

    for (const auto& quark : quarks) {
        intparams.quark=quark;
        if (pol == L)
        {
            // 1st line
            double I, Ierr, Iprob;
            intparams.contribution="Omega_L_const";
            Cuba(cubamethod, 2, integrand_ILdip_massive, &intparams, &I, &Ierr, &Iprob);
            result += SQR(quark.charge) * I;

            // 2nd line of 2103.14549 (166)
            intparams.contribution="ab";
            double Iab,Iaberr,Iabprob;
            Cuba(cubamethod, 3, integrand_ILdip_massive, &intparams, &Iab, &Iaberr, &Iabprob);
            intparams.contribution="cd";
            double Icd,Icderr,Icdprob;
            Cuba(cubamethod, 4, integrand_ILdip_massive, &intparams, &Icd, &Icderr, &Icdprob);
            result += SQR(quark.charge) * (Iab + Icd);
        }

    }

    // We have facotorized out \int d^2 b - note the normalization convention!
    // Define sigma_0 = 2 \int d^2 b
    
    // Correspondingly I need to include 1/2
    result *= 1./2.;

        
    // 2pi from overall angular integral
    return fac*result*2.0*M_PI;
}

/*
 * \sigma_qg
 * qg part of the NLO cross section
 * 
 * Longitudinal reference: (167) but instead of q^+, k^+ we integrate over z_i
 * Explicit expressoin is docs/NLO_DIS_cross_section_with_massive_quarks.pdf (13) 
*/
double NLODIS::Sigma_qg(double Q2, double xbj, Polarization pol)
{
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^3 below (from d^2x_{01} d^2x_{02} d^2b)
    double fac=4.0*NC*ALPHA_EM/std::pow(2.0*M_PI,3.0);
    IntegrationParams intparams;
    intparams.nlodis=this;
    intparams.Q2=Q2;
    intparams.xbj=xbj;
    intparams.pol=pol;

    double result=0;
    for (const auto& quark : quarks) {
        intparams.quark=quark;
        if (pol == L)
        {
            // Note (21): this contribution is split into 3 parts I_1, I_2 and I_3
            double I2,I2err,I2prob;
            intparams.contribution="I2"; //"fast" means that u integration is done analytically(?)
            Cuba(cubamethod, 6, integrand_ILqgunsub_massive, &intparams, &I2, &I2err, &I2prob);
            result += SQR(quark.charge) * I2;
            
        }
    }


    // We have facotorized out (not performed) \int d^2 b - note the normalization convention!
    // Define sigma_0 = 2 \int d^2 b
    
    // Correspondingly I need to include 1/2 to this result
    result *= 1./2.;

    // 2pi: overall integral over one angle
    return  fac*2.0*M_PI*result;

}

/*
 * LO Photon-proton cross section
    * Q2 [GeV^2]: photon virtuality
    * xbj: Bjorken-x
    * pol: photon polarization (T or L)
    * 
 */
double NLODIS::Photon_proton_cross_section_LO(double Q2, double xbj, Polarization pol)
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
    IntegrationParams params_struct = {this, Q2, xbj, 0, pol, nullptr};

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
 
 /*
  * Coordinate space coupling
  * TODO: implement flavor thresholds and C^2 scale
  */
 double NLODIS::Alphas(double r)
 {
    const double LambdaQCD = 0.241; // GeV
    const double b0 = (11.0*NC - 2.0*quarks.size())/(12.0*M_PI);
    double mu2 = 4.0/(r*r) + LambdaQCD*LambdaQCD;
    double as = 1.0/(b0*log(mu2/(LambdaQCD*LambdaQCD)));
    if (as > 0.7)
    {
        as = 0.7; // Freeze coupling
    }
    return as;
 }