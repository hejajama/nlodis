/*
 * NLO DIS code
 * Based on code originally written by H. Hänninen
 * Modified by J. Penttala, H. Mäntysaari, C. Casuga
 */

#ifndef _NLODIS_HPP_
#define _NLODIS_HPP_

#include <amplitudelib.hpp> // Dipole amplitude
#include "qcd.hpp"
#include <vector>
#include <string>
#include <gsl/gsl_integration.h>

enum Polarization
{
    T,
    L
};

enum Scheme
{
    UNSUB
};

enum Order
{
    LO,
    NLO
};



class NLODIS
{
    public:
        
        NLODIS(std::string bkdata);

        double F2(double Q2, double xbj);
        double Photon_proton_cross_section(double Q2, double xbj, Polarization pol);

         double Photon_proton_cross_section_LO(double Q2, double xbj, Polarization pol);

         /// NLO caluclation ingredients 
        double Sigma_dip(double Q2, double xbj, Polarization pol);

         ////


    
        void SetOrder(Order o) { order = o; }
        double GetMaxR() { return maxr; }
        
        double Alphas(double r);
        AmplitudeLib& GetDipole() { return dipole; }

    private:
    
        /* Integrand for LO cross section
        * |\Psi|^2 N(r): integrand for the LO cross section, to be integrated over r and z
        * Note: Jacobian of r is _NOT_ included here
        * The dipole-target cross section is sigma_0 \int d^2r dz Integrand_photon_target_LO
        * The factor \sigma_0 = 2 \int d^2 b is NOT included here
        */
        double Integrand_photon_target_LO(double r, double z, double x, double Q2, Polarization pol );
    
        AmplitudeLib dipole;
        Scheme scheme = UNSUB;
        std::vector<Quark> quarks;
        Order order = LO;
        double maxr = 99;
      
};

// Data structure used to store integration parameters
struct IntegrationParams {
        NLODIS* nlodis;
        double Q2;
        double xbj;
        double z;
        Polarization pol;
        gsl_integration_workspace* w_r;
        Quark quark;
        std::string contribution;
    };

inline double SQR(double x) { return x*x; }


// Helper functions in nlodishelper.cpp that need to be accessed outside
// e.g. for unit tests
int integrand_ILdip_massive(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata);
double ILdip_massive_Icd(double Q2, double z1, double x01sq, double mf, double xi, double x); 
double ILdip_massive_Iab(double Q2, double z1, double r, double mf, double xi);
double ILdip_massive_Omega_L_Const(double Q2, double z1, double r, double mf);

#endif 
