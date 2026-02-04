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

enum NcScheme
{
    FiniteNC,
    LargeNC
};

enum RunningCouplingScheme
{
    SMALLEST,
    PARENT
};


class NLODIS
{
    public:
        
        NLODIS(std::string bkdata);

        double F2(double Q2, double xbj);
        double FL(double Q2, double xbj);
        double Photon_proton_cross_section(double Q2, double xbj, Polarization pol);

         double Photon_proton_cross_section_LO(double Q2, double xbj, Polarization pol);

         /// NLO caluclation ingredients 
        double Sigma_dip(double Q2, double xbj, Polarization pol);
        double Sigma_qg(double Q2, double xbj, Polarization pol);

         ////


    
        void SetOrder(Order o) { order = o; }
        double GetMaxR() const { return maxr; }
        
        double Alphas(double r) const;
        AmplitudeLib& GetDipole() { return dipole; }

        double z2_lower_bound(double xbj, double Q2);

        double TripoleAmplitude(double x01, double x02, double x21, double Y); 

        double EvolutionRapidity(double xbj, double Q2, double z2);

        void SetNcScheme(NcScheme scheme_) { nc_scheme = scheme_; }
        void SetRunningCouplingScheme(RunningCouplingScheme rc_) { rc_scheme = rc_; }
        double RunningCouplinScale(double x01, double x02, double x21);

        void SetRunningCouplingC2Alpha(double c2) { C2_alpha = c2; }
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
        NcScheme nc_scheme = LargeNC;
        RunningCouplingScheme rc_scheme = SMALLEST;
        const double Q0sqr = 1; // Non-perturbative target scale, should match the one used in the NLO DIS fit!
        double C2_alpha = 10.0; // Scale factor in the coordinate space running coupling
      
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

// Sigma_qg longitudinal part helper
int integrand_ILqgunsub_massive(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata);
double ILNLOqg_massive_tripole_part_I2_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double G_integrand_simplified(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda, double y);
double ILNLOqg_massive_dipole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I3_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);
#endif 
