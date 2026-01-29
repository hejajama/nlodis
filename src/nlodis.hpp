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

    
        void SetOrder(Order o) { order = o; }

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
      
};


#endif 
