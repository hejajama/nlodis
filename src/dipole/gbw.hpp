

#pragma once

#include "dipoleamplitude.hpp"
#include <vector>
#include <string>




/*
 * @brief GBW-type dipole amplitude
 * 
 * This class implements a simple GBW-type dipole amplitude, for testing purposes
    * The amplitude is given by N(r,Y) = 1 - exp(- (r^2 Qs^2(Y))^gamma /4 )
    * where Qs^2(Y) = Qs0^2 * exp(lambda * Y)
    * The parameters are Qs0^2, lambda and gamma
    * The initial condition is at Y0 = ln(1/X0), and the amplitude is frozen for Y<Y0
 * This class is used for testing the NLO DIS code, and is not intended for phenomenological applications
 */

class GBWDipole : public Dipole
{
    public:
        GBWDipole(double Qs0sqr, double lambda, double gamma, double X0);
        GBWDipole() : GBWDipole(0.1, 0.3, 1.0, 1) {} // Default parameters

        /**
         * @brief Dipole amplitude at given rapidity dipole size r [GeV^-1] and evolution rapidity Y
         * 
         * @param r Dipole size [GeV^-1]
         * @param Y Evolution rapidity Y = ln 1/X0
         */
        double DipoleAmplitude(double r, double Y) const override;

        std::string GetString() const override;

    double X0() const override { return X0_; }
    
        
    private:
        double Qs0sqr_;
        double lambda_;
        double gamma_;
        double X0_;
};
