#include "gbw.hpp"
#include <cmath>
#include <sstream>

/*
 * GBW dipole amplitude
 * References: Phys.Rev.D 95 (2017) 1, 014030, e-Print: 1611.10100 [hep-ph]
 */

GBWDipole::GBWDipole(double Qs0sqr, double lambda, double gamma, double X0)
    : Qs0sqr_(Qs0sqr), lambda_(lambda), gamma_(gamma), X0_(X0) {}



/**
 * @brief Dipole amplitude N(r,Y)
 * 
 * N(r,Y) = 1 - exp(- (r^2 Qs^2(Y))^gamma /4 )
 * where Qs^2(Y) = Qs0^2 * exp(lambda * Y)
 * 
 * The amplitude is frozen for Y < Y0 = ln(1/X0)
 */
double GBWDipole::DipoleAmplitude(double r, double Y) const
{
    // Initial rapidity Y0 = ln(1/X0)
    double Y0 = std::log(1.0 / X0_);
    
    // Freeze amplitude for Y < Y0
    double Yeff = Y < Y0 ? Y0 : Y;
    
    // Saturation scale squared as a function of Y
    double Qs2_Y = Qs0sqr_ * std::exp(lambda_ * Yeff);
    
    // Argument of exponential: (r^2 Qs^2(Y))^gamma / 4
    double arg = std::pow(r * r * Qs2_Y, gamma_) / 4.0;
    
    // Dipole amplitude: N(r,Y) = 1 - exp(-arg)
    if (std::abs(arg) < 1e-7) {
        // For small arg, use expansion to avoid numerical issues
        return arg; // N ~ arg for small arg
    }
    else if (arg > 50) {
        // For large arg, N ~ 1 to avoid numerical issues with exp(-arg)
        return 1.0;
    }
    return 1.0 - std::exp(-arg);
}

std::string GBWDipole::GetString() const
{
    std::stringstream ss;
    ss << "GBW dipole amplitude: N(r,Y) = 1 - exp(-(r^2 Q_s^2(Y))^gamma/4)\n"
       << "  Q_s^2(Y) = Q_s0^2 * exp(lambda * Y)\n"
       << "  Q_s0^2 = " << Qs0sqr_ << " GeV^2\n"
       << "  lambda = " << lambda_ << "\n"
       << "  gamma = " << gamma_ << "\n"
       << "  X0 = " << X0_ << "\n"
       << "  Initial rapidity Y0 = ln(1/X0) = " << std::log(1.0 / X0_);
    
    return ss.str();
}
