#include "dipoleamplitude.hpp"
#include "vector.hpp"

double Dipole::DipoleAmplitude(Vec r, Vec b, double Y) const
{
    // Ignore b dependence for now, just call the r-dependent version
    return DipoleAmplitude(r.Len(), Y);
}

double Dipole::DipoleAmplitude(double r, double b, double Y) const
{
    // Ignore b dependence for now, just call the r-dependent version
    return DipoleAmplitude(r, Y);
}