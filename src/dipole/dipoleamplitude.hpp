

/**
 * Virtual class for the dipole amplitude
 */

#pragma once

#include "vector.hpp"

class Dipole
{
    public:
    virtual ~Dipole() = default;
    /**
     * @brief Dipole-target scattering amplitude 
     * 
     * @param r Dipole size [GeV^-1]
     * @param Y Evolution rapidity
     * 
     */
        virtual double DipoleAmplitude(double r, double Y) const = 0;
        
        virtual double DipoleAmplitude(Vec r, Vec b, double Y) const;
        /** 
         * Support for b-dependent dipoles, the derived class does not implement these,
         * b-dependence is ignored
         */
        virtual double DipoleAmplitude(double r, double b, double Y) const;
        virtual double SaturationScale(double Y, double Ns) const = 0;

        virtual double X0() const = 0;

        /**
         * @brief Initialize dipole amplitude interpolation at given rapidity Y
         * 
         * If the implementation does not need this, this does not do anything
         * 
         * @param Y Evolution rapidity
         */
        virtual void InitializeInterpolation(double Y) {};

};