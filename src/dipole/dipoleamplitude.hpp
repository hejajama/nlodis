

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

        /**
         * Saturation scale
         *
         * Solve saturation scale defined as N(r^2=2/Q_s^2) = N_s
         * 
         * @return Saturation scale in GeV^2
         */
        virtual double SaturationScale(double Y, double Ns) const;

        /**
         * @brief X0 parameter: initial condition of fit
         * 
         * The initial condition is taken at evolution rapidity Y = ln 1/X0
         * 
         * TODO: it has not been tested that this code works for fits with X0 != 1
         * 
         */
        virtual double X0() const = 0;

        /**
         * @brief Initialize dipole amplitude interpolation at given rapidity Y
         * 
         * If the implementation does not need this, this does not do anything
         * 
         * @param Y Evolution rapidity
         */
        virtual void InitializeInterpolation(double Y) {};

        /**
         * @brief Informative string about the dipole amplitude
         */
        virtual std::string GetString() const;

        /**
         * @brief Minimum dipole size where dipole amplitude is available
         * 
         * One can assume that N(r<r_min)=0
         * 
         * @return Minimum dipole size [GeV^-1]
         */
        virtual double MinR() const { return 1e-10; }

        /**
         * @brief Maximum dipole size where dipole amplitude is available
         * 
         * One can assume that N(r>r_max)=1
         * 
         * @return Maximum dipole size [GeV^-1]
         */
        virtual double MaxR() const { return 1e3; }

};