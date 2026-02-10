/*
 * BK evolved dipole ampltiude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2020-2025
 */

#pragma once

#include "interpolation.hpp"
#include "../dipoleamplitude.hpp"
#include <vector>
#include <string>



enum AMPLITUDE_INTERPOLATION_METHOD
{
    SPLINE_LINEAR,      // spline in r, linear in y
    LINEAR_LINEAR       // linear in r and y
};

/**
 * Amplitude class
 *
 * Stores the loaded dipole amplitude (loaded by DataFile class), and
 * evaluates the amplitude at given r,Y.
 */
class BKDipole : public Dipole
{
    public:
        /**
         * Loader.
         * 
         * Constructor, loads the solution for the BK equation from
         * the given file and stores the result in memory.
         */
        BKDipole(std::string datafile);
    
        /* Load amplitude from the give data tables
         *
         * data[i][j] is dipole amplitude at rapidity yvals[i] for dipole size rvals[i]
         */
        BKDipole(std::vector< std::vector< double > > data, std::vector<double> yvals_, std::vector<double> rvals_);
    
        ~BKDipole();
    
        /**
         * Dipole amplitude at given rapidity Y
         *
         * Evaluates the dipole amplitude at given dipole size r and at
         * given rapidity Y (initial condition is at Y0, dipole is frozen in 0<Y<Y0).
         *
         * @param r dipole size in GeV^-1
         * @param Y rapidity at which the dipole is evaluated
         */
        double DipoleAmplitude(double r, double Y) const;



        /**
         * Initialize interpolation at given rapidity Y
         */
        void InitializeInterpolation(double Y);

        /**
         * Test if interpolator is initialized at given Bjorken-x
         */
        bool InterpolatorInitialized(double Y) const;

        /**
         * Create interpolator at given evolution rapidity nd return it
         */
        Interpolator* ConstructInterpolator(double Y) const;

    
        /**
         * Saturation scale
         *
         * Solve saturation scale defined as N(r^2=2/Q_s^2) = N_s
         * 
         * @return Saturation scale in GeV^2
         */
        double SaturationScale(double Y, double Ns) const;

        
        /**
         * Number or rapidity values in the BK solution
         */
        int YPoints() const;

        /**
         * Return ith rapidity value
         *
         * Can be used to avoid interpolation uncertainties 
         */
        double YValue(int yind) const;
         

        /**
         * Number of dipole sizes in the BK solution
         */
        int RPoints() const;

        /**
         * Minimum dipole size/rapidity
         */
        double MinR() const;
        double MinY() const;

        /**
         * Maximum dipole size/rapidity
         */
        double MaxR() const;
        double MaxY() const;
        
        /**
         * Set interpolation method
        */
        void SetInterpolationMethod(AMPLITUDE_INTERPOLATION_METHOD m) { interpolation_method = m;}

        /**
         * Specify whether an error message to stderr is printed for too
         * small/alrge dipoles 
         *
         * @return previous setting
         */

        bool SetOutOfRangeErrors(bool er);



        /**
         * Informative string about the dipole amplitude
         */
        std::string GetString() const;

	double X0() const { return 1; }
    
        
    private:
        // [yind][rind]
        std::vector< std::vector<double> > n;
        std::vector<double> yvals;
        std::vector<double> lnrvals;
        std::vector<double> rvals;
        Interpolator *interpolator;

        double interpolator_Y;  //! rapidity at which the interpolator is initialized

        double minr;
        double rmultiplier;
    
        AMPLITUDE_INTERPOLATION_METHOD interpolation_method;

        /**
         * Max dipole size for interpolation
         *
         * When interpolating in coordinate space,
         * force N(r>maxr_interpolate)=1 (avoid osciallatory 
         * artifacts from interpolation code)
         */
        double maxr_interpolate;
    
               
        int rpoints;    //! Number of dipole sizes

        bool out_of_range_errors = false;  //! If true, don't print "out of range" errors
        
        std::string info_string;
 
        std::string datafilename;   //! Name of the file where the amplitude is read, for error messages

        
};




const int INTERPOLATION_POINTS = 6;




const std::string AMPLITUDELIB_VERSION = "Bayesian NLO DIS fit (Casuga et al) 2025";
