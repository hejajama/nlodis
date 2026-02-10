/*
 * General purpose interpolation class
 * Wrapper for GSL interpolator
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#pragma once

#include <memory>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;


/**
 * Interpolation method
 */
enum class InterpolationMethod {
    SPLINE
};

/// Custom deleters for GSL resources
struct GslSplineDeleter {
    void operator()(gsl_spline* s) const noexcept {
        if (s) gsl_spline_free(s);
    }
};

struct GslInterpAccelDeleter {
    void operator()(gsl_interp_accel* a) const noexcept {
        if (a) gsl_interp_accel_free(a);
    }
};


class Interpolator
{
    public:
        /**
         * @brief Create interpolator from two arrays.
         *
         * Pointers to the arrays are saved and used when Interpolator::Initialize()
         * is called. Also these pointers can be asked from the class later (see GetXData()
         * and GetYData()). The arrays are not used for evaluating the interpolator.
         * @param x array of x coordinates
         * @param y array of y coordinates
         * @param p number of points in arrays
         */
        Interpolator(double* x, double* y, std::size_t p);

        /**
         * Create interpolator from two std::vectors.
         *
         * The given vectors are not saved or referenced later.
         */
        Interpolator(const std::vector<double> &x, const std::vector<double> &y);
        Interpolator(const Interpolator& inter);
        Interpolator() { };
        ~Interpolator();
        void Clear();
        Interpolator& operator=(const Interpolator& inter);
        /**
         * Evaluate interpolator f(x)
         */
        double Evaluate(double x) const;
        /**
         * Evaluate 1st derivative of the interpolated function f'(x)
         */
        double Derivative(double x) const;    
        /**
         * Evaluate 2nd derivative of the interpolated function f''(x)
         */
        double Derivative2(double x) const;   
        /**
         * Select interpolation method (spline or bspline)
         */
        void SetMethod(InterpolationMethod m);

        /**
         * Initialize interpolator.
         *
         * This is done automatically at the constructor. If
         * interpolation method is changed (see SetMethod()), this
         * must be evaluated.
         */
        int Initialize();

        /**
         * Smallest x value supported
         */
        constexpr double MinX() const noexcept;
        /**
         * Largest x-value supported
         */
        constexpr double MaxX() const noexcept;

        std::vector<double> GetXData() const noexcept;
        std::vector<double> GetYData() const noexcept;
        gsl_spline* GetGslSpline() const;
        std::size_t GetNumOfPoints() const noexcept;
        InterpolationMethod GetMethod() const noexcept;
        
        constexpr bool Freeze() const noexcept;

        /**
         * Set what to do when interpolated function is evaluated
         * outside the interpolation region.
         *
         * @param f true if pre-saved values are returned
         */
        void SetFreeze(bool f);
        /**
         * Set value that is returned if interpolator is evaluated at
         * x<minx
         */
        void SetUnderflow(double min);
        /**
         * Set value that is returned if interpolator is evaluated at
         * x>maxx
         */
        void SetOverflow(double max);
        double UnderFlow() const noexcept;
        double OverFlow() const noexcept;

        /**
         * Set interpolation range (maxx), default is upper limit of the
         * xdata
         */
        void SetMaxX(double x);
        /**
         * Set interpolation range (minx), default is lower  limit of the
         * xdata
         */
        void SetMinX(double x);

        /**
         * Choose wether an error will be printed to stderr if interpolator
         * is evaluated outside the interpolation range.
         */
        void SetOutOfRangeErrors(bool set);

    private:
        /// Check that x values are monotonically increasing, throw on error
        void ValidateMonotonicIncreasing() const;
        InterpolationMethod method;
        std::vector<double> xdata, ydata;
        double minx, maxx;
        std::size_t points;
        bool ready;
        
        bool freeze;		// true if return freeze_under/overflow if
        double freeze_underflow;	// asked to evaluate interpolation
		double freeze_overflow;	// outside the spesified range
        
        // spline with RAII cleanup via unique_ptr
        std::unique_ptr<gsl_interp_accel, GslInterpAccelDeleter> acc;
        std::unique_ptr<gsl_spline, GslSplineDeleter> spline;
        static const int k=4;
        static const int ncoeffs = 12;
        static const int nbreak = ncoeffs-k+2;

        bool out_of_range_errors;


};


