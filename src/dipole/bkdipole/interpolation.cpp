/*
 * AmplitudeLib tools
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2015
 */

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation.hpp"

// This is also defined in config.hpp, but if this class is used standalone
// it is more safe to not to include config.hpp but define this here.
#ifndef LINEINFO
    #define LINEINFO __FILE__ << ":" << __LINE__
#endif

typedef unsigned int uint;
using std::isinf;
using std::isnan;

/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator::Initialize()
{
    int status=0;
    out_of_range_errors = true;
    if (ready)
    {
        // unique_ptr will automatically clean up old resources
        acc.reset();
        spline.reset();
        ready=false;
    }
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            acc.reset(gsl_interp_accel_alloc());
            spline.reset(gsl_spline_alloc(gsl_interp_cspline, points));
            status = gsl_spline_init(spline.get(), xdata.data(), ydata.data(), points);
            break;
    }
    ready=true;
    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        return -1;
    }
    return 0;   //ok, there is no error handling at the moment...
}


double Interpolator::Evaluate(double x) const
{
    if (isnan(x) or isinf(x))
    {
        throw std::invalid_argument("Interpolator::Evaluate: invalid argument x=" + std::to_string(x));
    }
    
    if (!ready)
    {
        throw std::logic_error("Interpolator not initialized. Call Initialize() first.");
    }

    if (x<minx or x>maxx)
    {
		if (freeze)
		{
			if (x<minx) return freeze_underflow;
			else return freeze_overflow;
		}
		if (x < 0.9999*minx or x > 1.00001*maxx)	// if not true, no need to display error
        {
            //if (out_of_range_errors)
            //    cerr << "x=" << x << " is not within limits [" << minx << ", " << maxx << "], forcing "
            //        << "it in that interval! " << LINEINFO << endl;
        }
        throw std::out_of_range("x=" + std::to_string(x) + " is not within limits [" + std::to_string(minx) + ", " + std::to_string(maxx) + "]");

    }
    
    double res; int status;
    res=0;
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            status = gsl_spline_eval_e(spline.get(), x, acc.get(), &res);
            if (status)
            {
                throw std::runtime_error("Interpolation failed: " + std::string(gsl_strerror(status)) + ", x=" + std::to_string(x));
            }
            break;
        default:
            throw std::logic_error("Invalid interpolation method");
    }

    if (isnan(res) or isinf(res))
    {
        throw std::runtime_error("Interpolation produced invalid result at x=" + std::to_string(x));
    }

    
    return res;   
}

double Interpolator::Derivative(double x) const
{
    double res=0; int status=0;
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            status = gsl_spline_eval_deriv_e(spline.get(), x, acc.get(), &res);
            break;
        default:
            throw std::logic_error("Derivative not implemented for this interpolation method");
    }
    if (status)
        throw std::runtime_error("Derivative evaluation failed at x=" + std::to_string(x));

    return res;
}

double Interpolator::Derivative2(double x) const
{
    double res; int status=0;
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            status = gsl_spline_eval_deriv2_e(spline.get(), x, acc.get(), &res);
            break;
        default:
            throw std::logic_error("2nd derivative not implemented for this interpolation method");
    }

    if (status)
    {
        throw std::runtime_error("2nd derivative evaluation failed at x=" + std::to_string(x));
    }
    return res;

}

Interpolator::Interpolator(double *x, double *y, std::size_t p)
{
    points=p;
    xdata.assign(x, x + p);
    ydata.assign(y, y + p);
    minx = xdata.front();
    maxx = xdata.back();
    method = InterpolationMethod::SPLINE;

    ready=false;
    freeze=false;
    freeze_underflow = ydata.front();
    freeze_overflow = ydata.back();

    ValidateMonotonicIncreasing();

    Initialize();
}

void Interpolator::ValidateMonotonicIncreasing() const
{
    for (std::size_t i = 1; i < xdata.size(); ++i)
    {
        if (xdata[i-1] >= xdata[i])
        {
            throw std::invalid_argument("Grid points are not monotonically increasing at index " + std::to_string(i));
        }
    }
}

Interpolator& Interpolator::operator=(const Interpolator& inter)
{
    if (this == &inter)
        return *this;
    Clear();
    points = inter.points;
    xdata = inter.xdata;
    ydata = inter.ydata;
    minx = inter.minx;
    maxx = inter.maxx;
    method = inter.method;
    freeze = inter.freeze;
    freeze_underflow = inter.freeze_underflow;
    freeze_overflow = inter.freeze_overflow;
    out_of_range_errors = inter.out_of_range_errors;
    Initialize();
    return *this;
}

Interpolator::Interpolator(const std::vector<double> &x, const std::vector<double> &y)
{
    points = x.size();
    xdata = x;
    ydata = y;
    minx = xdata.front();
    maxx = xdata.back();

    method = InterpolationMethod::SPLINE;
    ready=false;
    freeze=false;
    freeze_underflow = ydata.front();
    freeze_overflow = ydata.back();

    ValidateMonotonicIncreasing();

    Initialize();
}

void Interpolator::SetMethod(InterpolationMethod m)
{
    method = m;
    
}

void Interpolator::Clear()
{
    acc.reset();
    spline.reset();
}

Interpolator::~Interpolator()
{
    Clear();

}

std::vector<double> Interpolator::GetXData() const noexcept
{
    return xdata;
}
std::vector<double> Interpolator::GetYData() const noexcept
{
    return ydata;
}
std::size_t Interpolator::GetNumOfPoints() const noexcept
{
    return points;
}
InterpolationMethod Interpolator::GetMethod() const noexcept
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator::Interpolator(const Interpolator& inter)
{
    points = inter.points;
    xdata = inter.xdata;
    ydata = inter.ydata;
    minx = inter.minx;
    maxx = inter.maxx;
    method = inter.method;
    freeze = inter.freeze;
    freeze_underflow = inter.freeze_underflow;
    freeze_overflow = inter.freeze_overflow;
    out_of_range_errors = inter.out_of_range_errors;
    
    ready=false;
    Initialize();
}

gsl_spline* Interpolator::GetGslSpline() const
{
    if (!spline)
    {
        throw std::logic_error("GSL spline is not initialized");
    }
    return spline.get();
}

constexpr double Interpolator::MinX() const noexcept
{
	return minx;
}

constexpr double Interpolator::MaxX() const noexcept
{
	return maxx;
}


constexpr bool Interpolator::Freeze() const noexcept
{
	return freeze;
}
void Interpolator::SetFreeze(bool f)
{
	freeze=f;
}
void Interpolator::SetUnderflow(double min)
{
	freeze_underflow=min;
}
 void Interpolator::SetOverflow(double max)
 {
	 freeze_overflow=max;
 }
double Interpolator::UnderFlow() const noexcept
{
	 return freeze_underflow;
}
double Interpolator::OverFlow() const noexcept
{
	return freeze_overflow;
}

void Interpolator::SetMaxX(double x)
{
    maxx=x;
}

void Interpolator::SetMinX(double x)
{
    minx=x;
}

void Interpolator::SetOutOfRangeErrors(bool er)
{
    out_of_range_errors=er;
}
