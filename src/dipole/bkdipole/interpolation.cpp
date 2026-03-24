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
            if (interpolate_log) {
                std::vector<double> log_xdata(points), log_ydata(points);
                for (std::size_t i = 0; i < points; ++i) {
                    log_xdata[i] = std::log(xdata[i]);
                    log_ydata[i] = std::log(ydata[i]);
                }
                status = gsl_spline_init(spline.get(), log_xdata.data(), log_ydata.data(), points);
            } else {
                status = gsl_spline_init(spline.get(), xdata.data(), ydata.data(), points);
            }
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

    if (interpolate_log && x <= 0.0)
    {
        throw std::invalid_argument("Interpolator::Evaluate: x must be strictly positive for log interpolation, got x=" + std::to_string(x));
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
            if (interpolate_log) {
                double log_x = std::log(x);
                double log_res = 0;
                status = gsl_spline_eval_e(spline.get(), log_x, acc.get(), &log_res);
                if (status)
                {
                    throw std::runtime_error("Interpolation failed: " + std::string(gsl_strerror(status)) + ", x=" + std::to_string(x));
                }
                res = std::exp(log_res);
            } else {
                status = gsl_spline_eval_e(spline.get(), x, acc.get(), &res);
                if (status)
                {
                    throw std::runtime_error("Interpolation failed: " + std::string(gsl_strerror(status)) + ", x=" + std::to_string(x));
                }
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
    if (interpolate_log && x <= 0.0)
    {
        throw std::invalid_argument("Interpolator::Derivative: x must be strictly positive for log interpolation, got x=" + std::to_string(x));
    }

    double res=0; int status=0;
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            if (interpolate_log) {
                // f(x) = exp(g(log(x))), so f'(x) = f(x) * g'(log(x)) / x
                double log_x = std::log(x);
                double log_f = 0, g_prime = 0;
                int s1 = gsl_spline_eval_e(spline.get(), log_x, acc.get(), &log_f);
                int s2 = gsl_spline_eval_deriv_e(spline.get(), log_x, acc.get(), &g_prime);
                if (s1 || s2)
                    throw std::runtime_error("Derivative evaluation failed at x=" + std::to_string(x));
                res = std::exp(log_f) * g_prime / x;
            } else {
                status = gsl_spline_eval_deriv_e(spline.get(), x, acc.get(), &res);
            }
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
    if (interpolate_log && x <= 0.0)
    {
        throw std::invalid_argument("Interpolator::Derivative2: x must be strictly positive for log interpolation, got x=" + std::to_string(x));
    }

    double res = 0; int status=0;
    switch(method)
    {
        case InterpolationMethod::SPLINE:
            if (interpolate_log) {
                // f(x) = exp(g(log(x)))
                // f''(x) = f(x)/x^2 * [ (g'(u))^2 + g''(u) - g'(u) ]  where u=log(x)
                double log_x = std::log(x);
                double log_f = 0, g_prime = 0, g_prime2 = 0;
                int s1 = gsl_spline_eval_e(spline.get(), log_x, acc.get(), &log_f);
                int s2 = gsl_spline_eval_deriv_e(spline.get(), log_x, acc.get(), &g_prime);
                int s3 = gsl_spline_eval_deriv2_e(spline.get(), log_x, acc.get(), &g_prime2);
                if (s1 || s2 || s3)
                    throw std::runtime_error("2nd derivative evaluation failed at x=" + std::to_string(x));
                res = std::exp(log_f) / (x * x) * (g_prime * g_prime + g_prime2 - g_prime);
            } else {
                status = gsl_spline_eval_deriv2_e(spline.get(), x, acc.get(), &res);
            }
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
    interpolate_log = false;
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
    for (std::size_t i = 0; i < xdata.size(); ++i)
    {
        if (isnan(xdata[i]) || isinf(xdata[i]))
        {
            throw std::invalid_argument("Grid x value is not finite at index " + std::to_string(i));
        }
    }

    for (std::size_t i = 1; i < xdata.size(); ++i)
    {
        if (xdata[i-1] >= xdata[i])
        {
            throw std::invalid_argument("Grid points are not monotonically increasing at index " + std::to_string(i));
        }
    }
}

void Interpolator::ValidateLogCompatible() const
{
    for (std::size_t i = 0; i < xdata.size(); ++i)
    {
        if (xdata[i] <= 0.0)
            throw std::invalid_argument("Log interpolation requires all x values to be strictly positive; x[" + std::to_string(i) + "]=" + std::to_string(xdata[i]));
    }
    for (std::size_t i = 0; i < ydata.size(); ++i)
    {
        if (ydata[i] <= 0.0)
            throw std::invalid_argument("Log interpolation requires all y values to be strictly positive; y[" + std::to_string(i) + "]=" + std::to_string(ydata[i]));
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
    interpolate_log = inter.interpolate_log;
    Initialize();
    return *this;
}

Interpolator::Interpolator(const std::vector<double> &x, const std::vector<double> &y, bool interpolate_log)
{
    this->interpolate_log = interpolate_log;
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
    if (interpolate_log)
        ValidateLogCompatible();

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
    interpolate_log = inter.interpolate_log;
    
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
