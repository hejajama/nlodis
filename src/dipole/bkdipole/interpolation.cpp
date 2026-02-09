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
        // Interpolator is already initialized: clear it (GSL part)
        // and re-initialize
        switch(method)
        {
            case INTERPOLATE_SPLINE:
                if (spline != NULL)
                {
                    gsl_spline_free(spline);
                    spline=NULL;
                }
                if (acc != NULL)
                {
                    gsl_interp_accel_free(acc);
                    acc=NULL;
                }
                break;
        }
        ready=false;
        
    }
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc(gsl_interp_cspline, points);
            status = gsl_spline_init(spline, xdata.data(), ydata.data(), points);
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


double Interpolator::Evaluate(double x)
{
    if (isnan(x) or isinf(x))
    {
        cerr << "Trying to evaluate interpolator with x=" << x << " at " << LINEINFO << endl;
        exit(1);
    }
    
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator::Initialize()?" << endl;
        return 0;
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
    
    double res, yerr; int status;
    res=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_e(spline, x, acc, &res);
            if (status)
            {
                cerr << "Interpolation failed at " << LINEINFO << ", error " << gsl_strerror(status)
                 << " (" << status << "), x=" << x << ", minx=" << xdata[0]
                 << ", maxx=" << xdata[points-1] << ", result=" << res << endl;
                 exit(1);
            }
            break;
        default:
            cerr << "Interpolation method is invalid! " << LINEINFO << endl;
            exit(1);
    }

    if (isnan(res) or isinf(res))
    {
        cerr << "Interpolation at x=" << x << " gives " << res << endl;
		return 0;
        exit(1);
    }

    
    return res;   
}

double Interpolator::Derivative(double x)
{
    double res=0; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv_e(spline, x, acc, &res);
            break;
        default:
            cerr << "Derivative is not implemented for this interpolation method!" << " " << LINEINFO << endl;
            exit(1);
    }
    if (status)
        cerr << "An error occurred while evaluating the derivative at x=" << x
        << " result " << res << " " << LINEINFO << endl;

    return res;
}

double Interpolator::Derivative2(double x)
{
    double res; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv2_e(spline, x, acc, &res);
            break;
    }

    if (status)
    {
        cerr << "2nd derivative interpolation failed at x=" << x <<", result "
        << res << " " << LINEINFO << endl;
    }
    return res;

}

Interpolator::Interpolator(double *x, double *y, int p)
{
    points=p;
    xdata.assign(x, x + p);
    ydata.assign(y, y + p);
    minx = xdata.front();
    maxx = xdata.back();
    method = INTERPOLATE_SPLINE;

    ready=false;
    freeze=false;
    freeze_underflow = ydata.front();
    freeze_overflow = ydata.back();

    for (int i=0; i<p; i++)
    {
        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
                exit(1);
            }
        }
    }

    Initialize();
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

    for (uint i=0; i<x.size(); i++)
    {
        xdata[i]=x[i];
        ydata[i]=y[i];

        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
                exit(1);
            }
        }
    }
    method = INTERPOLATE_SPLINE;
    ready=false;
    freeze=false;
    freeze_underflow = ydata.front();
    freeze_overflow = ydata.back();

    Initialize();
}

void Interpolator::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
    
}

void Interpolator::Clear()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            if (spline != NULL)
            {
                gsl_spline_free(spline);
                spline=NULL;
            }
            if (acc != NULL)
            {
                gsl_interp_accel_free(acc);
                acc=NULL;
            }
            break;
    }
}

Interpolator::~Interpolator()
{
    Clear();

}

std::vector<double> Interpolator::GetXData()
{
    return xdata;
}
std::vector<double> Interpolator::GetYData()
{
    return ydata;
}
int Interpolator::GetNumOfPoints() const
{
    return points;
}
INTERPOLATION_METHOD Interpolator::GetMethod() const
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
    if (spline == NULL)
    {
        cerr << "GSL spline is not initialized! " << LINEINFO << endl;
        exit(1);
    }
    return spline;
}

double Interpolator::MinX()
{
	return minx;
}

double Interpolator::MaxX()
{
	return maxx;
}


bool Interpolator::Freeze()
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
double Interpolator::UnderFlow()
{
	 return freeze_underflow;
}
double Interpolator::OverFlow()
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
