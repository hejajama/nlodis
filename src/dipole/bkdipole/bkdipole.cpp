/*
 * BK evolved dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>,2020-2025
 */

#include "../dipoleamplitude.hpp"
#include "bkdipole.hpp"
#include "datafile.hpp"
#include <string>
#include <cmath>

#include <string>
#include <sstream>
#include <stdexcept>


#include <algorithm>

using std::isinf;
using std::isnan;

int FindIndex(double val, const std::vector<double>& vec);

inline double SQR(double x) { return x*x; }
#define LINEINFO __FILE__ << ":" << __LINE__

/*
 * Load data from a given file
 */
BKDipole::BKDipole(std::string datafile)
{
    // Read BK solution
    datafilename = datafile;
    DataFile data(datafile);
    auto [n_data, y_data] = data.GetData();
    n = n_data;
    yvals = y_data;

    // Initialize
    interpolation_method = SPLINE_LINEAR;
    
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();

    for (int i = 0; i < rpoints; i++)
    {
        double tmpr = minr * std::pow(rmultiplier, i);
        lnrvals.push_back(std::log(tmpr));
        rvals.push_back(tmpr);
    }

    interpolator_Y = -1;  // if >=0, interpolator is initialized

    std::stringstream ss;
    ss << "Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " initial rapidity " << MinY() << ", maximum rapidity "
        << yvals[yvals.size() - 1]
        << " Q_{s,0}^2(initial rapidity) = " << SaturationScale(0, 1 - std::exp(-0.5)) << " GeV^2 [ N(r^2=2/Q_s^2) = " << 1 - std::exp(-0.5) << "]";
    info_string = ss.str();
}

/*
 * Constructor, dipole data is given in an array
 *
 */
BKDipole::BKDipole(std::vector<std::vector<double>> data, std::vector<double> yvals_, std::vector<double> rvals_)
{
    interpolation_method = SPLINE_LINEAR;
    minr = rvals_[0];
    rmultiplier = rvals_[1] / rvals_[0];
    rpoints = rvals_.size();
    
    rvals = std::move(rvals_);
    n = std::move(data);
    yvals = std::move(yvals_);
    
    for (int i = 0; i < rpoints; i++)
    {
        double tmpr = std::log(minr * std::pow(rmultiplier, i));
        double tmpr2 = rvals[i];
        if (std::abs(std::exp(tmpr) / tmpr2 - 1.0) > 0.01)
        {
            throw std::runtime_error("Dipole sizes do not form a logarithmic grid! tmpr " + std::to_string(tmpr) + " r2 " + std::to_string(rvals[i]));
        }
        lnrvals.push_back(tmpr);
    }
    
    interpolator_Y = -1;
    
    std::stringstream ss;
    ss << "#BKDipole initialized , minr: " << minr
       << " maxr: " << MaxR() << " rpoints: " << rpoints << " Y0=" << MinY() << ", max Y=" << MaxY()
       << " Q_{s,0}^2(Y=0) = " << SaturationScale(0, 1. - std::exp(-0.5)) << " GeV^2 [ N(r^2=2/Q_s^2, Y=0) = " << 1 - std::exp(-0.5) << "]";
    info_string = ss.str();
    
}

/*
 * Release reserved memory
 */

BKDipole::~BKDipole()
{
    if (interpolator_Y>=0)
    {
        delete interpolator;
    }
    /*delete[] tmpnarray;
    delete[] tmprarray;*/
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 */
double BKDipole::DipoleAmplitude(double r, double y)  const
{
    if (isnan(r) or isinf(r))
    {
	throw std::invalid_argument("r=" + std::to_string(r) + ", y=" + std::to_string(y) + " at BKDipole::DipoleAmplitude()");    
    }
    
    if (r < MinR() or r > MaxR() )
    {
        if (out_of_range_errors)
            cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
                << " asked r=" << r << " " << LINEINFO
                << endl;
        if (r<MinR()) r=MinR()*1.000001; else if (r>MaxR()) r=MaxR()*0.999999;
    }
    
    // If Y0>0, the dipole is frozen in [0,Y0]
    if (y>=0 and y<MinY())
        y=MinY();
    
    
    if (y<MinY() or y>MaxY() )
    {
        //if (out_of_range_errors)
            throw std::runtime_error("y must be between limits [" + std::to_string(0) + ", " + std::to_string(yvals[yvals.size()-1]) + "], asked y=" + std::to_string(y) + " datafile " + datafilename);
    }

    
    // Use already initialized interpolator
    if (InterpolatorInitialized(y))
    {
        double result=0;

        // Can't interpolate (too large dipole), return 1.0
        if (r >= maxr_interpolate and maxr_interpolate>0)
        {
            return 1.0;
        }

        result = interpolator->Evaluate(r);
        
        if (result>1) return 1;              // limit N(r)<=1
        if (result<0) return 0;              // limit N >= 0
        return result;

    }
    
    int rind = FindIndex(r, rvals);
    int yind = FindIndex(y, yvals);
    
    // Check interpolation method, fastest is linear in r and y
    if (interpolation_method == LINEAR_LINEAR)
    {
        int rind2 = rind+1;
        int yind2 = yind+1;
        if (rind2 >= rvals.size()) return 1.0;
        if (yind2 >= yvals.size()) return n[yind][rind]; // this should never happen
        double r1_ = rvals[rind];
        double r2_ = rvals[rind2];
        double y1_ = yvals[yind];
        double y2_ = yvals[yind2];
        double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (n[yind][rind])*(r2_ - r)*(y2_ - y) + (n[yind][rind2])*(r - r1_)*(y2_ - y) + (n[yind2][rind])*(r2_ - r)*(y - y1_) + (n[yind2][rind2])*(r - r1_)*(y - y1_) );
        
        // Handle possible floating point errors
        if (bilin>1) return 1;
        if (bilin<0) return 0;
        return bilin;
        
    }

    // Initialize new interpolator and use it
    // Note: This is relatively slow. User should
    // usually InitializeInterpolation and then
    // multiple evaluations are fast at the same Y
    int interpolation_points = INTERPOLATION_POINTS;

    int interpolation_start, interpolation_end;
    if (rind - interpolation_points/2 < 0)
    {
        interpolation_start=0;
        interpolation_end=interpolation_points;
    }
    else if (rind + interpolation_points/2 > RPoints()-1 )
    {
        interpolation_end = RPoints()-1;
        interpolation_start = RPoints()-interpolation_points-3;
    }
    else
    {
        interpolation_start = rind - interpolation_points/2;
        interpolation_end = rind + interpolation_points/2;
    }

    int interpo_points = interpolation_end - interpolation_start + 1;
    std::vector<double> tmparray(interpo_points);
    std::vector<double> tmpxarray(interpo_points);
    for (int i = interpolation_start; i <= interpolation_end; i++)
    {
        tmpxarray[i - interpolation_start] = rvals[i];
        tmparray[i - interpolation_start] = n[yind][i];

        // Interpolate in y if possible
        if (yind < static_cast<int>(yvals.size() - 1))
        {
            tmparray[i - interpolation_start] = n[yind][i]
                + (y - yvals[yind]) * (n[yind + 1][i] - n[yind][i])
                / (yvals[yind + 1] - yvals[yind]);
        }
    }

    Interpolator interp(tmpxarray, tmparray);
    interp.Initialize();
    double result = interp.Evaluate(r);

    if (result < 0) return 0;
    if (result > 1) return 1;

    return result;

}

bool BKDipole::InterpolatorInitialized(double Y) const
{
    // Check Y=0 separately
    if (std::abs(Y) < 0.001 and std::abs(interpolator_Y) < 0.001)
        return true;
    
    if (interpolator_Y >= 0 and std::abs(Y - interpolator_Y)/std::min(Y, interpolator_Y) < 0.001)
        return true;
    else
        return false;
}




/*
 * Initializes dipole interpolation at given rapidity
 */
void BKDipole::InitializeInterpolation(double Y)
{
    // Between 0<Y<MinY() use MinY
    if (Y>=0 and Y<=MinY()) Y=MinY();
    
	if (Y < 0)
	{
		throw std::runtime_error("Asked to initialize interpolator at negative Y=" + std::to_string(Y));
	}

	if (Y>MaxY())
	{
		throw std::runtime_error("Asked to initialize interpolator with too large y=" + std::to_string(Y) + ", maxy=" + std::to_string(MaxY()));
	}
	
    if (InterpolatorInitialized(Y)) return;    // Already done it
    
    // Remove old interpolator if constructed
    if (interpolator_Y>=0)
    {
        interpolator_Y=-1;
        delete interpolator;
    }
    
    std::vector<double> tmpnarray(rpoints);
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = rvals[i];
        if (i==0) tmpr*=1.0001; if (i==rpoints-1) tmpr*=0.9999;
        tmpnarray[i] = DipoleAmplitude(tmpr, Y);
    }
    interpolator = new Interpolator(rvals, tmpnarray);

    interpolator->Initialize();
    interpolator_Y = Y;
    
    maxr_interpolate=-1;
    int iter=0;
    
    // Find rmax s.t. N(r) \approx 1 at r>rmax
    double step=2; double prevr=0.1;
    const int MAXITER=40;
    for (double r=0.01; r<=MaxR(); r+=step)
    {
        if (DipoleAmplitude(r,Y)>=0.99999)
        {
            if (step<1e-2)
            {
                maxr_interpolate = r;
                break;
            }
            step /= 1.5;
            r=prevr;
        }
        prevr=r;
        iter++;
        if (iter > MAXITER)
        {
            maxr_interpolate=-1;
            break;  // Didn't find, dont force any upper limit
        }
    }
}



/*
 * Initializes interpolator and returns it
 */
Interpolator* BKDipole::ConstructInterpolator(double Y) const
{
    std::vector<double> tmpnvals;
    for (int i=0; i<rpoints; i++)
    {
        tmpnvals.push_back(DipoleAmplitude(rvals[i], Y));
    }
    Interpolator* inter = new Interpolator(rvals, tmpnvals);
    inter->Initialize();

    return inter;
        
}


int BKDipole::RPoints() const
{
    return rpoints;
}
double BKDipole::MinR() const
{
    return minr;
}

double BKDipole::MaxR() const
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int BKDipole::YPoints() const
{
    return yvals.size();
}

double BKDipole::MaxY() const
{
    return yvals[yvals.size()-1];
}

double BKDipole::MinY() const
{
    return yvals[0];
}


bool BKDipole::SetOutOfRangeErrors(bool er)
{
    bool outofrange = out_of_range_errors;
    out_of_range_errors=er;
    return outofrange;    
}


std::string BKDipole::GetString() const
{
	return info_string;
}



/*
 * Return ith rapidity value
 */
double BKDipole::YValue(int yind) const
{
    if (yind >= YPoints() or yind<0)
    {
        throw std::runtime_error("Asked rapidity of index " + std::to_string(yind) +", but the maximum index is " + std::to_string(YPoints()-1) );
        exit(1);
    }
    return yvals[yind];
}


int FindIndex(double val, const std::vector<double>& vec)
{
    if (vec.empty() || val < vec.front())
        return -1;

    auto it = std::upper_bound(vec.begin(), vec.end(), val);
    return static_cast<int>(it - vec.begin()) - 1;
}