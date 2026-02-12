#include "dipoleamplitude.hpp"
#include "vector.hpp"
#include <string>
#include <memory>
#include <stdexcept>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>


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

std::string Dipole::GetString() const
{
    return "No info string implemented for the dipole amplitude";
}


/*
 * Saturation scale
 * Defined as N(r^2=2/Q_s^2) = Ns
 * 
 * Returns Q_s^2 in GeV^2
 */

struct SatscaleSolverHelper
{
    double Y;
    double Ns;
    const Dipole* N;
};

double SatscaleSolverHelperf(double r, void* p)
{
    auto* par = static_cast<SatscaleSolverHelper*>(p);
    return par->N->DipoleAmplitude(r, par->Y) - par->Ns;
}


double Dipole::SaturationScale(double Y, double Ns) const
{
    SatscaleSolverHelper helper{Y, Ns, this};
    constexpr int MAX_ITER = 1000;
    constexpr double ROOTFINDACCURACY = 0.00001;
    
    gsl_function f{};
    f.function = &SatscaleSolverHelperf;
    f.params = &helper;

    auto solver_deleter = [](gsl_root_fsolver* s) { gsl_root_fsolver_free(s); };
    auto solver = std::unique_ptr<gsl_root_fsolver, decltype(solver_deleter)>(
        gsl_root_fsolver_alloc(gsl_root_fsolver_bisection),
        solver_deleter
    );
    
    if (!solver) {
        throw std::runtime_error("Failed to allocate root solver");
    }
    
    gsl_root_fsolver_set(solver.get(), &f, MinR() * 1.0001, MaxR() * 0.999);
    
    int iter = 0;
    int status = GSL_CONTINUE;
    double min = 0, max = 0;
    
    while (status == GSL_CONTINUE && iter < MAX_ITER)
    {
        iter++;
        gsl_root_fsolver_iterate(solver.get());
        min = gsl_root_fsolver_x_lower(solver.get());
        max = gsl_root_fsolver_x_upper(solver.get());
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    
    }

    if (iter >= MAX_ITER) {
        throw std::runtime_error("Root finding failed for saturation scale at Y=" + std::to_string(Y));
    }

    return 2.0 / (gsl_root_fsolver_root(solver.get()) * gsl_root_fsolver_root(solver.get()));
}
