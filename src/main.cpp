#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <csignal>
#include "nlodis.hpp"
#include "datatypes.hpp"
#include <gsl/gsl_errno.h>
#include "dipole/bkdipole/bkdipole.hpp"
#include "gitsha1.h"
#include "debug.h"


struct ModelConfig {
    std::string name;
    std::string datafile;
    double c2_alpha;
    double charm_mass;
    double sigma0_2;
    RunningCouplingScheme rc_scheme;
};

static RunningCouplingScheme ParseRunningCouplingScheme(const std::string& s) {
    if (s == "PARENT") return RunningCouplingScheme::PARENT;
    if (s == "SMALLEST") return RunningCouplingScheme::SMALLEST;
    throw std::runtime_error("Invalid rc_scheme: " + s);
}

static void PrintUsage(const char* prog) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog << " [--no-header] "
              << "name datafile c2_alpha charm_mass sigma0_2 rc_scheme carlisle "
              << "[name datafile c2_alpha charm_mass sigma0_2 rc_scheme carlisle ...]\n";
}

int main(int argc, char* argv[]) {
#ifdef DEBUG
    EnableFloatingPointExceptions();
#endif

    std::cout << "# NLODIS code, git commit " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << std::endl;
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    (void)old_handler;

    bool print_header = true;
    std::vector<ModelConfig> configs;

    int i = 1;
    if (i < argc && std::string(argv[i]) == "--no-header") {
        print_header = false;
        ++i;
    }

    const int fields = 6;
    int remaining = argc - i;
    if (remaining == 0 || remaining % fields != 0) {
        PrintUsage(argv[0]);
        return 1;
    }
    
    ModelConfig cfg;
    while (i < argc) {
        cfg.name = argv[i++];
        cfg.datafile = argv[i++];
        cfg.c2_alpha = std::stod(argv[i++]);
        cfg.charm_mass = std::stod(argv[i++]);
        cfg.sigma0_2 = std::stod(argv[i++]);
        cfg.rc_scheme = ParseRunningCouplingScheme(argv[i++]);
    }

    if (print_header) {
        std::cout << "fit,x,Q2,F2,FL" << std::endl;
    }

    std::vector<double> Q2vals = {1.5, 4.5, 10, 45, 100, 200};

    NLODIS dis;
    for (double x = 0.01; x >= 1e-5; x /= 2.0) {
        dis.SetDipole(std::make_unique<BKDipole>(cfg.datafile));
        dis.SetRunningCouplingC2(cfg.c2_alpha);
        dis.SetRunningCouplingScheme(cfg.rc_scheme);
        dis.SetOrder(Order::NLO);
        dis.SetQuarkMass(Quark::Type::C, cfg.charm_mass);
        dis.SetQuarkMass(Quark::Type::LIGHT, 0.005);
        dis.SetProtonTransverseArea(cfg.sigma0_2, Unit::MB);
        dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH);
        

        for (const auto& Q2 : Q2vals) {
            double FT,FL;
            try {
                FT = dis.FT(Q2, x);
                FL = dis.FL(Q2, x);
            } catch (const std::exception& e) {
                std::cerr << "Error computing FT and FL for x=" << x << ", Q2=" << Q2 << ": " << e.what() << std::endl;
                continue;
            }
            double F2 = FT + FL;
            std::cout << cfg.name << "," << x << "," << Q2 << "," << F2 << "," << FL << std::endl;
        }
    }

    return 0;
}