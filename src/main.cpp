#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <csignal>
#include <map>
#include "nlodis.hpp"
#include "datatypes.hpp"
#include <gsl/gsl_errno.h>
#include "dipole/bkdipole/bkdipole.hpp"
#include "gitsha1.h"
#include "integration.hpp"
#include "debug.h"

struct ModelConfig {
    std::string name;
    std::string datafile;
    double c2_alpha;
    double charm_mass;
    double sigma0_2;
    RunningCouplingScheme rc_scheme;
    Order order;
};

static RunningCouplingScheme ParseRunningCouplingScheme(const std::string& s) {
    if (s == "PARENT") return RunningCouplingScheme::PARENT;
    if (s == "SMALLEST") return RunningCouplingScheme::SMALLEST;
    throw std::runtime_error("Invalid rc_scheme: " + s);
}

static Order ParseOrder(const std::string& s) {
    if (s == "LO") return Order::LO;
    if (s == "NLO") return Order::NLO;
    throw std::runtime_error("Invalid order: " + s);
}

static void PrintUsage(const char* prog) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog << " [--no-header] \\\n";
    std::cerr << "    --name <name> \\\n";
    std::cerr << "    --datafile <path> \\\n";
    std::cerr << "    --C2 <c2_alpha> \\\n";
    std::cerr << "    --charm_mass <value> [GeV] \\\n";
    std::cerr << "    --proton_area <value> [mb] \\\n";
    std::cerr << "    --rc_scheme <PARENT|SMALLEST> \\\n";
    std::cerr << "    [--order <LO|NLO>]\n";
    std::cerr << "    [--runmode <F2FL_GRID|HERA_FL>]\n";
}

enum RunMode {
    F2FL_GRID,
    HERA_FL
};

int main(int argc, char* argv[]) {
#ifdef DEBUG
    //EnableFloatingPointExceptions();
#endif

    std::cout << "# NLODIS code, git commit " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << std::endl;
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    (void)old_handler;

    bool print_header = true;
    RunMode runmode = RunMode::HERA_FL;
    ModelConfig cfg;
    cfg.order = Order::NLO;

    std::map<std::string, bool> seen;

    int i = 1;
    if (i < argc && std::string(argv[i]) == "--no-header") {
        print_header = false;
        ++i;
    }

    NLODIS dis;

    while (i < argc) {
        std::string flag(argv[i++]);
        
        if (i >= argc) {
            std::cerr << "Error: " << flag << " requires an argument\n";
            PrintUsage(argv[0]);
            return 1;
        }
        
        std::string value(argv[i++]);
        
        if (flag == "--name") {
            cfg.name = value;
            seen["name"] = true;
        } else if (flag == "--datafile") {
            cfg.datafile = value;
            seen["datafile"] = true;
        } else if (flag == "--C2") {
            cfg.c2_alpha = std::stod(value);
            seen["C2"] = true;
        } else if (flag == "--charm_mass") {
            cfg.charm_mass = std::stod(value);
            seen["charm_mass"] = true;
        } else if (flag == "--proton_area") {
            cfg.sigma0_2 = std::stod(value);
            seen["proton_area"] = true;
        } else if (flag == "--rc_scheme") {
            cfg.rc_scheme = ParseRunningCouplingScheme(value);
            seen["rc_scheme"] = true;
        } else if (flag == "--order") {
            cfg.order = ParseOrder(value);
            seen["order"] = true;
        }
        else if (flag == "--no_header")
            print_header = false;
        else if (flag == "--runmode") {
            if (value == "F2FL_GRID") {
                runmode = RunMode::F2FL_GRID;
            } else if (value == "HERA_FL") {
                runmode = RunMode::HERA_FL;
            } else {
                std::cerr << "Invalid run mode: " << value << "\n";
                PrintUsage(argv[0]);
                return 1;
            }
        }
        else if (flag == "--mcintpoints")
        {
            dis.SetMCIntegrationPoints(static_cast<int>(std::stod(value)));
        }
        else if (flag == "--cubamethod")
        {
            // This is a hidden flag to set the Cuba integration method, for testing purposes
            dis.SetMCIntegrationMethod(value);
        }

        else {
            std::cerr << "Unknown flag: " << flag << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }
    
    std::vector<std::string> required = {"name", "datafile", "C2", "charm_mass", "proton_area", "rc_scheme"};
    for (const auto& req : required) {
        if (!seen[req]) {
            std::cerr << "Error: missing required argument --" << req << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (print_header) {
        if (runmode == RunMode::F2FL_GRID)
            std::cout << "fit,x,Q2,F2,FL" << std::endl;
        else if (runmode == RunMode::HERA_FL)
            std::cout << "fit,x,Q2,FL" << std::endl;
    }
 
    dis.SetDipole(std::make_unique<BKDipole>(cfg.datafile));
    dis.SetRunningCouplingC2(cfg.c2_alpha);
    dis.SetRunningCouplingScheme(cfg.rc_scheme);
    dis.SetOrder(cfg.order);
    dis.SetQuarkMass(Quark::Type::C, cfg.charm_mass);
    dis.SetQuarkMass(Quark::Type::LIGHT, 0.005);
    dis.SetProtonTransverseArea(cfg.sigma0_2, Unit::MB);
    dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH);
    dis.PrintConfiguration("# ");

    if (runmode == RunMode::F2FL_GRID) {
        for (double x = 0.01; x >= 1e-5; x /= 2.0) {
            std::vector<double> Q2vals = {1.5, 4.5, 10, 45, 100, 200};
            for (const auto& Q2 : Q2vals) {
                double FT, FL;
                try {
                    FT = dis.FT(Q2, x);
                    FL = dis.FL(Q2, x);
                } catch (const std::exception& e) {
                    std::cerr << "Error computing FT/FL for x=" << x << ", Q2=" << Q2
                            << ": " << e.what() << std::endl;
                    continue;
                }
                double F2 = FT + FL;
                std::cout << cfg.name << "," << x << "," << Q2 << "," << F2 << "," << FL << std::endl;
            }
        }
    } else if (runmode == RunMode::HERA_FL) {
        std::vector<std::tuple<double, double>> data_points = {
             {1.5, 0.000028},
            {2.0, 0.000043},
            {2.5, 0.000059},
            {3.5, 0.000088},
            {5.0, 0.000129},
            {6.5, 0.000169},
            {8.5, 0.000224},
            {12.0, 0.000319},
            {15.0, 0.000402},
            {20.0, 0.000540},
            {25.0, 0.000687},
            {35.0, 0.000958},
            {45.0, 0.001210},
            {60.0, 0.001570},
            {90.0, 0.002430},
            {120.0, 0.003030},
            {150.0, 0.004020},
            {200.0, 0.005410},
            {250.0, 0.007360},
            {346.0, 0.009860},
            //{636.0, 0.018400}  
    };

        for (const auto& [Q2, x] : data_points) {
            double FL;
            try {
                FL = dis.FL(Q2, x);
            } catch (const std::exception& e) {
                std::cerr << "Error computing FL for x=" << x << ", Q2=" << Q2
                        << ": " << e.what() << std::endl;
                continue;
            }
            std::cout << cfg.name << "," << x << "," << Q2 << "," << FL << std::endl;
        }
    }

    return 0;
}