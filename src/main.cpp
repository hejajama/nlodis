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
#include "dipole/bkdipole/bkdipole.hpp"
#include "gitsha1.h"
#include "integration.hpp"

struct ModelConfig {
    std::string name="nlodis";
    std::string datafile;
    double c2_alpha;
    double charm_mass;
    double light_mass = 0.01; 
    double sigma0_2;
    RunningCouplingScheme rc_scheme;
    int nf=-1;
    Order order;
    double maxr=9999;
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
    std::cerr << "    --light_mass <value> [GeV] \\\n";
    std::cerr << "    --proton_area <value> [mb] \\\n";
    std::cerr << "    --rc_scheme <PARENT|SMALLEST> \\\n";
    std::cerr << "    [--order <LO|NLO>]\n";
    std::cerr << "    [--runmode <F2FL_GRID|HERA_FL>]\n";
    std::cerr << "    [--epsrel <value>] (relative precision for MC integration, default 0.01)\n";
    std::cerr << "    [--mcintpoints <value>] (number of points for MC integration, default 1e6)\n";
    std::cerr << "    [--cubamethod <method>] (Cuba integration method, default 'vegas')\n";
    std::cerr << "    [--nf <value>] (number of active flavors for running coupling, default determined from quark list)\n";
}

enum RunMode {
    F2FL_GRID,
    HERA_FL
};

int main(int argc, char* argv[]) {

    std::cout << "# NLODIS code, git commit " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << std::endl;

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
        } else if (flag == "--light_mass") {
            cfg.light_mass = std::stod(value);
            seen["light_mass"] = true;
        }
        else if (flag == "--epsrel")
        {
            dis.SetMCIntegrationEpsRel(std::stod(value));
            seen["epsrel"] = true;
        }
        else if (flag == "--nf")
        {
            cfg.nf = std::stoi(value);
            seen["nf"] = true;
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
            // Note: only vegas is currently supported
            dis.SetMCIntegrationMethod(value);
        }
        else if (flag == "--maxr")
        {
            cfg.maxr=std::stod(value);
            dis.SetMaxR(std::stod(value));
        }

        else {
            std::cerr << "Unknown flag: " << flag << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }
    
    std::vector<std::string> required = {"datafile", "C2", "charm_mass", "proton_area", "rc_scheme"};
    for (const auto& req : required) {
        if (!seen[req]) {
            std::cerr << "Error: missing required argument --" << req << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (print_header) {
        if (runmode == RunMode::F2FL_GRID)
            std::cout << "fit,x,Q2,F2 light,FL light,F2 charm,FL charm,F2_LO" << std::endl;
        else if (runmode == RunMode::HERA_FL)
            std::cout << "fit,x,Q2,FL" << std::endl;
    }
 
    dis.SetDipole(std::make_unique<BKDipole>(cfg.datafile));
    dis.SetRunningCouplingC2(cfg.c2_alpha);
    dis.SetRunningCouplingScheme(cfg.rc_scheme);
    dis.SetOrder(cfg.order);
    dis.SetQuarkMass(Quark::Type::C, cfg.charm_mass);
    dis.SetQuarkMass(Quark::Type::LIGHT, cfg.light_mass);
    dis.SetProtonTransverseArea(cfg.sigma0_2, Unit::MB);
    dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH);
    dis.SetDipole(std::make_unique<BKDipole>(cfg.datafile));
    dis.SetActiveFlavors(cfg.nf);



    if (auto* bk = dynamic_cast<BKDipole*>(&dis.GetDipole())) {
        bk->SetInterpolationMethod(LINEAR_LINEAR);   // fastest and typically accurate enough
    } else {
        throw std::runtime_error("Dipole is not BKDipole");
    }
    dis.PrintConfiguration("# ");

    Quark light(Quark::Type::LIGHT, cfg.light_mass);
    Quark charm(Quark::Type::C, cfg.charm_mass);

    std::vector<double> Q2vals = { 4.5, 45, 100};
    if (runmode == RunMode::F2FL_GRID) {
        for (const auto& Q2 : Q2vals) {
            double minx=1e-5;
            for (double x = minx; x <= 1e-2; x *= 1.5) {
                double FT_light, FL_light,FT_charm,FL_charm;
                dis.SetQuarks({light});
                try {
                    FT_light = dis.FT(Q2, x);
                    FL_light = dis.FL(Q2, x);

                    dis.SetQuarks({charm});
                    FT_charm = dis.FT(Q2,x);
                    FL_charm = dis.FL(Q2,x);

                } catch (const std::exception& e) {
                    std::cerr << "Error computing FT/FL for x=" << x << ", Q2=" << Q2
                            << ": " << e.what() << std::endl;
                    continue;
                }
                

                double F2_light = FT_light + FL_light;
                double F2_charm = FT_charm + FL_charm;
                dis.SetOrder(Order::LO);
                dis.SetQuarks({light, charm});
                double loF2 = dis.F2(Q2, x);
                dis.SetOrder(Order::NLO);
                std::cout << cfg.name << "," << x << "," << Q2 << "," << F2_light << "," << FL_light << "," << F2_charm <<","<< FL_charm <<"," << loF2 << std::endl;
              
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