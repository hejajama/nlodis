#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "nlodis.hpp"
#include "datatypes.hpp"
#include <gsl/gsl_errno.h>
#include "dipole/bkdipole/bkdipole.hpp"
#include "gitsha1.h"

using namespace std;

/// Model configuration parameters
struct ModelConfig {
    std::string name;
    std::string datafile;
    double c2_alpha;
    double charm_mass;
    double sigma0_2;
    RunningCouplingScheme rc_scheme;
    double carlisle;
};

int main(int argc, char* argv[]) {

    cout << "# NLODIS code, git commit " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << endl;
    // Suppress GSL error handler for underflow errors during integration
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

    std::vector<ModelConfig> configs{
        {"KCBK parent", "/Users/hejajama/code/nlodisfit_bayesian/data/pd/bk_map.dat", 663, 1.4, 20.7, RunningCouplingScheme::PARENT,0.549267},
        {"KCBK smallest", "/Users/hejajama/code/nlodisfit_bayesian/data/balsd/bk_map.dat", 1.7, 1.25, 8.75, RunningCouplingScheme::SMALLEST,0},
        {"NLOBK MV smallest", "/Users/hejajama/Downloads/mv_bk.dat", 23, 1.04, 23.5, RunningCouplingScheme::SMALLEST,0},
        {"NLOBK MVgamma smallest", "/Users/hejajama/Downloads/mvgam_bk.dat", 1314.9257306, 1.2049379, 22.9017918, RunningCouplingScheme::SMALLEST,0},
        {"NLOBK MVgamma parent", "/Users/hejajama/Downloads/pd_nlo_bk.dat", std::pow(10, 3.88), 1.20, 24.3, RunningCouplingScheme::PARENT,0.554365}
    };
    
    double Q2 = 4.5;
    double xbj = 3.2e-3;
    double y =  1.3896E-02;
    double exp = 5.7189E-01;
    double experr =  0.009723674103;


    //cout <<"Computing at Q^2=" << Q2 << " GeV^2 and xbj=" << xbj << endl;

    NLODIS dis;
    dis.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/Downloads/mvgam_bk.dat"));
    // Running coupling scale
    dis.SetRunningCouplingC2(1314.9257306); 
    // The distance scale is set by the smallest dipole size
    dis.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    // Perform NLO calculation
    dis.SetOrder(Order::NLO);
    // Set charm quark mass 
    dis.SetQuarkMass(Quark::Type::C, 1.2049379);
// Set the proton transverse area
    dis.SetProtonTransverseArea(22.9017918, Unit::MB);
    dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH);
     dis.PrintConfiguration("# ");
     double F2 = dis.F2(Q2, xbj);
     cout << "F2 = " << F2 << endl;

   /* for (double mq = 0.2; mq >= 0.001; mq /= 1.5)
    {
        dis.SetQuarkMass(Quark::Type::LIGHT, mq);
        double FT = dis.FT(Q2, xbj);
        double FL = dis.FL(Q2, xbj);
        double F2 = FT + FL;
        //double FT = dis.FT(Q2, xbj);
        double sigmar = F2 - y*y/(1+SQR(1-y))*FL;
        //cout << "mq = " << mq << " GeV: F2 = " << F2 << ", FL = " << FL << ", FT = " << FT << ", sigma_r = " << sigmar << endl;
        cout << mq << " " << sigmar << endl;
    }
    return 0;
     
*/
   
    /*
dis.SetQuarkMass(Quark::Type::U, 0.01);
    dis.SetQuarkMass(Quark::Type::D, 0.01);
    dis.SetQuarkMass(Quark::Type::S, 0.01);
  
  //  Quark light(Quark::Type::LIGHT, 0.01);
    //dis.SetQuarks({light, Quark(Quark::Type::C, 1.2049379)});


    dis.PrintConfiguration();

    cout << "== EXP sigma_r(Q^2=" << Q2 << ", x=" << xbj << ") = " << exp << " +/- " << experr << endl;
    double fy = y*y/(1+SQR(1-y));
    double FT = dis.FT(Q2, xbj);
    double FL = dis.FL(Q2, xbj);
    double F2 = FT+FL;
    double sigmar = F2 - fy*FL;
    cout <<  " F2 = " << F2 << ", sigma_r = " << sigmar << endl;

    double photonproton_to_F2 = Q2 / (4.0 * M_PI * M_PI * Constants::AlphaEM);
    double sigmadip = photonproton_to_F2*(dis.Sigma_dip_d2b(Q2, xbj, Polarization::T) + dis.Sigma_dip_d2b(Q2, xbj, Polarization::L));
    double sigma_qg = photonproton_to_F2*(dis.Sigma_qg_d2b(Q2, xbj, Polarization::T) + dis.Sigma_qg_d2b(Q2, xbj, Polarization::L));
    double sigma_IC = photonproton_to_F2*(dis.Photon_proton_cross_section_LO_d2b(Q2, dis.GetDipole().X0(), Polarization::T) + dis.Photon_proton_cross_section_LO_d2b(Q2, dis.GetDipole().X0(), Polarization::L));
    cout << "Sigma_dip = " << sigmadip*22.9017*2.568 << endl;
    cout << "Sigma_qg = " << sigma_qg*22.9017*2.568 << endl;
    cout << "Sigma_IC = " << sigma_IC*22.9017*2.568 << endl;

    cout << sigma_IC << endl;

*/
        /*
                                                
    for (const auto& cfg : configs) {
        //NLODIS dis;
        dis.SetDipole(std::make_unique<BKDipole>(cfg.datafile));
        dis.SetRunningCouplingC2(cfg.c2_alpha);
        dis.SetRunningCouplingScheme(cfg.rc_scheme);
        dis.SetOrder(Order::NLO);
        dis.SetQuarkMass(Quark::Type::C, cfg.charm_mass);
        dis.SetProtonTransverseArea(cfg.sigma0_2, Unit::MB);

        double FT = dis.FT(Q2, xbj);
        double FL = dis.FL(Q2, xbj);
        double F2 = FT+FL;
        double sigmar = F2 - fy*FL;
        cout << cfg.name << " F2 = " << F2 << ", sigma_r = " << sigmar << "," << " -- Carlisle sigma_r = " << cfg.carlisle << endl;

        //dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::FREEZE);
        //dis.PrintConfiguration(); 
        //double F2_result_sharp = dis.F2(Q2, xbj);
        //dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH); 
        //double F2_result_smooth = dis.F2(Q2, xbj);
        //cout << cfg.name << " F2 = " << F2_result_sharp << " (FREEZE), " << F2_result_smooth << " (SMOOTH)" <<  " -- Carlisle " << cfg.carlisle << endl;
    }

    /*
    cout <<" == KCBK, parent ==" << endl;
    NLODIS kcbk_parent;
    kcbk_parent.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/code/nlodisfit_bayesian/data/pd/bk_map.dat"));
    kcbk_parent.SetRunningCouplingC2Alpha(663);
    kcbk_parent.SetRunningCouplingScheme(RunningCouplingScheme::PARENT);
    kcbk_parent.SetOrder(Order::NLO);
    kcbk_parent.SetQuarkMass(Quark::Type::C, 1.4);
    kcbk_parent.SetProtonTransverseArea(20.7, Unit::MB); // Set \sigma_0/2 = 20.7 mb
    double Q2=4.5;
    double xbj=3e-3;
    //double FL_kcbk_parent = kcbk_parent.FL(Q2, xbj);
    //double FT_kcbk_parent = kcbk_parent.FT(Q2, xbj);
    double F2_kcbk_parent = kcbk_parent.F2(Q2, xbj);
    //std::cout << "KCBK parent dipole FL(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FL_kcbk_parent*2*simga0_2_kcbk << std::endl;
    //std::cout << "KCBK parent dipole FT(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FT_kcbk_parent*2*simga0_2_kcbk << std::endl;
    std::cout << "KCBK parent dipole F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_kcbk_parent << std::endl;
    cout << endl;

    cout << " == KCBK, smallest ==" << endl;
    NLODIS kcbk_smallest;
    kcbk_smallest.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/code/nlodisfit_bayesian/data/balsd/bk_map.dat"));
    kcbk_smallest.SetRunningCouplingC2Alpha(1.7);
    kcbk_smallest.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    kcbk_smallest.SetOrder(Order::NLO);
    kcbk_smallest.SetQuarkMass(Quark::Type::C, 1.25);
    kcbk_smallest.SetProtonTransverseArea(8.75, Unit::MB);
    double F2_kcbk_smallest = kcbk_smallest.F2(Q2, xbj);
    std::cout << "KCBK smallest dipole F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_kcbk_smallest << std::endl;
    cout << endl;

    cout << "== NLOBK, smallest, MV ==" << endl;

    NLODIS mv;
    mv.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/Downloads/mv_bk.dat"));
    mv.SetRunningCouplingC2Alpha(23);
    mv.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    mv.SetOrder(Order::NLO);
    mv.SetQuarkMass(Quark::Type::C, 1.04);
    mv.SetProtonTransverseArea(23.5, Unit::MB); // Set \sigma_0/2 = 23.5 mb
    //double FL_mv = mv.FL(Q2, xbj);
    double F2_mv = mv.F2(Q2, xbj);
    std::cout << "NLOBK MV smallest F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_mv << std::endl;
    //std::cout << "MV smallest FL(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FL_mv*2*sigma0_2_mv << std::endl;
    cout << endl;

    cout << "== NLOBK, smallest, mvgamma ==" << endl;
    NLODIS mvgamma;
    mvgamma.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/Downloads/mvgam_bk.dat"));
    mvgamma.SetRunningCouplingC2Alpha(std::pow(10,2.9));
    mvgamma.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    mvgamma.SetQuarkMass(Quark::Type::C, 1.16 );
    mvgamma.SetOrder(Order::NLO);
    mvgamma.SetProtonTransverseArea(22.5, Unit::MB); // Set \sigma_0/2 = 22.5 mb
    double F2_mvgamma = mvgamma.F2(Q2, xbj);
    std::cout << "NLOBK MV gamma smallest F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_mvgamma << std::endl;
    cout << endl;

    cout << "== NLOBK, parent, mvgamma ==" << endl;
    NLODIS mvgamma_parent;
    mvgamma_parent.SetDipole(std::make_unique<BKDipole>("/Users/hejajama/Downloads/pd_nlo_bk.dat"));
    mvgamma_parent.SetRunningCouplingC2Alpha(std::pow(10.0, 3.88));
    mvgamma_parent.SetRunningCouplingScheme(RunningCouplingScheme::PARENT);
    mvgamma_parent.SetQuarkMass(Quark::Type::C, 1.2);
    mvgamma_parent.SetOrder(Order::NLO);
    mvgamma_parent.SetProtonTransverseArea(24.3, Unit::MB); //
    double F2_mvgamma_parent = mvgamma_parent.F2(Q2, xbj);
    std::cout << "NLOBK MV gamma parent F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_mvgamma_parent << std::endl;
    cout << endl;
    */

    /*
    // H1 data points for FL
    std::vector<std::pair<double, double>> q2_x_pairs = {
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

    cout << "Q^2 xbj FL_KCBK_parentdipole FL_MV FL_MVgamma" << endl;
    for (const auto& [q2, x] : q2_x_pairs) {
        // Use q2 and x here
        double FL_kcbk = kcbk_parentdipole.FL(q2, x)*2*simga0_2_kcbk;
        double FL_mv_val = mv.FL(q2, x)*2*sigma0_2_mv;
        double FL_mvgamma_val = mvgamma.FL(q2, x)*2*sigma0_2_mvgamma;

        cout << q2 << " " << x << " " << FL_kcbk << " " << FL_mv_val << " " << FL_mvgamma_val << endl;
    }

    /*NLODIS dis(argv[1]);
    dis.SetRunningCouplingC2Alpha(std::stod(argv[2]));
    double sigma0_2 = std::stod(argv[3])*2.568; // Convert mb to GeV^-2
    dis.SetOrder(Order::NLO);
    double Q2=10;
    double xbj=2e-4;

    

    cout << "#sigma(gamma+p;Q^2="<< Q2 << ",xbj="<< xbj << ",pol=L)=" << endl;
    double res = dis.FL(Q2, xbj);
    cout << res*sigma0_2*2 << endl;
*/
    return 0;
}