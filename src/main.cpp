#include <iostream>
#include <string>
#include <vector>
#include "nlodis.hpp"
#include "datafile.hpp"
#include "datatypes.hpp"
#include <gsl/gsl_errno.h>

using namespace std;

int main(int argc, char* argv[]) {
    
    // Suppress GSL error handler for underflow errors during integration
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

    NLODIS kcbk_parentdipole("/Users/hejajama/code/nlodisfit_bayesian/data/pd/bk_map.dat");
    kcbk_parentdipole.SetRunningCouplingC2Alpha(663);
    kcbk_parentdipole.SetRunningCouplingScheme(PARENT);
    kcbk_parentdipole.SetOrder(NLO);
    kcbk_parentdipole.SetQuarkMass(Quark::Type::C, 1.4);
    kcbk_parentdipole.SetProtonTransverseArea(20.7, MB); // Set \sigma_0/2 = 20.7 mb
    double Q2=4.5;
    double xbj=3e-3;
    //double FL_kcbk_parent = kcbk_parentdipole.FL(Q2, xbj);
    //double FT_kcbk_parent = kcbk_parentdipole.FT(Q2, xbj);
    double F2_kcbk_parent = kcbk_parentdipole.F2(Q2, xbj);
    //std::cout << "KCBK parent dipole FL(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FL_kcbk_parent*2*simga0_2_kcbk << std::endl;
    //std::cout << "KCBK parent dipole FT(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FT_kcbk_parent*2*simga0_2_kcbk << std::endl;
    std::cout << "KCBK parent dipole F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_kcbk_parent << std::endl;
      
    NLODIS mv("/Users/hejajama/Downloads/mv_bk.dat");
    mv.SetRunningCouplingC2Alpha(23);
    mv.SetRunningCouplingScheme(SMALLEST);
    mv.SetOrder(NLO);
    mv.SetQuarkMass(Quark::Type::C, 1.04);
    mv.SetProtonTransverseArea(23.5, MB); // Set \sigma_0/2 = 23.5 mb
    //double FL_mv = mv.FL(Q2, xbj);
    double F2_mv = mv.F2(Q2, xbj);
    std::cout << "NLOBK MV smallest F2(Q2="<< Q2 << ", xbj="<< xbj << ") = " << F2_mv << std::endl;
    //std::cout << "MV smallest FL(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FL_mv*2*sigma0_2_mv << std::endl;

    /*
    NLODIS mvgamma("/Users/hejajama/Downloads/mvgam_bk.dat");
    mvgamma.SetRunningCouplingC2Alpha(std::pow(10,2.9));
    mvgamma.SetRunningCouplingScheme(SMALLEST);
    mvgamma.SetOrder(NLO);
    double sigma0_2_mvgamma = 22.5*2.5684624;
    ///double FL_mvgamma = mvgamma.FL(Q2, xbj);
    //std::cout << "MVgamma smallest FL(Q2="<< Q2 << ", xbj="<< xbj << ") = " << FL_mvgamma*2*sigma0_2_mvgamma << std::endl;

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
    double sigma0_2 = std::stod(argv[3])*2.5684624; // Convert mb to GeV^-2
    dis.SetOrder(NLO);
    double Q2=10;
    double xbj=2e-4;

    

    cout << "#sigma(gamma+p;Q^2="<< Q2 << ",xbj="<< xbj << ",pol=L)=" << endl;
    double res = dis.FL(Q2, xbj);
    cout << res*sigma0_2*2 << endl;
*/
    return 0;
}