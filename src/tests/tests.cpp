#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include "../nlodis.hpp"
#include "../integration.hpp"

#include <amplitudelib.hpp>

const std::string gbw_datafile = "gbw.dat";

/*
 * LO cross section tests 
 */

 TEST(LEADING_LOG_STRUCTURE_FUNCTIONS)
 {

    NLODIS dis("gbw.dat");

    double Q2 = 10.0; // GeV^2
    double xbj = 1;  // Test case datafile is generated such that it starts from x0=1
    
    dis.SetOrder(LO);

    double f2_ic=dis.F2(Q2,xbj);
    ASSERT_ALMOST_EQUAL(f2_ic, 0.0202973, 0.001);

    // Smaller x, tests interpolation and computation of the evolution rapidity
    double f2 = dis.F2(Q2, 1e-4);
    ASSERT_ALMOST_EQUAL(f2, 0.120141, 0.001);
 }


 
/*
 * Test that the GBW dipole behaves as expected
 */

TEST(GBW_DIPOLE_SATSCALE)
{
    // file gbw.dat corresponds to the GBW dipole amplitude
    // N(r,Y) = 1 - exp(-r^2 Q_s^2(Y)/4)
    // with Q_s^2(Y) = 1.0*exp(-lambda*Y) [GeV^2] and lambda=1/3
    AmplitudeLib N(gbw_datafile);

    double Ns=1-std::exp(-0.5);

    // Test initial condition
    ASSERT_ALMOST_EQUAL(N.SaturationScale(0, Ns), 1.0, 1e-3);

    // Test evolution rapidity that is not part of the grid
    double Y=std::log(0.01/0.001);
    // Test first without iniitalizing interpolation at fixed Y
    ASSERT_ALMOST_EQUAL(N.SaturationScale(Y, Ns), std::exp(1./3.*Y), 1e-3);
    // Then faster versoin with pre-initialized interpolation
    N.InitializeInterpolation(Y);
    ASSERT_ALMOST_EQUAL(N.SaturationScale(Y, Ns), std::exp(1./3.*Y), 1e-3);

}

TEST_MAIN()
