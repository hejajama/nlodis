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
    Quark u; u.type = Quark::U; u.mass = 0.14; u.charge = 2.0/3.0;
    Quark d; d.type = Quark::D; d.mass = 0.14; d.charge = -1.0/3.0;
    Quark s; s.type = Quark::S; s.mass = 0.14;  s.charge = -1.0/3.0;
    Quark c; c.type = Quark::C; c.mass = 1.4;   c.charge = 2.0/3.0;
    std::vector<Quark> quark_list = {u,d,s,c};
    dis.SetQuarks(quark_list);

    double Q2 = 10.0; // GeV^2
    double xbj = 1;  // Test case datafile is generated such that it starts from x0=1
    
    dis.SetOrder(LO);

    

    double f2_ic=dis.F2(Q2,xbj);
    ASSERT_ALMOST_EQUAL(f2_ic, 0.0202973+0.00543532, 0.0005); // light + charm

    // Smaller x, tests interpolation and computation of the evolution rapidity
    double f2 = dis.F2(Q2, 1e-4);
    ASSERT_ALMOST_EQUAL(f2, 0.120141+0.0538897, 0.001); // light +charm
 }

/*
 * Tripole i.e. 1-S_{012}
 * S_{012} defined in https://arxiv.org/pdf/2007.01645 eq 13
 */
TEST(TRIPOLE_AMPLITUDE)
{
    NLODIS dis("gbw.dat");

    dis.SetNcScheme(LargeNC);
    double x01 = 1.0;
    double x02 = 2.0;
    double x21 = std::sqrt( SQR(x01) + SQR(x02) - 2.0*x01*x02*0.5 ); // angle cos=0.5
    double Y = 2.0;

    double tripole_largeNC = dis.TripoleAmplitude(x01, x02, x21, Y);
    double N02 = dis.GetDipole().DipoleAmplitude(x02, Y);
    double N12 = dis.GetDipole().DipoleAmplitude(x21, Y);
    double expected_largeNC = 1.0 - (1.0 - N02)*(1.0 - N12);
    ASSERT_ALMOST_EQUAL(tripole_largeNC, expected_largeNC, 1e-6);
    dis.SetNcScheme(FiniteNC);
    double tripole_finiteNC = dis.TripoleAmplitude(x01, x02, x21, Y);
    double N01 = dis.GetDipole().DipoleAmplitude(x01, Y);

    double expected_finiteNC = NC/(2.0*CF)*((1.0 - N02)*(1.0 - N12) - 1.0/(NC*NC)*(1.0 - N01));
    ASSERT_ALMOST_EQUAL(tripole_finiteNC, expected_finiteNC, 1e-6);
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
