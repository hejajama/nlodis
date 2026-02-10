#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include <dipole/vector.hpp>
#include <dipole/dipoleamplitude.hpp>
#include <gsl/gsl_rng.h>
#include "../nlodis.hpp"
#include "../integration.hpp"

#include "dipole/bkdipole/bkdipole.hpp"

const std::string gbw_datafile = "gbw.dat";

/*
 * LO cross section tests 
 */

 TEST(LEADING_LOG_STRUCTURE_FUNCTIONS)
 {

    BKDipole gbwdatafile(gbw_datafile);
    NLODIS dis;
    dis.SetDipole(std::make_unique<BKDipole>(gbw_datafile));
    Quark u; u.type = Quark::U; u.mass = 0.14; u.charge = 2.0/3.0;
    Quark d; d.type = Quark::D; d.mass = 0.14; d.charge = -1.0/3.0;
    Quark s; s.type = Quark::S; s.mass = 0.14;  s.charge = -1.0/3.0;
    Quark c; c.type = Quark::C; c.mass = 1.4;   c.charge = 2.0/3.0;
    std::vector<Quark> quark_list = {u,d,s,c};
    dis.SetQuarks(quark_list);

    // Set \sigma_0/2 = 1 so that results match the setup used to compute reference values
    dis.SetProtonTransverseArea(1);

    double Q2 = 10.0; // GeV^2
    double xbj = 1;  // Test case datafile is generated such that it starts from x0=1
    
    dis.SetOrder(Order::LO);

    

    double f2_ic=dis.F2(Q2,xbj);
    // Reference: factor 2 is the optical theorem 2 not included when the refernece value
    // was computed
    double ref = (0.0202973+0.00543532)*2; // light + charm + bottom
    
    ASSERT_ALMOST_EQUAL(f2_ic, ref, ref/100); // light + charm

    // Smaller x, tests interpolation and computation of the evolution rapidity
    double f2 = dis.F2(Q2, 1e-4);
    ref = 2.0*(0.120141+0.0538897);
    ASSERT_ALMOST_EQUAL(f2, ref, ref/100); // light +charm
 }

/*
 * Tripole i.e. 1-S_{012}
 * S_{012} defined in https://arxiv.org/pdf/2007.01645 eq 13
 */
TEST(TRIPOLE_AMPLITUDE)
{
    NLODIS dis;
    BKDipole gbwdatafile("gbw.dat");
    dis.SetDipole(std::make_unique<BKDipole>(gbwdatafile));

    dis.SetNcScheme(NcScheme::LargeNC);
    double x01 = 1.0;
    double x02 = 2.0;
    double x21 = std::sqrt( SQR(x01) + SQR(x02) - 2.0*x01*x02*0.5 ); // angle cos=0.5
    double Y = 2.0;

    double tripole_largeNC = dis.TripoleAmplitude(x01, x02, x21, Y);
    double N02 = dis.GetDipole().DipoleAmplitude(x02, Y);
    double N12 = dis.GetDipole().DipoleAmplitude(x21, Y);
    double expected_largeNC = 1.0 - (1.0 - N02)*(1.0 - N12);
    ASSERT_ALMOST_EQUAL(tripole_largeNC, expected_largeNC, 1e-6);
    dis.SetNcScheme(NcScheme::FiniteNC);
    double tripole_finiteNC = dis.TripoleAmplitude(x01, x02, x21, Y);
    double N01 = dis.GetDipole().DipoleAmplitude(x01, Y);

    double expected_finiteNC = 1.-Constants::NC/(2.0*Constants::CF)*((1.0 - N02)*(1.0 - N12) - 1.0/(Constants::NC*Constants::NC)*(1.0 - N01));
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
    BKDipole N(gbw_datafile);

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

TEST(vector_class) {
    gsl_rng *global_rng = gsl_rng_alloc(gsl_rng_default);
    Vec v1(1,2);
    Vec v2(5,-1);
    double eps=1e-7;
    ASSERT_ALMOST_EQUAL(v1*v2, 3,eps);
    ASSERT_ALMOST_EQUAL((v1+v2).GetX(), 6, eps)
    ASSERT_ALMOST_EQUAL((v1+v2*(-4)).GetY(), 2+(-1)*(-4), eps)
    ASSERT_ALMOST_EQUAL(v1.Len(), std::sqrt(1*1+2*2),eps);



    gsl_rng_free(global_rng);
}

TEST_MAIN()
