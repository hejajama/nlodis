#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include <dipole/vector.hpp>
#include <dipole/dipoleamplitude.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include "../nlodis.hpp"
#include "../integration.hpp"

#include "dipole/bkdipole/bkdipole.hpp"
#include "dipole/gbw.hpp"

const std::string gbw_datafile = "gbw.dat";

/*
 * LO cross section tests 
 */

 TEST(LEADING_LOG_STRUCTURE_FUNCTIONS)
 {

    BKDipole gbwdatafile(gbw_datafile);
    NLODIS dis;
    dis.SetDipole(std::make_unique<BKDipole>(gbw_datafile));
    Quark u(Quark::Type::U, 0.14); ; u.type = Quark::U; u.mass = 0.14;
    Quark d(Quark::Type::D, 0.14); d.type = Quark::D; d.mass = 0.14; 
    Quark s(Quark::Type::S, 0.14); s.type = Quark::S; s.mass = 0.14; 
    Quark c(Quark::Type::C, 1.4); c.type = Quark::C; c.mass = 1.4;  
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

    // Same, but instead of using u+d+s, use effective light quark datatype
    Quark light(Quark::Type::LIGHT, 0.14);
    dis.SetQuarks({light, c});
    f2 = dis.F2(Q2, 1e-4);
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

TEST(DipoleAmplitude_GBW)
{
    // file gbw.dat corresponds to the GBW dipole amplitude
    // N(r,Y) = 1 - exp(-r^2 Q_s^2(Y)/4)
    // with Q_s^2(Y) = 1.0*exp(lambda*Y) [GeV^2] and lambda=1/3
    BKDipole N(gbw_datafile);


    // Test dipole amplitude values 
    auto gbw_amplitude = [](double r, double Y) {
        double Qs2 = std::exp(Y / 3.0);  // Q_s^2(Y) = exp(-lambda*Y) with lambda=1/3
        return 1.0 - std::exp(-r * r * Qs2 / 4.0);
    };
    ASSERT_ALMOST_EQUAL(N.DipoleAmplitude(0.1, 0.0), gbw_amplitude(0.1, 0.0), gbw_amplitude(0.1, 0.0)/100);
    ASSERT_ALMOST_EQUAL(N.DipoleAmplitude(2.0, 1.0), gbw_amplitude(2.0, 1.0), gbw_amplitude(2.0, 1.0)/100);
    N.InitializeInterpolation(2.05);
    ASSERT_ALMOST_EQUAL(N.DipoleAmplitude(15.0, 2.05), gbw_amplitude(15.0, 2.05), gbw_amplitude(15.0, 2.05)/100);


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


TEST(integration_methods_suave_vegas)
{
    // Compute F2 and FL using both the Suave and Vegas methods, and check that the results are consistent within the expected accuracy of the integration method.
    
    // Avoid especially underflow errors when evaluating K_n(r*eps) 
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    
    NLODIS dis;
    GBWDipole gbw;
    dis.SetDipole(std::make_unique<GBWDipole>());
    Quark light(Quark::Type::LIGHT, 0.14);
    Quark charm(Quark::Type::C, 1.4);
    dis.SetQuarks({light, charm});
    dis.SetProtonTransverseArea(1);
    dis.SetOrder(Order::NLO);
    dis.SetMCIntegrationPoints(6e5);

    const double Q2vals[] = {1,10,50};
    for (auto Q2 : Q2vals) {
        double xbj = 1e-4;
        dis.SetMCIntegrationMethod("suave");
        double FT_suave = dis.FT(Q2, xbj);
        double FL_suave = dis.FL(Q2, xbj);
        dis.SetMCIntegrationMethod("vegas");
        double FT_vegas = dis.FT(Q2, xbj);
        double FL_vegas = dis.FL(Q2, xbj);

        // This test does not require a very good accuracy in order to run fast, 
        // if needed, this can be made stricter if the number of mc integration points is increased also
        // More MC points is needed towards higher Q^2
        const double relacc_goal = 0.15;

        cout << "Q2=" << Q2 << " GeV^2, FL difference " << (FL_suave-FL_vegas)/std::max(std::abs(FL_suave), std::abs(FL_vegas)) <<" FT difference " << (FT_suave-FT_vegas)/std::max(std::abs(FT_suave), std::abs(FT_vegas)) << endl;
        ASSERT_ALMOST_EQUAL(FT_suave, FT_vegas, std::max(std::abs(FT_suave), std::abs(FT_vegas))*relacc_goal); 
        ASSERT_ALMOST_EQUAL(FL_suave, FL_vegas, std::max(std::abs(FL_suave), std::abs(FL_vegas))*relacc_goal);
        
        //
    }

}

TEST_MAIN()

