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

/**
 * NLO cross section tests
 * Reference numbers are computed using the code published in April 2026
 */
TEST(NLO_structure_functions)
{
    Quark light(Quark::Type::LIGHT, 0.14);
    Quark charm(Quark::Type::C, 1.4);
    std::vector<Quark> quark_list = {light, charm};
    NLODIS dis;
    dis.SetQuarks(quark_list);
    dis.SetDipole(std::make_unique<GBWDipole>(1,0.3, 1, 1)); // (double Qs0sqr, double lambda, double gamma, double X0);
    dis.SetMCIntegrationPoints(5e6);
    dis.SetMCIntegrationEpsRel(1e-3);
    dis.SetOrder(Order::NLO);
    dis.SetProtonTransverseArea(10, Unit::MB); // Set \sigma_0/2
    dis.SetRunningCouplingC2(2);
    dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::SMOOTH);
    dis.SetActiveFlavors(4);
    dis.SetNcScheme(NcScheme::FiniteNC);

    double Q2 = 10.0; // GeV^2
    double xbj = 1e-4;  

    double FL = dis.FL(Q2, xbj);
    double FT = dis.FT(Q2, xbj);
    double expected_FL = 1.22661;
    double expected_FT = 5.72033;

    ASSERT_ALMOST_EQUAL(FL, expected_FL, std::abs(expected_FL/1e2));    
    ASSERT_ALMOST_EQUAL(FT, expected_FT, std::abs(expected_FT/1e2));



}

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
    BKDipole gbw(gbw_datafile);
    dis.SetDipole(std::make_unique<BKDipole>(gbw));

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

    // NLO structure functions

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

        //cout << "Q2=" << Q2 << " GeV^2, FL difference " << (FL_suave-FL_vegas)/std::max(std::abs(FL_suave), std::abs(FL_vegas)) <<" FT difference " << (FT_suave-FT_vegas)/std::max(std::abs(FT_suave), std::abs(FT_vegas)) << endl;
        ASSERT_ALMOST_EQUAL(FT_suave, FT_vegas, std::max(std::abs(FT_suave), std::abs(FT_vegas))*relacc_goal); 
        ASSERT_ALMOST_EQUAL(FL_suave, FL_vegas, std::max(std::abs(FL_suave), std::abs(FL_vegas))*relacc_goal);
        
    }


    // Test integration method used to compute 2 dimensional integrals
    // Makes sure that the cuhre algorithm is used correclty by evaluating the same result
    // using also the standard vegas method
    // These are evaluated in NLODIS::Sigma_dip_d2b
    for (auto Q2 : Q2vals) {
        IntegrationParams intparams{dis};
        intparams.Q2=Q2;
        intparams.xbj=0.001;
        intparams.pol=Polarization::L;

        // Contribution from 
        intparams.quark=light;
        double I, Ierr, Iprob;
        intparams.contribution="Omega_L_const";
        Cuba("cuhre", 2, integrand_dip_massive, &intparams, &I, &Ierr, &Iprob, dis.GetMCIntegrationConfig());
        double I_vegas, Ierr_vegas, Iprob_vegas;
        Cuba("vegas", 2, integrand_dip_massive, &intparams, &I_vegas, &Ierr_vegas, &Iprob_vegas, dis.GetMCIntegrationConfig());
        //cout << "Dipole contribution integrand difference " << (I-I_vegas)/std::max(std::abs(I), std::abs(I_vegas)) << endl;
        ASSERT_ALMOST_EQUAL(I, I_vegas, std::max(std::abs(I), std::abs(I_vegas))*0.01);

        intparams.pol=Polarization::T;
        // Contribution from 
        intparams.quark=charm;
        intparams.contribution="T0";
        Cuba("cuhre", 2, integrand_dip_massive, &intparams, &I, &Ierr, &Iprob, dis.GetMCIntegrationConfig());
        Cuba("vegas", 2, integrand_dip_massive, &intparams, &I_vegas, &Ierr_vegas, &Iprob_vegas, dis.GetMCIntegrationConfig());
        //cout << "Dipole contribution integrand difference " << (I-I_vegas)/std::max(std::abs(I), std::abs(I_vegas)) << endl;
        ASSERT_ALMOST_EQUAL(I, I_vegas, std::max(std::abs(I), std::abs(I_vegas))*0.01);
    }
}



/*
 * Test EvolutionRapidity function
 * Y = log(W^2 * z2 / Q0^2) where W^2 = Q^2/xbj
 */
TEST(EVOLUTION_RAPIDITY)
{
    NLODIS dis;
    
    double Q2 = 10.0;
    double xbj = 0.01;
    double z2 = 0.1;
    
    // Expected: Y = log((Q2/xbj) * z2 / Q0^2)
    // Q0^2 = 1 (from nlodis.hpp)
    double W2 = Q2 / xbj;
    double expected_Y = log(W2 * z2 / 1.0);
    
    double Y = dis.EvolutionRapidity_qqg(xbj, Q2, z2);
    
    ASSERT_ALMOST_EQUAL(Y, expected_Y, 1e-10);
    
    // Test another case
    Q2 = 5.0;
    xbj = 0.001;
    z2 = 0.5;
    W2 = Q2 / xbj;
    expected_Y = log(W2 * z2 / 1.0);
    Y = dis.EvolutionRapidity_qqg(xbj, Q2, z2);
    
    ASSERT_ALMOST_EQUAL(Y, expected_Y, 1e-10);
}

/*
 * Test RunningCouplinScale function
 */
TEST(RUNNING_COUPLING_SCALE)
{
    NLODIS dis;
    
    double x01 = 1.0;
    double x02 = 2.0;
    double x21 = 1.5;
    
    // Test SMALLEST scheme 
    dis.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    double scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x01, 1e-10);
    
    // Test PARENT scheme
    dis.SetRunningCouplingScheme(RunningCouplingScheme::PARENT);
    scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x01, 1e-10);
    
    // Test with different values
    x01 = 0.5;
    x02 = 3.0;
    x21 = 0.1;

    dis.SetRunningCouplingScheme(RunningCouplingScheme::SMALLEST);
    scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x21, 1e-10);
    
    dis.SetRunningCouplingScheme(RunningCouplingScheme::PARENT);
    scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x01, 1e-10);

   
}

/*
 * Test Alphas (coordinate space coupling)
 * Formula: alphas = 1/(b0 * log(mu^2/Lambda^2))
 * where mu^2 = 4C^2/r^2 + Lambda^2
 */
TEST(ALPHAS_COORDINATE_SPACE)
{
    NLODIS dis;
    double C2=8;
    dis.SetRunningCouplingC2(C2); 
    
    const double LambdaQCD = 0.241; // GeV
    const int Nf = 4; // number of flavors (u,d,s,c)
    const double b0 = (11.0*Constants::NC - 2.0*Nf)/(12.0*M_PI);

    
    // Test at r = 2 GeV^-1 
    double r = 2.0;
    double mu2 = 4.0*C2/(r*r);
    double expected_as = 1.0/(b0*log(mu2/(LambdaQCD*LambdaQCD)));
    dis.SetRunningCouplingIRScheme(RunningCouplingIRScheme::FREEZE);
    double as = dis.Alphas(r);
    ASSERT_ALMOST_EQUAL(as, expected_as, 1e-6);
    
    // Test freezing at large r (should be capped at 0.7)
    r = 100.0;
    as = dis.Alphas(r);
    ASSERT_ALMOST_EQUAL(as, 0.7, 1e-10);

    // Test user-defined freeze cap
    dis.SetRunningCouplingMaxAlphaS(0.5);
    as = dis.Alphas(r);
    ASSERT_ALMOST_EQUAL(as, 0.5, 1e-10);
}



TEST_MAIN()



