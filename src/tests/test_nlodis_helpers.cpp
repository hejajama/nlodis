#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include "../nlodis.hpp"
#include <gsl/gsl_sf_dilog.h>

using namespace std;

/*
 * Test EvolutionRapidity function
 * Y = log(W^2 * z2 / Q0^2) where W^2 = Q^2/xbj
 */
TEST(EVOLUTION_RAPIDITY)
{
    NLODIS dis("gbw.dat");
    
    double Q2 = 10.0;
    double xbj = 0.01;
    double z2 = 0.1;
    
    // Expected: Y = log((Q2/xbj) * z2 / Q0^2)
    // Q0^2 = 1 (from nlodis.hpp)
    double W2 = Q2 / xbj;
    double expected_Y = log(W2 * z2 / 1.0);
    
    double Y = dis.EvolutionRapidity(xbj, Q2, z2);
    
    ASSERT_ALMOST_EQUAL(Y, expected_Y, 1e-10);
    
    // Test another case
    Q2 = 5.0;
    xbj = 0.001;
    z2 = 0.5;
    W2 = Q2 / xbj;
    expected_Y = log(W2 * z2 / 1.0);
    Y = dis.EvolutionRapidity(xbj, Q2, z2);
    
    ASSERT_ALMOST_EQUAL(Y, expected_Y, 1e-10);
}

/*
 * Test RunningCouplinScale function
 */
TEST(RUNNING_COUPLING_SCALE)
{
    NLODIS dis("gbw.dat");
    
    double x01 = 1.0;
    double x02 = 2.0;
    double x21 = 1.5;
    
    // Test SMALLEST scheme 
    dis.SetRunningCouplingScheme(SMALLEST);
    double scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x01, 1e-10);
    
    // Test PARENT scheme
    dis.SetRunningCouplingScheme(PARENT);
    scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x01, 1e-10);
    
    // Test with different values
    x01 = 0.5;
    x02 = 3.0;
    x21 = 0.1;
    
    dis.SetRunningCouplingScheme(SMALLEST);
    scale = dis.RunningCouplinScale(x01, x02, x21);
    ASSERT_ALMOST_EQUAL(scale, x21, 1e-10);
    
    dis.SetRunningCouplingScheme(PARENT);
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
    NLODIS dis("gbw.dat");
    double C2=8;
    dis.SetRunningCouplingC2Alpha(C2); 
    
    const double LambdaQCD = 0.241; // GeV
    const int Nf = 4; // number of flavors (u,d,s,c)
    const double b0 = (11.0*NC - 2.0*Nf)/(12.0*M_PI);

    
    // Test at r = 2 GeV^-1 
    double r = 2.0;
    double mu2 = 4.0*C2/(r*r);
    double expected_as = 1.0/(b0*log(mu2/(LambdaQCD*LambdaQCD)));
    double as = dis.Alphas(r);
    ASSERT_ALMOST_EQUAL(as, expected_as, 1e-6);
    
    // Test freezing at large r (should be capped at 0.7)
    r = 100.0;
    as = dis.Alphas(r);
    ASSERT_ALMOST_EQUAL(as, 0.7, 1e-10);
}


/*
 * Test ILdip_massive_Omega_L_Const
 * Check basic properties and some edge cases
 */
TEST(ILDIP_MASSIVE_OMEGA_L_CONST)
{
    double Q2 = 10.0;
    double z = 0.5;
    double r = 2.0;
    double mf = 1.3;
    
    double result = ILdip_massive_Omega_L_Const(Q2, z, r, mf);
    
    // Should be finite
    ASSERT_TRUE(std::isfinite(result));
    
    // Test at very small r (should be approx 0)
    r = 1e-10;
    result = ILdip_massive_Omega_L_Const(Q2, z, r, mf);
    ASSERT_ALMOST_EQUAL(result, 0.0, 1e-15);
    
    // Test symmetry under z -> 1-z
    z = 0.3;
    r = 2.0;
    double result_z = ILdip_massive_Omega_L_Const(Q2, z, r, mf);
    double result_1minusz = ILdip_massive_Omega_L_Const(Q2, 1.0-z, r, mf);
    
    // The function should be symmetric
    ASSERT_ALMOST_EQUAL(result_z, result_1minusz, fabs(result_z)*1e-10);
}


