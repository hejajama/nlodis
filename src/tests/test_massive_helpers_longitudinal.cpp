
#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include "../nlodis.hpp"
#include "../integration.hpp"
#include <cuba.h>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

struct Inthleperams
{
    double r;
    double z;
    double mf;

    IntegrationParams userdata;
};

const string cubamethod = "cuhre";

TEST(I_cd_massive_L)
{
    double r=2;
    double z=0.3;
    double Q2=10;
    double mf=1.3;
    double xi=0.6;
    double x=0.1;

    

    double Icd = ILdip_massive_Icd(Q2, z, r, mf, xi, x);

    double expected_Icd =-5.986346671214665e-05;

    
    ASSERT_ALMOST_EQUAL(Icd, expected_Icd, std::abs(std::min(Icd, expected_Icd)/1e2));

    double Iab = ILdip_massive_Iab(Q2, z, r, mf, xi);
    double expected_Iab = 0.0014655116231699;
    ASSERT_ALMOST_EQUAL(Iab, expected_Iab, std::abs(std::min(Iab, expected_Iab)/1e2));
}

/* \mathcal{G}
 * One integration done analytically, see (69) of the note
 */
TEST(Gintegrated)
{
    int a = 2, b = 1;
    double Qbar = std::sqrt(2.0), mf = 1.3, x2 = std::sqrt(0.4), x3 = std::sqrt(0.5);
    double omega = 0.7, lambda = 1.2;

    struct GParams {
        int a, b;
        double Qbar, mf, x2, x3, omega, lambda;
    } params = {a, b, Qbar, mf, x2, x3, omega, lambda};

    gsl_function F;
    F.function = [](double y, void* p) -> double {
        auto* gp = static_cast<GParams*>(p);
        return G_integrand_simplified(gp->a, gp->b, gp->Qbar, gp->mf, gp->x2, gp->x3, gp->omega, gp->lambda, y);
    };
    F.params = &params;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_integration_qags(&F, 0.0, 1.0, 0, 1e-4, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    double expected = 1.44491135902142;

    ASSERT_ALMOST_EQUAL(result, expected, std::min(result, expected)/1e2);

}

// ============================================================
// ILdip_massive_Iab tests  (2103.14549 eq. 113)
// ============================================================

/*
 * Iab z <-> 1-z symmetry.
 * The front factor 4 Q^2 (z(1-z))^2 is manifestly symmetric.
 * Under z -> 1-z: bessel_arg_2 and bessel_arg_3 are exchanged so
 * (2 b1 - b2 - b3) is unchanged, and the kernel
 *   1/xi * (-2 log(xi)/(1-xi) + (1+xi)/2)
 * does not depend on z.  Hence the integrand is symmetric for all xi.
 */
TEST(ILdip_massive_Iab_z_symmetry)
{
    const double Q2 = 10.0;
    const double r   = 1.5;
    const double mf  = 1.3;

    const double z_values[]  = {0.2, 0.3, 0.4};
    const double xi_values[] = {0.1, 0.3, 0.5, 0.7, 0.9};

    for (double z : z_values)
        for (double xi : xi_values)
        {
            double val_z   = ILdip_massive_Iab(Q2, z,       r, mf, xi);
            double val_1mz = ILdip_massive_Iab(Q2, 1.0 - z, r, mf, xi);
            ASSERT_ALMOST_EQUAL(val_z, val_1mz, fabs(val_z) * 1e-10);
        }
}

/*
 * At z = 0.5 the two sub-leading Bessel arguments are equal (bessel_arg_2 = bessel_arg_3),
 * so the combination (2 b1 - b2 - b3) simplifies to 2(b1 - b2).
 * We verify this equality of arguments analytically and then check that the
 * function evaluates to the simplified closed form.
 */
TEST(ILdip_massive_Iab_z_half_b2_equals_b3)
{
    const double Q2 = 10.0;
    const double r   = 1.5;
    const double mf  = 1.3;
    const double z   = 0.5;
    const double xi  = 0.3;

    double kappa_z = sqrt(z * (1.0 - z) * Q2 + SQR(mf));

    // bessel_arg_2 = sqrt(kappa_z^2 + (1-z)*xi/(1-xi)*mf^2)*r
    // bessel_arg_3 = sqrt(kappa_z^2 +     z *xi/(1-xi)*mf^2)*r
    // At z=0.5 both coefficients of mf^2 are identical -> args must be equal
    double arg2 = sqrt(SQR(kappa_z) + (1.0 - z) * xi / (1.0 - xi) * SQR(mf)) * r;
    double arg3 = sqrt(SQR(kappa_z) +        z  * xi / (1.0 - xi) * SQR(mf)) * r;
    ASSERT_ALMOST_EQUAL(arg2, arg3, 1e-14);

    // Simplified form:  front * b1 * kernel * 2*(b1 - b2)
    double front   = 4.0 * Q2 * SQR(z * (1.0 - z));
    double b1      = gsl_sf_bessel_K0(kappa_z * r);
    double b2      = gsl_sf_bessel_K0(arg2);
    double kernel  = 1.0 / xi * (-2.0 * log(xi) / (1.0 - xi) + (1.0 + xi) / 2.0);
    double expected = front * b1 * kernel * 2.0 * (b1 - b2);

    double computed = ILdip_massive_Iab(Q2, z, r, mf, xi);
    ASSERT_ALMOST_EQUAL(computed, expected, fabs(expected) * 1e-10);
}

/*
 * Iab positivity.
 * Physical reasoning:
 *  - Front factor 4 Q^2 (z(1-z))^2 >= 0.
 *  - kernel = 1/xi*(-2 log xi/(1-xi) + (1+xi)/2) > 0 for xi in (0,1)
 *    because log(xi) < 0 so the first term is positive.
 *  - K0 is monotonically decreasing; bessel_arg_2/3 > kappa_z*r,
 *    so b2, b3 < b1  =>  (2 b1 - b2 - b3) > 0.
 * Therefore Iab >= 0 for all physical Q2, z, r, mf, xi.
 */
TEST(ILdip_massive_Iab_positivity)
{
    const double Q2 = 10.0;
    const double r   = 1.5;
    const double mf  = 1.3;

    const double z_vals[]  = {0.2, 0.3, 0.5, 0.7, 0.8};
    const double xi_vals[] = {0.1, 0.3, 0.5, 0.7, 0.9};

    for (double z : z_vals)
        for (double xi : xi_vals)
        {
            double result = ILdip_massive_Iab(Q2, z, r, mf, xi);
            ASSERT_TRUE(result >= 0.0);
        }
}

// ============================================================
// ILdip_massive_Icd tests  (2103.14549 eq. 114, 115)
// ============================================================

/*
 * Icd z <-> 1-z symmetry.
 * Under z -> 1-z: C^L_m(z, x, xi) and C^L_m(1-z, x, xi) are exchanged,
 * and simultaneously kappa(z,x,xi) and kappa(1-z,x,xi) are exchanged.
 * The front factor (z(1-z))^2 is symmetric.  The two Bessel-K0 difference
 * terms therefore swap, leaving the sum invariant.
 */
TEST(ILdip_massive_Icd_z_symmetry)
{
    const double Q2 = 10.0;
    const double r   = 1.5;
    const double mf  = 1.3;
    const double xi  = 0.4;
    const double x   = 0.3;

    const double z_vals[] = {0.2, 0.3, 0.4};
    for (double z : z_vals)
    {
        double val_z   = ILdip_massive_Icd(Q2, z,       r, mf, xi, x);
        double val_1mz = ILdip_massive_Icd(Q2, 1.0 - z, r, mf, xi, x);
        ASSERT_ALMOST_EQUAL(val_z, val_1mz, fabs(val_z) * 1e-10);
    }
}

/*
 * Icd finiteness at a range of kinematic points.
 */
TEST(ILdip_massive_Icd_finiteness)
{
    const double Q2    = 10.0;
    const double mf    = 1.3;

    const double r_vals[]  = {0.5, 1.0, 2.0, 5.0};
    const double z_vals[]  = {0.2, 0.5, 0.8};
    const double xi_vals[] = {0.1, 0.5, 0.9};
    const double x_vals[]  = {0.1, 0.5, 0.9};

    for (double r  : r_vals)
    for (double z  : z_vals)
    for (double xi : xi_vals)
    for (double x  : x_vals)
    {
        double result = ILdip_massive_Icd(Q2, z, r, mf, xi, x);
        ASSERT_TRUE(std::isfinite(result));
    }
}




// ============================================================
// L_dip tests  (2103.14549 eq. 101)
// ============================================================

/*
 * L_dip massless limit: in the massless case gamma -> 1. There would be
 * a divergence that is subtracted from L_dip.
 * Hence in the massless limit L_dip should vanish
 */
TEST(L_dip_massless_limit_multiple_z)
{
    const double Q2 = 1000.0; // large Q^2 -> deep massless region
    const double mf  = 1e-5;  // effectively massless

    const double z_values[] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9};
    for (double z : z_values)
    {
        double result   = L_dip(Q2, z, mf);
        double expected = 0.0; // L_dip should vanish in the massless limit
        ASSERT_ALMOST_EQUAL(result, expected, 1e-6);
    }
}

/*
 * L_dip z <-> 1-z symmetry: the four dilogarithm arguments come in pairs
 * that are exchanged under z -> 1-z, so L_dip is symmetric.
 */
TEST(L_dip_z_symmetry)
{
    const double Q2 = 10.0;
    const double mf  = 1.3;

    const double z_values[] = {0.1, 0.2, 0.3, 0.4};
    for (double z : z_values)
    {
        double val_z   = L_dip(Q2, z,       mf);
        double val_1mz = L_dip(Q2, 1.0 - z, mf);
        ASSERT_ALMOST_EQUAL(val_z, val_1mz, fabs(val_z) * 1e-10);
    }
}

TEST(Transverse_I_Omega_L_massless_limits)
{
    /// Massless limit test
    // https://arxiv.org/pdf/2204.02486 eq (113)


    double Q=std::sqrt(10);
    double z=0.3;
    double mf=0.001; // almost massless

    ASSERT_ALMOST_EQUAL(OmegaT_V(Q,z,mf), 0, 1e-5);
    
    

}



// ============================================================
// OmegaL_V tests  (2103.14549 eq. 100)
// ============================================================

/*
 * OmegaL_V massless limit: eq. (100) in 2103.14549 implies OmegaL_V -> 0
 * as mf -> 0 (gamma -> 1).  In that limit every log-term cancels pairwise:
 *   1/(2z)[log(1-z) + 1*log(2/(2-2z))] + 1/(2(1-z))[log(z) + 1*log(2/(2z))]
 *   = 1/(2z)*0 + 1/(2(1-z))*0  = 0
 */
TEST(OmegaL_V_massless_limit)
{
    const double Q2 = 10.0;
    const double mf = 1e-6; // effectively massless

    const double z_values[] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9};
    for (double z : z_values)
    {
        double result = OmegaL_V(Q2, z, mf);
        // Relative tolerance loosened because the limit converges slowly in gamma
        ASSERT_ALMOST_EQUAL(result, 0.0, 5e-4);
    }
}

/*
 * OmegaL_V z <-> 1-z symmetry
 * From eq. (100), the two log-blocks are exchanged under z -> 1-z so the sum
 * is invariant.
 */
TEST(OmegaL_V_z_symmetry)
{
    const double Q2 = 10.0;
    const double mf  = 1.3; // charm quark mass

    const double z_values[] = {0.1, 0.2, 0.3, 0.4};
    for (double z : z_values)
    {
        double val_z    = OmegaL_V(Q2, z,       mf);
        double val_1mz  = OmegaL_V(Q2, 1.0 - z, mf);
        ASSERT_ALMOST_EQUAL(val_z, val_1mz, fabs(val_z) * 1e-10);
    }
}



// ============================================================
// ILdip_massive_Omega_L_Const test  (2103.14549 eq. 166, 2nd line)
// ============================================================

/*
 * Exponential decay at large r.
 * ILdip_massive_Omega_L_Const contains K0(kappa_z * r)^2.  Since
 * K0(x) ~ sqrt(pi/(2x)) * exp(-x) for large x, the integrand decreases
 * super-exponentially, so the value at r=30 must be smaller than at r=20.
 */
TEST(ILdip_massive_Omega_L_Const_large_r_decay)
{
    const double Q2 = 4.0;
    const double z  = 0.5;
    const double mf = 1.3;

    double val_r20 = ILdip_massive_Omega_L_Const(Q2, z, 20.0, mf);
    double val_r30 = ILdip_massive_Omega_L_Const(Q2, z, 30.0, mf);
    ASSERT_TRUE(fabs(val_r30) < fabs(val_r20));
}

/*
 * Omega_L_Const is symmetric under z <-> 1-z (the prefactor (z(1-z))^2
 * and kappa_z^2 = z(1-z)Q^2 + mf^2 are both symmetric, and the remaining
 * kinematic function OmegaL_V + L_dip is symmetric by the tests above).
 */
TEST(ILdip_massive_Omega_L_Const_z_symmetry)
{
    const double Q2 = 10.0;
    const double r   = 2.0;
    const double mf  = 1.3;

    const double z_vals[] = {0.1, 0.2, 0.3, 0.4};
    for (double z : z_vals)
    {
        double val_z   = ILdip_massive_Omega_L_Const(Q2, z,       r, mf);
        double val_1mz = ILdip_massive_Omega_L_Const(Q2, 1.0 - z, r, mf);
        ASSERT_ALMOST_EQUAL(val_z, val_1mz, fabs(val_z) * 1e-10);
    }
}

// ============================================================
// G_integrand_simplified tests  (2103.14549 eq. 163 / note eq. 69)
// ============================================================

/*
 * Analytical verification for (a=1, b=1).
 *
 * With a=1, b=1 the formula reduces to:
 *   G^(1,1)(y) = (2/y) * K_0( sqrt((y*lambda*mf^2 + Qbar^2 + mf^2)*(y*x3^2 + omega*x2^2)/y) )
 *
 * because:
 *   y^(-0.5*(2-a+b)) = y^(-1)
 *   2^(a+b-1)        = 2^1  = 2
 *   omega^(b-1)      = omega^0 = 1
 *   fraction^(0.5*(a+b-2)) = fraction^0 = 1
 *   K_{a+b-2}        = K_0
 *
 * We evaluate this independently and compare with the function.
 */
TEST(G_integrand_simplified_a1b1_analytical)
{
    const int    a      = 1, b = 1;
    const double Qbar   = sqrt(2.0);
    const double mf     = 1.3;
    const double x2     = sqrt(0.4);
    const double x3     = sqrt(0.5);
    const double omega  = 0.7;
    const double lambda = 1.2;
    const double y      = 0.5;

    double num   = y * lambda * SQR(mf) + SQR(Qbar) + SQR(mf); // = 4.704
    double inner = y * SQR(x3) + omega * SQR(x2);               // = 0.53
    double expected = (2.0 / y) * gsl_sf_bessel_K0(sqrt(num * inner / y));

    double computed = G_integrand_simplified(a, b, Qbar, mf, x2, x3, omega, lambda, y);
    ASSERT_ALMOST_EQUAL(computed, expected, fabs(expected) * 1e-10);
}

/*
 * G_integrand_simplified is non-negative for all y in (0,1).
 * K_nu(x) > 0 for x > 0, all nu, and all prefactors are positive.
 */
TEST(G_integrand_simplified_positivity)
{
    const double Qbar   = sqrt(3.0);
    const double mf     = 0.14;
    const double x2     = sqrt(0.5);
    const double x3     = sqrt(0.5);
    const double omega  = 0.5;
    const double lambda = 0.5;

    // Test (a=1, b=1) and (a=1, b=2) cases
    const int ab_pairs[2][2] = {{1, 1}, {1, 2}};
    const double y_vals[] = {0.01, 0.1, 0.5, 0.9};

    for (auto& ab : ab_pairs)
        for (double y : y_vals)
        {
            double result = G_integrand_simplified(ab[0], ab[1], Qbar, mf, x2, x3, omega, lambda, y);
            ASSERT_TRUE(std::isfinite(result));
            ASSERT_TRUE(result >= 0.0);
        }
}

// ============================================================
// ILNLOqg permutation symmetry (quark-antiquark exchange)
// ============================================================
//
// Physical origin: the longitudinal photon wavefunction is symmetric under
// exchange of the quark and antiquark.  In the qqg impact factors this means
// z1 <-> z0 = 1-z1-z2  combined with  x02^2 <-> x21^2
// (the gluon position labels re-assign accordingly).
//
// Note: x20x21 = -1/2*(x01^2 - x21^2 - x02^2) is symmetric under x02<->x21,
//       so it does not change sign.
//
// All three parts (I1, I2, I3 at y_t1=y_t2) of sigma_qg^L must satisfy this.

/*
 * I1 permutation symmetry.
 */
TEST(ILNLOqg_I1_quark_antiquark_symmetry)
{
    const double Q2 = 10.0;
    const double mf  = 1.3;
    const double z1  = 0.3;
    const double z2  = 0.2;
    const double z0  = 1.0 - z1 - z2;

    const double x01sq = 1.0;
    const double x02sq = 0.6;
    const double x21sq = 0.8;

    double I1_orig = ILNLOqg_massive_tripole_part_I1(Q2, mf, z1, z2, x01sq, x02sq, x21sq);
    double I1_perm = ILNLOqg_massive_tripole_part_I1(Q2, mf, z0, z2, x01sq, x21sq, x02sq);

    ASSERT_ALMOST_EQUAL(I1_orig, I1_perm, fabs(I1_orig) * 1e-10);
}

/*
 * I1 permutation symmetry: second set of kinematics (charm quark, wider dipoles).
 */
TEST(ILNLOqg_I1_quark_antiquark_symmetry_charm)
{
    const double Q2 = 4.0;    // Q ~ 2 GeV, relevant for charm production
    const double mf  = 1.4;
    const double z1  = 0.4;
    const double z2  = 0.15;
    const double z0  = 1.0 - z1 - z2;

    const double x01sq = 0.5;
    const double x02sq = 1.2;
    const double x21sq = 0.9;

    double I1_orig = ILNLOqg_massive_tripole_part_I1(Q2, mf, z1, z2, x01sq, x02sq, x21sq);
    double I1_perm = ILNLOqg_massive_tripole_part_I1(Q2, mf, z0, z2, x01sq, x21sq, x02sq);

    ASSERT_ALMOST_EQUAL(I1_orig, I1_perm, fabs(I1_orig) * 1e-10);
}

/*
 * I2 permutation symmetry.
 * The single additional integral variable y_t is shared by both the k- and
 * l-contributions, so it is unaffected by the quark-antiquark exchange.
 */
TEST(ILNLOqg_I2_quark_antiquark_symmetry)
{
    const double Q2  = 10.0;
    const double mf   = 1.3;
    const double z1   = 0.3;
    const double z2   = 0.2;
    const double z0   = 1.0 - z1 - z2;
    const double y_t  = 0.5;

    const double x01sq = 1.0;
    const double x02sq = 0.5;
    const double x21sq = 0.7;

    double I2_orig = ILNLOqg_massive_tripole_part_I2(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t);
    double I2_perm = ILNLOqg_massive_tripole_part_I2(Q2, mf, z0, z2, x01sq, x21sq, x02sq, y_t);

    ASSERT_ALMOST_EQUAL(I2_orig, I2_perm, fabs(I2_orig) * 1e-8);
}

/*
 * I3 permutation symmetry at y_t1 = y_t2.
 * When the two integration variables coincide, all cross-terms int_k * int_l
 * become symmetric (commutativity), so the full I3 is invariant under
 * z1 <-> z0 combined with x02^2 <-> x21^2.
 */
TEST(ILNLOqg_I3_quark_antiquark_symmetry_equal_yt)
{
    const double Q2   = 10.0;
    const double mf    = 1.3;
    const double z1    = 0.3;
    const double z2    = 0.2;
    const double z0    = 1.0 - z1 - z2;
    const double y_t   = 0.4;   // y_t1 = y_t2

    const double x01sq = 1.0;
    const double x02sq = 0.5;
    const double x21sq = 0.7;

    double I3_orig = ILNLOqg_massive_tripole_part_I3(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_t);
    double I3_perm = ILNLOqg_massive_tripole_part_I3(Q2, mf, z0, z2, x01sq, x21sq, x02sq, y_t, y_t);

    ASSERT_ALMOST_EQUAL(I3_orig, I3_perm, fabs(I3_orig) * 1e-8);
}

// ============================================================
// ILNLOqg_massive_dipole_uvsub sign test
// ============================================================

/*
 * The UV-subtraction term for the qg dipole contribution is defined with
 * explicit minus signs (cf. code terms -1/x02^2 * ... * K0^2,
 * -1/x21^2 * ... * K0^2).  All other factors (exponential suppression,
 * K0^2, z^2, etc.) are non-negative, so the total must be <= 0.
 * Physically, this term subtracts the UV-singular part to render the
 * remaining integral finite.
 */
TEST(ILNLOqg_dipole_uvsub_is_nonpositive)
{
    const double Q2 = 10.0;
    const double mf  = 1.3;
    const double z1  = 0.3;
    const double z2  = 0.1;

    // Small x02/x21 relative to x01 so exp(-x/x01/e) ~ 1
    const double x01sq = 1.0;
    const double x02sq = 0.02;
    const double x21sq = 0.02;

    double result = ILNLOqg_massive_dipole_uvsub(Q2, mf, z1, z2, x01sq, x02sq, x21sq);
    ASSERT_TRUE(result <= 0.0);
}

/*
 * UV-subtraction decays when x02, x21 >> x01 because the exponential factor
 * exp(-x02/x01/e) -> 0.  In that limit ILNLOqg_massive_dipole_uvsub -> 0.
 */
TEST(ILNLOqg_dipole_uvsub_large_x02_decay)
{
    const double Q2 = 10.0;
    const double mf  = 1.3;
    const double z1  = 0.3;
    const double z2  = 0.1;
    const double x01sq = 0.1;

    // For very large x02 the exp suppression kills the contribution
    double uvsub_large = ILNLOqg_massive_dipole_uvsub(Q2, mf, z1, z2, x01sq, 100.0, 100.0);
    ASSERT_ALMOST_EQUAL(uvsub_large, 0.0, 1e-30);
}
