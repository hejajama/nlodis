#include "unit_test_framework.hpp"
#include <cmath>
#include <algorithm>
#include "../nlodis.hpp"

using namespace std;

// ============================================================
// Transverse dipole helper tests (massive)
//
// Reference formulas: https://arxiv.org/pdf/2204.02486
// (implementation contains some modifications due to numerical stability)
// ============================================================

TEST(OmegaT_V_decomposition_unsymmetric_parts)
{
    // Omega_T^V(z) = Omega_T^V,unsym(z) + Omega_T^V,unsym(1-z)
    // See https://arxiv.org/pdf/2204.02486
    const double Q = sqrt(10.0);
    const double mf = 1.3;
    const double z = 0.31;

    const double full = OmegaT_V(Q, z, mf);
    const double rebuilt = OmegaT_V_unsymmetric(Q, z, mf) + OmegaT_V_unsymmetric(Q, 1.0 - z, mf);

    ASSERT_ALMOST_EQUAL(full, rebuilt, std::max(1e-13, fabs(full) * 1e-11));
}

TEST(OmegaT_N_decomposition_and_antisymmetry)
{
    // Omega_T^N(z) = Omega_T^N,unsym(z) - Omega_T^N,unsym(1-z)
    // and Omega_T^N(1-z) = -Omega_T^N(z).
    // See https://arxiv.org/pdf/2204.02486
    const double Q = sqrt(10.0);
    const double mf = 1.3;
    const double z = 0.27;

    const double full = OmegaT_N(Q, z, mf);
    const double rebuilt = OmegaT_N_unsymmetric(Q, z, mf) - OmegaT_N_unsymmetric(Q, 1.0 - z, mf);
    ASSERT_ALMOST_EQUAL(full, rebuilt, std::max(1e-13, fabs(full) * 1e-11));

    const double anti = OmegaT_N(Q, 1.0 - z, mf);
    ASSERT_ALMOST_EQUAL(full, -anti, std::max(1e-13, fabs(full) * 1e-11));
}

TEST(IT_V1_and_IT_VMS1_decompositions)
{
    // I_T^(V1/VMS1) are defined as sums of unsymmetric pieces in z and 1-z.
    // See https://arxiv.org/pdf/2204.02486
    const double Q = sqrt(8.0);
    const double mf = 1.4;
    const double z = 0.36;
    const double r = 1.2;
    const double xi = 0.41;

    const double v1 = IT_V1(Q, z, mf, r, xi);
    const double v1_rebuilt = IT_V1_unsymmetric(Q, z, mf, r, xi)
        + IT_V1_unsymmetric(Q, 1.0 - z, mf, r, xi);
    ASSERT_ALMOST_EQUAL(v1, v1_rebuilt, std::max(1e-12, fabs(v1) * 1e-10));

    const double vms1 = IT_VMS1(Q, z, mf, r, xi);
    const double vms1_rebuilt = IT_VMS1_unsymmetric(Q, z, mf, r, xi)
        + IT_VMS1_unsymmetric(Q, 1.0 - z, mf, r, xi);
    ASSERT_ALMOST_EQUAL(vms1, vms1_rebuilt, std::max(1e-12, fabs(vms1) * 1e-10));
}

TEST(IT_V2_VMS2_N_decompositions)
{
    // I_T^(V2/VMS2) are symmetric combinations and I_T^N is antisymmetric.
    // See https://arxiv.org/pdf/2204.02486
    const double Q = sqrt(10.0);
    const double mf = 1.3;
    const double z = 0.33;
    const double r = 1.7;
    const double y_chi = 0.43;
    const double y_u = 0.62;

    const double v2 = IT_V2(Q, z, mf, r, y_chi, y_u);
    const double v2_rebuilt = IT_V2_unsymmetric(Q, z, mf, r, y_chi, y_u)
        + IT_V2_unsymmetric(Q, 1.0 - z, mf, r, y_chi, y_u);
    ASSERT_ALMOST_EQUAL(v2, v2_rebuilt, std::max(1e-12, fabs(v2) * 1e-10));

    const double vms2 = IT_VMS2(Q, z, mf, r, y_chi, y_u);
    const double vms2_rebuilt = IT_VMS2_unsymmetric(Q, z, mf, r, y_chi, y_u)
        + IT_VMS2_unsymmetric(Q, 1.0 - z, mf, r, y_chi, y_u);
    ASSERT_ALMOST_EQUAL(vms2, vms2_rebuilt, std::max(1e-12, fabs(vms2) * 1e-10));

    const double n = IT_N(Q, z, mf, r, y_chi, y_u);
    const double n_rebuilt = IT_N_unsymmetric(Q, z, mf, r, y_chi, y_u)
        - IT_N_unsymmetric(Q, 1.0 - z, mf, r, y_chi, y_u);
    ASSERT_ALMOST_EQUAL(n, n_rebuilt, std::max(1e-12, fabs(n) * 1e-10));
}

TEST(ITdip_massive_terms_z_symmetry)
{
    // The transverse dipole ingredients are symmetric under z <-> 1-z.
    // See https://arxiv.org/pdf/2204.02486 (transverse dipole terms)
    const double Q2 = 10.0;
    const double mf = 1.3;
    const double x01sq = 1.8;
    const double xi = 0.37;
    const double y_chi = 0.52;
    const double y_u = 0.58;

    const double zvals[] = {0.21, 0.31, 0.44};
    for (double z : zvals)
    {
        const double zc = 1.0 - z;

        const double t0 = ITdip_massive_0(Q2, z, x01sq, mf);
        const double t0c = ITdip_massive_0(Q2, zc, x01sq, mf);
        ASSERT_ALMOST_EQUAL(t0, t0c, std::max(1e-12, fabs(t0) * 1e-10));

        const double t1 = ITdip_massive_1(Q2, z, x01sq, mf, xi);
        const double t1c = ITdip_massive_1(Q2, zc, x01sq, mf, xi);
        ASSERT_ALMOST_EQUAL(t1, t1c, std::max(1e-12, fabs(t1) * 1e-10));

        const double t2 = ITdip_massive_2(Q2, z, x01sq, mf, y_chi, y_u);
        const double t2c = ITdip_massive_2(Q2, zc, x01sq, mf, y_chi, y_u);
        ASSERT_ALMOST_EQUAL(t2, t2c, std::max(1e-12, fabs(t2) * 1e-10));
    }
}

TEST(IT_V2_unsymmetric_endpoint_finiteness)
{
    // Numerical-stability regression guard: evaluate near y_u endpoints where u=(1-y_u)/y_u can be extreme.
    const double Q = sqrt(12.0);
    const double mf = 1.3;
    const double z = 0.37;
    const double r = 2.1;

    const double y_chi_vals[] = {1e-6, 0.1, 0.5, 0.9, 0.999999};
    const double y_u_vals[] = {1e-6, 1e-4, 1e-2, 0.2, 0.5, 0.9, 0.999, 0.999999};

    for (double y_chi : y_chi_vals)
    for (double y_u : y_u_vals)
    {
        const double v = IT_V2_unsymmetric(Q, z, mf, r, y_chi, y_u);
        ASSERT_TRUE(std::isfinite(v));
    }
}

// ============================================================
// Transverse qqg helper tests (massive)
// ============================================================

TEST(ITNLOqg_I1_quark_antiquark_symmetry)
{
    // Under quark-antiquark exchange: z1 <-> z0=1-z1-z2 and x02^2 <-> x21^2.
    // See https://arxiv.org/pdf/2204.02486
    const double Q2 = 9.0;
    const double mf = 1.3;
    const double z1 = 0.31;
    const double z2 = 0.24;
    const double z0 = 1.0 - z1 - z2;

    const double x01sq = 1.3;
    const double x02sq = 0.6;
    const double x21sq = 0.9;

    const double orig = ITNLOqg_massive_tripole_part_I1(Q2, mf, z1, z2, x01sq, x02sq, x21sq);
    const double perm = ITNLOqg_massive_tripole_part_I1(Q2, mf, z0, z2, x01sq, x21sq, x02sq);

    ASSERT_ALMOST_EQUAL(orig, perm, std::max(1e-11, fabs(orig) * 1e-9));
}

TEST(ITNLOqg_I2_fast_quark_antiquark_symmetry)
{
    // Same exchange symmetry for the fast-I2 ingredient.
    // See https://arxiv.org/pdf/2204.02486
    const double Q2 = 9.0;
    const double mf = 1.3;
    const double z1 = 0.31;
    const double z2 = 0.24;
    const double z0 = 1.0 - z1 - z2;
    const double y_t = 0.42;

    const double x01sq = 1.3;
    const double x02sq = 0.6;
    const double x21sq = 0.9;

    const double orig = ITNLOqg_massive_tripole_part_I2(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t);
    const double perm = ITNLOqg_massive_tripole_part_I2(Q2, mf, z0, z2, x01sq, x21sq, x02sq, y_t);

    ASSERT_ALMOST_EQUAL(orig, perm, std::max(1e-11, fabs(orig) * 1e-8));
}

TEST(ITNLOqg_I3_fast_quark_antiquark_symmetry_equal_yt)
{
    // For y_t1 = y_t2, the fast-I3 ingredient should satisfy the same exchange symmetry.
    // See https://arxiv.org/pdf/2204.02486
    const double Q2 = 9.0;
    const double mf = 1.3;
    const double z1 = 0.31;
    const double z2 = 0.24;
    const double z0 = 1.0 - z1 - z2;
    const double y_t = 0.35;

    const double x01sq = 1.3;
    const double x02sq = 0.6;
    const double x21sq = 0.9;

    const double orig = ITNLOqg_massive_tripole_part_I3(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_t);
    const double perm = ITNLOqg_massive_tripole_part_I3(Q2, mf, z0, z2, x01sq, x21sq, x02sq, y_t, y_t);

    ASSERT_ALMOST_EQUAL(orig, perm, std::max(1e-11, fabs(orig) * 1e-8));
}

TEST(ITNLOqg_dipole_uvsub_large_daughter_decay)
{
    // UV-subtraction terms include exp(-x2^2/(x01^2 euler)) suppression.
    // At very large daughter dipoles this tends to zero.
    // See https://arxiv.org/pdf/2204.02486
    const double Q2 = 10.0;
    const double mf = 1.3;
    const double z1 = 0.3;
    const double z2 = 0.2;

    const double x01sq = 0.2;
    const double large = ITNLOqg_massive_dipole_uvsub(Q2, mf, z1, z2, x01sq, 100.0, 100.0);

    ASSERT_ALMOST_EQUAL(large, 0.0, 1e-30);
}
