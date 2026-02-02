
#include "unit_test_framework.hpp"
#include <cmath>
#include <string>
#include "../nlodis.hpp"
#include "../integration.hpp"
#include <cuba.h>

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