
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
    double xi=0.3;
    double x=0.1;

    double Icd = ILdip_massive_Icd(Q2, z, r, mf, xi, x);

    double expected_Icd =2.9343444440607052e-05;

    
    ASSERT_ALMOST_EQUAL(Icd, expected_Icd, std::min(Icd, expected_Icd)/1e2);
}