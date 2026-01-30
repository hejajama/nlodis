
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