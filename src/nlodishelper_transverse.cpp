#include "nlodis.hpp"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_bessel.h>

// Various helper functions defined in this file 
// after they are referenced 
// See note docs/NLO_DIS_cross_section_with_massive_quarks.pdf for definitions
double OmegaT_V_unsymmetric( double Q, double z, double mf );
double L_dip( double Q2, double z, double mf ); // defined in nlodishleper_longitudinal.cpp
double OmegaT_N_unsymmetric( double Q, double z, double mf );
double IT_V1_unsymmetric( double Q, double z, double mf, double r, double xi );
double IT_VMS1_unsymmetric( double Q, double z, double mf, double r, double xi );
double IT_V2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );
double IT_VMS2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );
double IT_N_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );
double IT_V1(double Q, double z, double mf, double r, double xi);
double IT_VMS1(double Q, double z, double mf, double r, double xi);
double IT_V2(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_VMS2(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_N(double Q, double z, double mf, double r, double y_chi, double y_u);

// Eq (41) in the note docs/NLO_DIS_cross_section_with_massive_quarks.pdf
double ITdip_massive_0(double Q2, double z1, double x01sq, double mf) {

    double x01 = sqrt( x01sq );
    double Q = sqrt(Q2);

    double kappa_z = sqrt( z1*(1.0-z1)*Q2 + SQR(mf) );
    double term1 = SQR( kappa_z * gsl_sf_bessel_K1( x01 *kappa_z ) ) * ( ( SQR(z1) + SQR(1.0-z1) ) * ( 5.0/2.0 - SQR(M_PI)/3.0 + SQR( log(z1/(1.0-z1)) ) + OmegaT_V_unsymmetric(Q, z1, mf) + OmegaT_V_unsymmetric(Q, 1.0-z1, mf) + L_dip(Q,z1,mf)  ) 
        + (2.0*z1-1.0)/2.0 * (OmegaT_N_unsymmetric(Q, z1, mf)-OmegaT_N_unsymmetric(Q, 1.0-z1, mf) ) );
    double term2 = SQR( mf * gsl_sf_bessel_K0( x01 *kappa_z ) ) * ( 3.0 -SQR(M_PI)/3.0 + SQR(log(z1/(1.0-z1))) + OmegaT_V_unsymmetric(Q, z1, mf) + OmegaT_V_unsymmetric(Q, 1.0-z1, mf) + L_dip( Q, z1, mf )  );

    double res= term1 + term2;

    return res;
}


// Note Eq (32)
double OmegaT_V_unsymmetric( double Q, double z, double mf ) {
    // The first part of Omega^T_V(gamma; z) function that appears in the transverse NLOdip part.
    // Omega^T_V(gamma; z) = Omega^T_V_unsymmetric(gamma; z) + Omega^T_V_unsymmetric(gamma; 1-z)
    double gamma = sqrt( 1.0 + 4.0 * SQR(mf/Q) );
    double res = (1.0+1.0/(2.0*z)) * ( log(1.0-z) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) )
    - 1.0/(2.0*z) * ( (z + 0.5 ) * (1.0-gamma) + SQR(mf)/SQR(Q) ) * log( ( z*(1.0-z)*SQR(Q)+SQR(mf) )/SQR(mf) ) ;

    return res;
}

// Note eq (31)
double OmegaT_N_unsymmetric( double Q, double z, double mf ) {
    // The first part of Omega^T_V(gamma; z) function that appears in the transverse NLOdip part.
    // Omega^T_N(gamma; z) = Omega^T_N_unsymmetric(gamma; z) - Omega^T_N_unsymmetric(gamma; 1-z)
    double gamma = sqrt( 1.0 + 4.0 * SQR(mf/Q) );
    double res = (1.0+z-2.0*SQR(z))/z * ( log(1.0-z) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) ) 
    + (1.0-z)/z * ( (0.5+z)*(gamma-1.0) - SQR(mf)/SQR(Q) ) * log(  (z*(1.0-z)*SQR(Q)+SQR(mf))/SQR(mf)) ;

    return res;
}

// Note eq (42)
double ITdip_massive_1(double Q2, double z1, double x01sq, double mf, double xi)  {
    // One additional integral: xi from 0 to 1.

    double x01 = sqrt( x01sq );
    double Q = sqrt(Q2);

    double kappa_z = sqrt( z1*(1.0-z1)*Q2 + SQR(mf) );

    double term1 = kappa_z * gsl_sf_bessel_K1( x01 * kappa_z) * ( SQR(z1) + SQR(1.0-z1) ) * IT_V1(Q, z1, mf, x01, xi)  ;
    double term2 = SQR(mf) * gsl_sf_bessel_K0( x01 * kappa_z ) * IT_VMS1(Q, z1, mf, x01, xi) ;

    double res= term1 + term2;

    return res;

}

double IT_V1(double Q, double z, double mf, double r, double xi) {
    // I^T_V1 = IT_V1_unsymmetric(z) + IT_V1_unsymmetric(1-z)
    return IT_V1_unsymmetric(Q, z, mf, r, xi) + IT_V1_unsymmetric(Q, 1.0-z, mf, r, xi);
}
// Note (31)
double IT_V1_unsymmetric( double Q, double z, double mf, double r, double xi ) {
    // The first part of the unintegrated I^T_V1 function that appears in the transverse NLOdip part.
    // I^T_V1 = I^T_V1_unsymmetric(z) + I^T_V1_unsymmetric(1-z)
    // Note that this has to be integrated over xi from 0 to 1.

    double kappa_z = sqrt( z*(1.0-z)*SQR(Q) + SQR(mf) );

    double term1 = 1.0/xi * ( 2.0*log(xi)/(1.0-xi) - (1.0+xi)/2.0 ) * ( sqrt( SQR(kappa_z) + xi/(1.0-xi) * (1.0-z) * SQR(mf)) * gsl_sf_bessel_K1( r*sqrt(SQR(kappa_z) + xi/(1.0-xi) *(1.0-z)* SQR(mf)) ) - kappa_z * gsl_sf_bessel_K1( r*kappa_z ) );
    double term2 = -( log(xi)/SQR(1.0-xi) + z/(1.0-xi) + z/2.0 ) * (1.0-z)*SQR(mf)/sqrt( SQR(kappa_z) + xi/(1.0-xi) * (1.0-z) *SQR(mf) ) * gsl_sf_bessel_K1( r*sqrt( SQR(kappa_z) + xi/(1.0-xi) * (1.0-z) *SQR(mf) ) );

    double res = term1 + term2;

    return res;

}

double IT_VMS1(double Q, double z, double mf, double r, double xi)
{
    // I^T_VMS1 = IT_VMS1_unsymmetric(z) + IT_VMS1_unsymmetric(1-z)
    return IT_VMS1_unsymmetric(Q, z, mf, r, xi) + IT_VMS1_unsymmetric(Q, 1.0-z, mf, r, xi);
}
double IT_VMS1_unsymmetric( double Q, double z, double mf, double r, double xi ) {
    // The first part of the unintegrated I^T_VMS1 function that appears in the transverse NLOdip part.
    // I^T_VMS1 = I^T_VMS1_unsymmetric(z) + I^T_VMS1_unsymmetric(1-z)
    // Note that this has to be integrated over xi from 0 to 1.

    double kappa_z = sqrt( z*(1.0-z)*SQR(Q) + SQR(mf) );

    double term1 = 1.0/xi * (2.0 * log(xi)/(1.0-xi) - (1.0+xi)/(2.0)) * ( gsl_sf_bessel_K0( r*sqrt( SQR(kappa_z) + xi/(1.0-xi) * (1.0-z) *SQR(mf) ) ) - gsl_sf_bessel_K0(r*kappa_z) );
    double term2 = ( -3.0/2.0 * (1.0-z)/(1.0-xi) + (1.0-z)/2.0 ) * gsl_sf_bessel_K0( r* sqrt(SQR(kappa_z) + xi/(1.0-xi) * (1-z) *SQR(mf)) );

    double res = term1 + term2;

    return res;

}


//// Note (43)
double ITdip_massive_2(double Q2, double z1, double x01sq, double mf, double y_chi, double y_u) {
    // Two additional integrals: y_chi and y_u, both from 0 to 1


    double x01 = sqrt( x01sq );
    double Q = sqrt(Q2);

    double kappa_z = sqrt( z1*(1.0-z1)*SQR(Q) + SQR(mf) );

    double term1 = kappa_z * gsl_sf_bessel_K1( x01 * kappa_z) * ( 
        ( SQR(z1) + SQR(1.0-z1) ) * IT_V2(Q, z1, mf, x01, y_chi, y_u))   
        +  (2.0*z1-1.0)/2.0 * IT_N(Q, z1, mf, x01, y_chi, y_u) ;
    double term2 = SQR(mf) * gsl_sf_bessel_K0( x01 * kappa_z ) * IT_VMS2(Q, z1, mf, x01, y_chi, y_u) ;

    double res= term1 + term2;

    return res;

}

double IT_V2(double Q, double z, double mf, double r, double y_chi, double y_u) {
    // I^T_V2 = IT_V2_unsymmetric(z) + IT_V2_unsymmetric(1-z)
    return IT_V2_unsymmetric(Q, z, mf, r, y_chi, y_u) + IT_V2_unsymmetric(Q, 1.0-z, mf, r, y_chi, y_u);
}

double IT_V2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_V2 function that appears in the transverse NLOdip part.
    // I^T_V2 = I^T_V2_unsymmetric(z) + I^T_V2_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*SQR(Q) + SQR(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*SQR(Q) + SQR(mf) );

    double term1 = - 1.0/(1.0-chi) * 1.0/(u*(u+1.0)) * SQR(mf)/SQR(kappa_chi) * ( 2.0*chi + SQR(  u/(u+1.0)) * 1.0/z * (z-chi) * (1.0-2.0*chi) ) * ( sqrt(SQR(kappa_z) + u*(1.0-z)/(1.0-chi)*SQR(kappa_chi)) * gsl_sf_bessel_K1(r*sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi)*SQR(kappa_chi) )) - kappa_z *gsl_sf_bessel_K1(r* kappa_z) );
    double term2 = -1.0/SQR(1.0-chi) * 1.0/(u+1.0) * (z-chi) * ( 1.0 - 2.0*u/(1.0+u)*(z-chi) + SQR(u/(u+1.0)) *1.0/z * SQR(z-chi) ) * SQR(mf)/sqrt( SQR(kappa_z) + u *(1.0-z)/(1.0-chi) * SQR(kappa_chi)) * gsl_sf_bessel_K1( r* sqrt( SQR(kappa_z) + u *(1.0-z)/(1.0-chi) * SQR(kappa_chi) ) );

    double jacobian = z / SQR(y_u);

    double res = jacobian * (term1 + term2);

    return res;


}

double IT_VMS2(double Q, double z, double mf, double r, double y_chi, double y_u) {
    // I^T_VMS2 = IT_VMS2_unsymmetric(z) + IT_VMS2_unsymmetric(1-z)
    return IT_VMS2_unsymmetric(Q, z, mf, r, y_chi, y_u) + IT_VMS2_unsymmetric(Q, 1.0-z, mf, r, y_chi, y_u);
}

double IT_VMS2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_VMS2 function that appears in the transverse NLOdip part.
    // I^T_VMS2 = I^T_VMS2_unsymmetric(z) + I^T_VMS2_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*SQR(Q) + SQR(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*SQR(Q) + SQR(mf) );

    double term1 = 1.0/(1.0-chi) * 1.0/SQR(u+1.0) * (-z - u/(1.0+u) * (z + u*chi)/z * (chi-(1.0-z))) * gsl_sf_bessel_K0(r * sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi) *SQR(kappa_chi) ));
    double term2 = 1.0/(u+1.0)/SQR(u+1.0) * ( SQR(kappa_z)/SQR(kappa_chi) * (1.0 + u * chi*(1.0-chi) / ( z*(1.0-z) )) - SQR(mf)/SQR(kappa_chi) * chi/(1.0-chi) * (2.0 *SQR(1.0+u)/u + u/(z*(1.0-z)) *SQR(z-chi) ) ) * ( gsl_sf_bessel_K0( r *sqrt( SQR(kappa_z) + u* (1.0-z)/(1.0-chi) *SQR(kappa_chi)) ) - gsl_sf_bessel_K0(r* kappa_z) );

    double jacobian = z / SQR(y_u);

    double res = jacobian * (term1 + term2);

    return res;

}

double IT_N(double Q, double z, double mf, double r, double y_chi, double y_u) {
    // I^T_N = IT_N_unsymmetric(z) - IT_N_unsymmetric(1-z)
    return IT_N_unsymmetric(Q, z, mf, r, y_chi, y_u) - IT_N_unsymmetric(Q, 1.0-z, mf, r, y_chi, y_u);
}

double IT_N_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_N function that appears in the transverse NLOdip part.
    // I^T_N = I^T_N_unsymmetric(z) - I^T_N_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*SQR(Q) + SQR(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*SQR(Q) + SQR(mf) );

    double term1 = 2.0*(1.0-z)/z * 1.0/SQR(u+1.0)/(u+1.0) * ( (2.0+u)*u*z + SQR(u)*chi ) * sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi) * SQR(kappa_chi) ) * gsl_sf_bessel_K1( r* sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi) *SQR(kappa_chi) ) );
    double term2 = 2.0*(1.0-z)/z * 1.0/SQR(u+1.0)/(u+1.0) * SQR(mf)/SQR(kappa_chi) * ( z/(1.0-z) + chi/(1.0-chi) * (u-2.0*z-2.0*u*chi) ) * ( sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi)*SQR(kappa_chi)) * gsl_sf_bessel_K1( r * sqrt( SQR(kappa_z) + u*(1.0-z)/(1.0-chi) *SQR(kappa_chi) ) ) - kappa_z * gsl_sf_bessel_K1( r* kappa_z ) );

    double jacobian = z / SQR(y_u);

    double res = jacobian * (term1 + term2);

    return res;

}