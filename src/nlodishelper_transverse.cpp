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


////////// Transverse qqg
//////////

// I1 = Note (54)

double ITNLOqg_massive_dipole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double Q = sqrt(Q2);
    return IT_dipole_jk_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_dipole_jkm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq );
}

double ITNLOqg_massive_tripole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){
    double Q = sqrt(Q2);
    return IT_tripole_jk_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_tripole_jkm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_tripole_F_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) +
           IT_tripole_Fm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq );
}


////////////////// Transverse I2

double ITNLOqg_massive_tripole_part_I2_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t){
    double Q = sqrt(Q2);
    return IT_tripole_jk_I2_fast(  Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t ) + 
           IT_tripole_jkm_I2_fast( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t  ) + 
           IT_tripole_F_I2_fast(   Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t  ) +
           IT_tripole_Fm_I2_fast(  Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t  );
}





// I_1, (54) of the note docs/NLO_DIS_cross_section_with_massive_quarks.pdf
// First part proportional to N_01
double IT_dipole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) * ( SQR(Qbar_j) + SQR(mf) ) / SQR(x2_j) 
                    * ( - exp( -SQR(x2_j) / ( x01sq *exp(M_EULER) ) ) * SQR(gsl_sf_bessel_K1( sqrt( x01sq * ( SQR(Qbar_j) + SQR(mf) ) ) ) ));
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) * ( SQR(Qbar_k) + SQR(mf) ) / SQR(x2_k) 
                    * ( - exp( -SQR(x2_k) / ( x01sq *exp(M_EULER) ) ) * SQR(gsl_sf_bessel_K1( sqrt( x01sq * ( SQR(Qbar_k) + SQR(mf) ) ) ) ));
   
    double res = term_j + term_k;

    return res;
}

// Part of (54) proportional to N_012
double IT_tripole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) * ( SQR(Qbar_j) + SQR(mf) ) / SQR(x2_j) 
                    * SQR(x3_j) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) * SQR( gsl_sf_bessel_K1( sqrt( SQR(x3_j) + omega_j * SQR(x2_j) ) * sqrt( SQR(Qbar_j) + SQR(mf) ) ) ) ;
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) * ( SQR(Qbar_k) + SQR(mf) ) / SQR(x2_k) 
                    * SQR(x3_k) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) * SQR( gsl_sf_bessel_K1( sqrt( SQR(x3_k) + omega_k * SQR(x2_k) ) * sqrt( SQR(Qbar_k) + SQR(mf) ) ) ) ;
   
    double res = term_j + term_k;

    return res;
}

double IT_dipole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2)) / SQR(x2_j) 
                    * ( - exp( -SQR(x2_j) / ( x01sq *exp(M_EULER) ) ) * SQR(gsl_sf_bessel_K0( sqrt( x01sq * ( SQR(Qbar_j) + SQR(mf) ) ) ) ));
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2)) / SQR(x2_k) 
                    * ( - exp( -SQR(x2_k) / ( x01sq *exp(M_EULER) ) ) * SQR(gsl_sf_bessel_K0( sqrt( x01sq * ( SQR(Qbar_k) + SQR(mf) ) ) ) ));
   
    double res = SQR(mf) * (term_j + term_k);

    return res;
}

double IT_tripole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2)) / SQR(x2_j) 
                    * SQR( gsl_sf_bessel_K0( sqrt( SQR(x3_j) + omega_j * SQR(x2_j) ) * sqrt( SQR(Qbar_j) + SQR(mf) ) ) );
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2)) / SQR(x2_k) 
                    * SQR( gsl_sf_bessel_K0( sqrt( SQR(x3_k) + omega_k * SQR(x2_k) ) * sqrt( SQR(Qbar_k) + SQR(mf) ) ) );
   
    double res = SQR(mf) * (term_j + term_k);

    return res;
}



double IT_tripole_F_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;

    double G22_sing_j = 1.0/SQR(x2_j) * sqrt( ( SQR(Qbar_j) + SQR(mf) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_j) + SQR(mf) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double G22_sing_k = 1.0/SQR(x2_k) * sqrt( ( SQR(Qbar_k) + SQR(mf) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_k) + SQR(mf) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double term_1 = 4.0/(z0+z2)/(z1+z2) * ( z2*SQR(z0-z1) * ( x2j_x3j * x2k_x3k - x2k_x3j * x2j_x3k ) 
                    - ( z1*(z0+z2) + z0*(z1+z2) )*(z0*(z0+z2)+z1*(z1+z2)) * x20x21 * x3j_x3k ) *
                    G22_sing_j * G22_sing_k;
    double term_2j = -(z0+z2)*z1*z2/SQR(z1+z2) * x2j_x3j * H_k * G22_sing_j;
    double term_2k = (z1+z2)*z0*z2/SQR(z0+z2) * x2k_x3k * H_j * G22_sing_k;
    double term_3j = -SQR(z0)*z1*z2/(z0+z2)/SQR(z0+z2) * x2j_x3j * H_j * G22_sing_j ;
    double term_3k = SQR(z1)*z0*z2/(z1+z2)/SQR(z1+z2) * x2k_x3k * H_k * G22_sing_k ;
    double term_4j = SQR(z0*z2)/(8.0*SQR(z0+z2)*SQR(z0+z2))*SQR(H_j);
    double term_4k = SQR(z1*z2)/(8.0*SQR(z1+z2)*SQR(z1+z2))*SQR(H_k);

    double res = 0.5 * (term_1 + term_2j + term_2k + term_3j + term_3k + term_4j + term_4k);

    return res;

}


double IT_tripole_Fm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );


    double G12_sing_j = 1.0/SQR(x2_j) * gsl_sf_bessel_K0( sqrt(( SQR(Qbar_j) + SQR(mf) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double G12_sing_k = 1.0/SQR(x2_k) * gsl_sf_bessel_K0( sqrt(( SQR(Qbar_k) + SQR(mf) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );


    double res = 0.5*SQR(mf) * (
        -1.0/(32.0*(z0+z2)*(z1+z2)) * ( (2.0*z0+z2)*(2.0*z1+z2)+SQR(z2) ) *x20x21 * 8.0 * G12_sing_j * 8.0 *G12_sing_k
    );

    return res;
}




////// Terms for I2
double IT_tripole_F_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;



    double int_22_bar_j = G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t) - G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, 0.0, y_t);
    double int_22_bar_k = G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t) - G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);

    double G22_sing_j = 1.0/SQR(x2_j) * sqrt( ( SQR(Qbar_j) + SQR(mf) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_j) + SQR(mf) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double G22_sing_k = 1.0/SQR(x2_k) * sqrt( ( SQR(Qbar_k) + SQR(mf) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_k) + SQR(mf) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double term_1 = 1.0/(4.0*(z0+z2)*(z1+z2)) * (z2*SQR(z0-z1)* ( x2j_x3j * x2k_x3k - x2k_x3j * x2j_x3k ) 
                    - (z1*(z0+z2)+z0*(z1+z2))* (z0*(z0+z2)+z1*(z1+z2))*x20x21 * x3j_x3k  )
                    * (
                        int_22_bar_k * G22_sing_j
                        +int_22_bar_j * G22_sing_k
                    );
    double term_2j = -(z0+z2)*z1*z2/(16.0*SQR(z1+z2)) * x2j_x3j * H_k *int_22_bar_j;
    double term_2k = (z1+z2)*z0*z2/(16.0*SQR(z0+z2)) * x2k_x3k * H_j *int_22_bar_k;
    double term_3j = -SQR(z0)*z1*z2/(16.0*(z0+z2)*SQR(z0+z2)) * x2j_x3j*H_j*int_22_bar_j;
    double term_3k = SQR(z1)*z0*z2/(16.0*(z1+z2)*SQR(z1+z2)) * x2k_x3k*H_k*int_22_bar_k;

    double res = 0.5 * (term_1 + term_2j + term_2k + term_3j + term_3k);

    return res;
}

double IT_tripole_Fm_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;


    double int_12_bar_j = G_integrand_simplified( 1, 2, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t) - G_integrand_simplified( 1, 2, Qbar_j, mf, x2_j, x3_j, omega_j, 0.0, y_t);
    double int_12_bar_k = G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t) - G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);

    double int_22_bar_j = G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t) - G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, 0.0, y_t);
    double int_22_bar_k = G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t) - G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);


    double int_21_j = G_integrand_simplified( 2, 1, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t);
    double int_21_k = G_integrand_simplified( 2, 1, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t);

    double int_11_j = G_integrand_simplified( 1, 1, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t);
    double int_11_k = G_integrand_simplified( 1, 1, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t);

    double G12_sing_j = 1.0/SQR(x2_j) * gsl_sf_bessel_K0( sqrt(( SQR(Qbar_j) + SQR(mf) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double G12_sing_k = 1.0/SQR(x2_k) * gsl_sf_bessel_K0( sqrt(( SQR(Qbar_k) + SQR(mf) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double G22_sing_j = 1.0/SQR(x2_j) * sqrt( ( SQR(Qbar_j) + SQR(mf) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_j) + SQR(mf) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double G22_sing_k = 1.0/SQR(x2_k) * sqrt( ( SQR(Qbar_k) + SQR(mf) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( SQR(Qbar_k) + SQR(mf) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) / ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_j) + SQR(mf) * (1.0+lambda_j) ) * ( SQR(x3_j) + omega_j * SQR(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) / ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( SQR(Qbar_k) + SQR(mf) * (1.0+lambda_k) ) * ( SQR(x3_k) + omega_k * SQR(x2_k) ) ) );

    double term_1j = -z0*z1*SQR(z2)/(16.0*(z0+z2)*SQR(z0+z2)) * x2j_x3j * int_21_j * 8.0 * G12_sing_j;
    double term_1k =  z0*z1*SQR(z2)/(16.0*(z1+z2)*SQR(z1+z2)) * x2k_x3k * int_21_k * 8.0 * G12_sing_k;
    double term_2 = -1.0/(32.0*(z0+z2)*(z1+z2))* ( (2.0*z0+z2)*(2.0*z1+z2) + SQR(z2) ) * x20x21 * (
         int_12_bar_k * 8.0 * G12_sing_j
        +int_12_bar_j * 8.0 * G12_sing_k
    );
    double term_3j = -SQR(z0*z2)/(16.0*(z0+z2)*SQR(z1+z2)) * x2j_x3k * int_21_k * 8.0 * G12_sing_j;
    double term_3k =  SQR(z1*z2)/(16.0*SQR(z0+z2)*(z1+z2)) * x2k_x3j * int_21_j * 8.0 * G12_sing_k;
    double term_4j = -z0*z1*SQR(z2)/(16.0*(z0+z2)*SQR(z0+z2)) * x2j_x3j * int_11_j * 16.0 * G22_sing_j ;
    double term_4k =  z0*z1*SQR(z2)/(16.0*(z1+z2)*SQR(z1+z2)) * x2k_x3k * int_11_k * 16.0 * G22_sing_k ;
    double term_5j = -(z0+z2)*SQR(z2)/(16.0*SQR(z1+z2)) * x2j_x3j * int_11_k * 16.0 * G22_sing_j ;
    double term_5k =  (z1+z2)*SQR(z2)/(16.0*SQR(z0+z2)) * x2k_x3k * int_11_j * 16.0 * G22_sing_k ;
    double term_6j = z0*z2*SQR(z2)/(4.0*SQR(z0+z2)*SQR(z0+z2))* H_j * int_11_j;
    double term_6k = z1*z2*SQR(z2)/(4.0*SQR(z1+z2)*SQR(z1+z2))* H_k * int_11_k;

    double res = 0.5*SQR(mf) * ( term_1j + term_1k + term_2 + term_3j + term_3k + term_4j + term_4k + term_5j + term_5k + term_6j + term_6k );

    return res;
}


double IT_tripole_jk_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );



    double int_22_bar_j = G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t) - G_integrand_simplified( 2, 2, Qbar_j, mf, x2_j, x3_j, omega_j, 0.0, y_t);
    double int_22_bar_k = G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t) - G_integrand_simplified( 2, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);


    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) 
                    * int_22_bar_j * SQR(x3_j)/8.0 * sqrt( ( SQR(Qbar_j) + SQR(mf) ) / (SQR(x3_j) + omega_j * SQR(x2_j) ) )
                    * gsl_sf_bessel_K1( sqrt( SQR(x3_j) + omega_j * SQR(x2_j) ) * sqrt( SQR(Qbar_j) + SQR(mf) ) )  ;
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) 
                    * int_22_bar_k * SQR(x3_k)/8.0 * sqrt( ( SQR(Qbar_k) + SQR(mf) ) / (SQR(x3_k) + omega_k * SQR(x2_k) ) )
                    * gsl_sf_bessel_K1( sqrt( SQR(x3_k) + omega_k * SQR(x2_k) ) * sqrt( SQR(Qbar_k) + SQR(mf) ) )  ;

    double res = term_j + term_k;

    return res;
}

double IT_tripole_jkm_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*SQR(z0+z2));
    double omega_k = z1*z2/(z0*SQR(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );


    double int_12_bar_j = G_integrand_simplified( 1, 2, Qbar_j, mf, x2_j, x3_j, omega_j, lambda_j, y_t) - G_integrand_simplified( 1, 2, Qbar_j, mf, x2_j, x3_j, omega_j, 0.0, y_t);
    double int_12_bar_k = G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t) - G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);

    double term_j = 1.0/SQR(z0+z2) * (2.0*z0*(z0+z2)+SQR(z2))
                    * int_12_bar_j/4.0 * gsl_sf_bessel_K0( sqrt( SQR(x3_j) + omega_j * SQR(x2_j) ) * sqrt( SQR(Qbar_j) + SQR(mf) ) )  ;
    double term_k = 1.0/SQR(z1+z2) * (2.0*z1*(z1+z2)+SQR(z2))
                    * int_12_bar_k/4.0 * gsl_sf_bessel_K0( sqrt( SQR(x3_k) + omega_k * SQR(x2_k) ) * sqrt( SQR(Qbar_k) + SQR(mf) ) )  ;

    double res = SQR(mf) * (term_j + term_k);

    return res;
}
