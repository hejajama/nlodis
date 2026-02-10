/* 
 * NLO DIS cross section
 *
 * This file contains some internal functions needed to evaluate the NLO DIS
 * impact factors
 *  
 * Functions are defined such that they are suitable for Cuba 
*/

#include "nlodis.hpp"
#include <gsl/gsl_sf_dilog.h>
#include "qcd.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>


/*
* Longitudinal photon, dipole contribution with massive quarks
* This term corresponds to second line of 2103.14549 (166)
* i.e. 4 Q^2 z^2(1-z)^2 alpha_s*Nc/pi*K_0(r*eps)*I_{v, part (c)+(d)}
* I_{v, part (c)+(d)} is given in (114)
 */

double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf, double xi, double x);



/*
 * https://arxiv.org/pdf/2103.14549 (114) 
*/
double ILdip_massive_Icd(double Q2, double z1, double r, double mf, double xi, double x) {
    // Note that this now uses the earlier parametrization (xi,x) instead of (chi,u). This should be better for numerics.
    double front_factor = 4.0*Q2*SQR(z1*(1.0-z1));

    double kappa_z = sqrt( z1*(1.0-z1)*Q2 + SQR(mf) );
    double bessel_inner_fun = kappa_z * r;
    double Icd_integrand = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        Icd_integrand = 0;
    }else{
        // 2103.14549 (89,97)
        double CLm1 = SQR(z1)*(1.0-xi)/(1.0-z1) * ( -SQR(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+z1*xi/(1.0-z1)) ) / (x*(1.0-xi)+xi/(1.0-z1)) ); // C^L_m(z,x,xi)
        double CLm2 = SQR(1.0-z1)*(1.0-xi)/(z1) * ( -SQR(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+(1.0-z1)*xi/(z1)) ) / (x*(1.0-xi)+xi/(z1)) );; // C^L_m(1-z,x,xi)
        double kappa1 = xi*SQR(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) ) * ( xi*(1.0-x) + x*(1.0-z1*(1.0-xi)/(1.0-z1)) ); // kappa(z,x,xi)
        double kappa2 = xi*SQR(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) ) * ( xi*(1.0-x) + x*(1.0-(1.0-z1)*(1.0-xi)/(z1)) );; // kappa(1-z,x,xi)
        Icd_integrand = gsl_sf_bessel_K0( bessel_inner_fun ) * SQR(mf)*
                     ( (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( r*sqrt( SQR(kappa_z)/(1.0-x) + kappa1 ) ) ) *
                        CLm1 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) * ( x/(1.0-x)*SQR(kappa_z) + kappa1 ) ) +
                        (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( r*sqrt( SQR(kappa_z)/(1.0-x) + kappa2 ) ) ) *
                        CLm2 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) * ( x/(1.0-x)*SQR(kappa_z) + kappa2 ) ));
    }   

    double dip_res = front_factor * Icd_integrand;
    return dip_res;
}

/*
 * https://arxiv.org/pdf/2103.14549 (113) 
*/
double ILdip_massive_Iab(double Q2, double z1, double r, double mf, double xi) {

    double front_factor = 4.0*Q2*SQR(z1*(1.0-z1));

    double kappa_z = sqrt( z1*(1.0-z1)*Q2 + SQR(mf) );
    double bessel_inner_fun = kappa_z * r;
    double bessel_arg_2=sqrt( SQR(kappa_z) + (1.0-z1)*xi/(1.0-xi) *SQR(mf) ) * r;
    double bessel_arg_3= sqrt( SQR(kappa_z) + z1*xi/(1.0-xi) *SQR(mf) ) * r;

    double b1 = gsl_sf_bessel_K0( bessel_inner_fun);
    double b2 = gsl_sf_bessel_K0( bessel_arg_2 );
    double b3 = gsl_sf_bessel_K0( bessel_arg_3 );


    double Iab_integrand = b1 * 1.0/xi * ( -2.0*log(xi)/(1.0-xi) + (1.0+xi)/2.0 ) *
                        (2.0*b1 - b2 - b3); 
      

    double dip_res = front_factor * Iab_integrand;
    return dip_res;
}

/*
 * Integrand on the 2nd line of (166) in https://arxiv.org/pdf/2103.14549
 */
double OmegaL_V( double Q2, double z, double mf );
double L_dip( double Q2, double z, double mf );
double ILdip_massive_Omega_L_Const(double Q2, double z1, double r, double mf) 
{
    double front_factor = 4.0*Q2*SQR(z1*(1.0-z1));
    double bessel_inner_fun = sqrt( Q2*z1*(1.0-z1) + SQR(mf))*r;
    double dip_res = 0;
    if (bessel_inner_fun < 1e-7){
        dip_res = 0;
    }else{
        dip_res = front_factor * SQR(gsl_sf_bessel_K0( bessel_inner_fun )) * 
        ( 5.0/2.0 - SQR(M_PI)/3.0 + SQR(log( z1/(1.0-z1) )) + OmegaL_V(Q2,z1,mf) + L_dip(Q2,z1,mf) );
    }   

    return dip_res;
}

// (100)
double OmegaL_V( double Q2, double z, double mf ) {
    // The Omega^L_V(gamma; z) function that appears in the longitudinal NLOdip part.
    double gamma = sqrt( 1.0 + 4.0 * SQR(mf)/Q2);
    double res = 1.0/(2.0*z) * ( log( 1.0-z ) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) ) 
                +1.0/(2.0*(1.0-z)) * ( log( z ) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*(1.0-z)) ))
                +1.0/( 4.0*z*(1.0-z) ) * ( gamma-1.0 + 2.0 * SQR(mf)/Q2 ) * log( (z*(1.0-z)*Q2 + SQR(mf) ) / SQR(mf) )  ;

    return res;
}
//(101)
double L_dip( double Q2, double z, double mf ) {
    // The L(gamma; z) function that appears in the transverse and longitudinal NLOdip part. 
    double gamma = sqrt( 1.0 + 4.0 * SQR(mf)/Q2);

    double res = gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * z) * ( 1.0 - gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * z) * ( 1.0 + gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * (1.0-z)) * ( 1.0 - gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * (1.0-z)) * ( 1.0 + gamma) ) );

    return res;
}


/***********************************************************
 * Longitudinal photon, qqg contribution
 */

 



/*
 * (22) in the note docs/NLO_DIS_cross_section_with_massive_quarks.pdf
 * 
 */


double ILNLOqg_massive_tripole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // sigma_qg divided into three parts. This part has no additional integrals. Only contains the part proportional to N_012

    double front_factor = 4.0*Q2;
    double Q = std::sqrt(Q2);

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*SQR(z0+z2));
    double omega_l = z1*z2/(z0*SQR(z1+z2));

    double x3_k = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );
    
    double term_k = 1.0/x02sq * SQR(z1) * ( 2.0*z0*(z0+z2) +SQR(z2) ) 
        * SQR( gsl_sf_bessel_K0( sqrt( SQR(Qbar_k) + SQR(mf) ) 
        * sqrt( SQR(x3_k) + omega_k * x02sq ) )  );
    double term_l = 1.0/x21sq * SQR(z0) * ( 2.0*z1*(z1+z2) +SQR(z2) ) 
        * SQR( gsl_sf_bessel_K0( sqrt( SQR(Qbar_l) + SQR(mf) ) 
        * sqrt( SQR(x3_l) + omega_l * x21sq ) )  );
    double term_kl = -2.0 * z0 *z1 * ( z0*(1.0-z0) + z1*(1.0-z1) ) * x20x21 / (x02sq * x21sq) 
        * gsl_sf_bessel_K0( sqrt( SQR(Qbar_k) + SQR(mf)) * sqrt(SQR(x3_k) + omega_k * x02sq) ) 
        * gsl_sf_bessel_K0( sqrt( SQR(Qbar_l) + SQR(mf)) * sqrt(SQR(x3_l) + omega_l * x21sq) );


    double res = front_factor * ( term_k + term_l + term_kl );
    return res;
}



double ILNLOqg_massive_dipole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // sigma_qg divided into three parts. This part has no additional integrals. Only contains the part proportional to N_01

    double front_factor = 4.0*Q2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Q = std::sqrt(Q2);
    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*SQR(z0+z2));
    double omega_l = z1*z2/(z0*SQR(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x3_k = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );



    double term_k = -1.0/x02sq * SQR(z1) * ( 2.0*z0*(z0+z2) +SQR(z2) ) * exp( -x02sq / x01sq / exp(M_EULER) ) * SQR( gsl_sf_bessel_K0( sqrt( SQR(Qbar_k) + SQR(mf) ) * sqrt(x01sq) )  );
    double term_l = -1.0/x21sq * SQR(z0) * ( 2.0*z1*(z1+z2) +SQR(z2) ) * exp( -x21sq / x01sq / exp(M_EULER) ) * SQR( gsl_sf_bessel_K0( sqrt( SQR(Qbar_l) + SQR(mf) ) * sqrt(x01sq) )  );

    double res = front_factor * ( term_k + term_l );
    return res;
}

/*
 * (23) in the note docs/NLO_DIS_cross_section_with_massive_quarks.pdf
 * But with one integral performed analytically, only y_t integral remaining
*/ 
double ILNLOqg_massive_tripole_part_I2_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t) {
    // sigma_qg divided into three parts. This part has one additional integral.

    double front_factor = 4.0*Q2;

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double Q = std::sqrt(Q2);
    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*SQR(z0+z2));
    double omega_l = z1*z2/(z0*SQR(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;
    double x2_k = sqrt( x02sq );
    double x2_l = sqrt( x21sq );

    double x3_k = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double int_12_bar_k = G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t)-G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t);
    double int_12_bar_l = G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l, y_t)-G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, 0.0, y_t);

    double term_k = SQR(z1) * ( 2.0*z0*(z0+z2) +SQR(z2) ) * 1.0/4.0 * int_12_bar_k
                    * gsl_sf_bessel_K0( sqrt( SQR(Qbar_k) + SQR(mf) ) *sqrt( SQR(x3_k) + omega_k * x02sq )   );
    double term_l = SQR(z0) * ( 2.0*z1*(z1+z2) +SQR(z2) ) * 1.0/4.0 * int_12_bar_l 
                    * gsl_sf_bessel_K0( sqrt( SQR(Qbar_l) + SQR(mf) ) *sqrt( SQR(x3_l) + omega_l * x21sq )   );
    double term_kl = -1.0/4.0 * z0*z1 * ( z0*(1.0-z0) + z1*(1.0-z1) ) * x20x21  * (
        1.0/( x21sq ) * int_12_bar_k * gsl_sf_bessel_K0( sqrt( SQR(Qbar_l) + SQR(mf) ) *sqrt( SQR(x3_l) + omega_l * x21sq ))
        + 1.0/( x02sq ) * int_12_bar_l * gsl_sf_bessel_K0( sqrt( SQR(Qbar_k) + SQR(mf) ) *sqrt( SQR(x3_k) + omega_k * x02sq ))
    );

    double res = front_factor * ( term_k + term_l + term_kl );
    return res;
}

/*
 * I3 for qqg contribution
 * Note 24 but some integrals performed analytically 
 * */

double ILNLOqg_massive_tripole_part_I3_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2) {
    // sigma_qg divided into three parts. This part has two additional integrals.

    double front_factor = 4.0*Q2;

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double Q = std::sqrt(Q2);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*SQR(z0+z2));
    double omega_l = z1*z2/(z0*SQR(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;
    double x2_k = sqrt( x02sq );
    double x2_l = sqrt( x21sq );

    double x3_k = sqrt( SQR(z0) / SQR(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( SQR(z1) / SQR(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double int_12_bar_k1 = G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t1)-G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t1);
    double int_12_bar_l1 = G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l, y_t1)-G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, 0.0, y_t1);
    double int_12_bar_k2 = G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t2)-G_integrand_simplified( 1, 2, Qbar_k, mf, x2_k, x3_k, omega_k, 0.0, y_t2);
    double int_12_bar_l2 = G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l, y_t2)-G_integrand_simplified( 1, 2, Qbar_l, mf, x2_l, x3_l, omega_l, 0.0, y_t2);


    double int_11_k1 = G_integrand_simplified( 1, 1, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t1);
    double int_11_l1 = G_integrand_simplified( 1, 1, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l, y_t1);
    double int_11_k2 = G_integrand_simplified( 1, 1, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k, y_t2);
    double int_11_l2 = G_integrand_simplified( 1, 1, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l, y_t2);




    double term_k = SQR(z1) * ( 2.0*z0*(z0+z2) +SQR(z2) ) * x02sq/64.0 * int_12_bar_k1 * int_12_bar_k2 ;
    double term_l = SQR(z0) * ( 2.0*z1*(z1+z2) +SQR(z2) ) * x21sq/64.0 * int_12_bar_l1 * int_12_bar_l2;
    double term_kl = -1.0/32.0 * z1*z0 * ( z1*(1.0-z1) + z0*(1.0-z0) ) * x20x21 * int_12_bar_k1 * int_12_bar_l2;
    double term_mf = SQR(mf)/16.0 * SQR(z2) * SQR(z2) * (  
        SQR(z1/(z0+z2)) * int_11_k1 * int_11_k2 + SQR( z0/(z1+z2) ) * int_11_l1 * int_11_l2 - 2.0 * z0/(z1+z2) * z1/(z0+z2) * int_11_k1 * int_11_l2
    );


    double res = front_factor * ( term_k + term_l + term_kl + term_mf );
    return res;
}


/*
 * G_(x)^{a,b) with u integrated analytically
 * See doc/NLO_DIS_cross_section_with_massive_quarks.pdf (69)
 * Unintegrated form is (163) in 2103.14549
*/
double G_integrand_simplified(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda, double y) {
    // Integrands of the functions G^(a;b)_(x) that appear in the qqg-part. Here the u-integral has been done, and only one integral remains.
    // The integration variable is y = omega * t

    double value = 1.0/pow( y, 0.5*(2.0-a+b) ) * pow( 2.0, a+b-1.0 ) * pow(omega, b-1.0) 
                    * pow(( y*lambda * SQR(mf) + SQR(Qbar) + SQR(mf) ) / ( y*SQR(x3) + omega*SQR(x2) ), 0.5*(a+b-2.0) )
                    * gsl_sf_bessel_Kn( a+b-2, sqrt(1.0/y * ( y*lambda * SQR(mf) + SQR(Qbar) + SQR(mf) ) * ( y*SQR(x3) + omega*SQR(x2) )) ) ;

    return value;
}

