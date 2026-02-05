#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_errno.h>

#include "nlodis.hpp"
#include "qcd.hpp"
#include "integration.hpp"

using namespace std;

static const std::string cubamethod = "suave";




/*
 * Structure function F2 
 *
 * Q2 [GeV^2]: photon virtuality
 * xbj: Bjorken-x
 * 
 */
double NLODIS::F2(double Q2, double xbj)
{
    double sigmaT = Photon_proton_cross_section(Q2, xbj, T);
    double sigmaL = Photon_proton_cross_section(Q2, xbj, L);

    return Q2 / (4.0 * M_PI * M_PI * ALPHA_EM) * (sigmaT + sigmaL);
}

/*
 * Structure function FL
 *
 * Q2 [GeV^2]: photon virtuality
 * xbj: Bjorken-x
 */
double NLODIS::FL(double Q2, double xbj)
{
    double sigmaL = Photon_proton_cross_section(Q2, xbj, L);
    return Q2 / (4.0 * M_PI * M_PI * ALPHA_EM) * sigmaL;
}


/*
 * qqg-target scattering amplitude
 * Ref. https://arxiv.org/pdf/2211.03504 (4)
 * In that notation, this is 1-S_{012}
 * 
 * x01, x02, x21: dipole sizes in GeV^-1
 * Y: evolution rapidity
 **/
double NLODIS::TripoleAmplitude(double x01, double x02, double x21, double Y) 
{
    double S01 = 1-dipole.DipoleAmplitude(x01, Y);
    double S02 = 1-dipole.DipoleAmplitude(x02, Y);
    double S12 = 1-dipole.DipoleAmplitude(x21, Y);

    if (nc_scheme == LargeNC)
    {
        return 1.0 - S02*S12;
    }
    else if (nc_scheme == FiniteNC)
    {
        return NC/(2.0*CF)*(S02*S12 - 1./SQR(NC)*S01);
    }
    else
    {
        throw std::runtime_error("NLODIS::TripoleAmplitude: unknown NC scheme");
    }
}

/*
 * Evolution rapidity (in the qqg contribution)
 *
 * Ref https://arxiv.org/pdf/2007.01645 eq (19)
 * xbj: Bjorken-x
 * Q2: photon virtuality in GeV^2
 * z2: gluon longitudinal momentum fraction
 */
double NLODIS::EvolutionRapidity(double xbj, double Q2, double z2) const
{
    double W2 = Q2 / xbj;
    return std::log(W2*z2/Q0sqr);
}

/*
 * Running coupling scale depending on the RC scheme used
 */
double NLODIS::RunningCouplinScale(double x01, double x02, double x21)
{
    if (rc_scheme == SMALLEST)
    {
        return std::min({x01, x02, x21});
    }
    else if (rc_scheme == PARENT)
    {
        return x01;
    }
    else
    {
        throw std::runtime_error("NLODIS::RunningCouplinScale: unknown running coupling scheme");
    }
}

/*
 * Photon-proton cross section [GeV^-2]
 *
 * To get the cross section, this has to be integrated over d^2r and multiplied by sigma_0 [in GeV^-2]
 * i.e. we replace 2\int d^2b -> sigma_0
 * 
 * */
double NLODIS::Photon_proton_cross_section(double Q2, double xbj, Polarization pol)
{
    if (scheme != UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }

    if (order==LO)
    {
        return Photon_proton_cross_section_LO(Q2, xbj, pol);
    }
    
    // NLO calculation

    double sigma_LO = Photon_proton_cross_section_LO(Q2, dipole.X0(), pol);
    double sigma_dip = Sigma_dip(Q2, xbj, pol);
    double  sigma_qg = Sigma_qg(Q2,xbj,pol);

    //cout <<"# Note: Sigma_LO: " << sigma_LO << " , Sigma_dip: " << sigma_dip << " , Sigma_qg: " << sigma_qg << endl;

    return sigma_LO + sigma_dip + sigma_qg;
}


/*
 * \sigma_dip
 * qq part of the NLO cross section
 * L polarization: https://arxiv.org/pdf/2103.14549 (166)
 * T polarization:
 */
double NLODIS::Sigma_dip(double Q2, double xbj, Polarization pol)
{
    double result=0;
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^2 below (from d^2x_{01} d^2b)
    double fac=4.0*NC*ALPHA_EM/SQR(2.0*M_PI); 
    

    for (const auto& quark : quarks) {
        IntegrationParams intparams;
        intparams.nlodis=this;
        intparams.Q2=Q2;
        intparams.xbj=xbj;
        intparams.pol=pol;

        intparams.quark=quark;
        if (pol == L)
        {
            // 1st line
            double I, Ierr, Iprob;
            intparams.contribution="Omega_L_const";
            Cuba(cubamethod, 2, integrand_dip_massive, &intparams, &I, &Ierr, &Iprob);
            result += I;

            // 2nd line of 2103.14549 (166)
            intparams.contribution="ab";
            double Iab,Iaberr,Iabprob;
            Cuba(cubamethod, 3, integrand_dip_massive, &intparams, &Iab, &Iaberr, &Iabprob);
            intparams.contribution="cd";

            double Icd,Icderr,Icdprob;
            Cuba(cubamethod, 4, integrand_dip_massive, &intparams, &Icd, &Icderr, &Icdprob);
            result += Iab + Icd;
        }
        else if (pol == T)
        {
            // T0 contribution
            intparams.contribution="T0";
            double IT0, IT0err, IT0prob;
            Cuba(cubamethod, 2, integrand_dip_massive, &intparams, &IT0, &IT0err, &IT0prob);
            result += IT0;

            // T1 contribution
            intparams.contribution="T1";
            double IT1, IT1err, IT1prob;
            Cuba(cubamethod, 3, integrand_dip_massive, &intparams, &IT1, &IT1err, &IT1prob);
            result += IT1;

            // T2 contribution
            intparams.contribution="T2";
            double IT2, IT2err, IT2prob;
            Cuba(cubamethod, 4, integrand_dip_massive, &intparams, &IT2, &IT2err, &IT2prob);
            result += IT2;
        }
        else
        {
            throw std::runtime_error("NLODIS::Sigma_dip: unknown polarization");
        }

        result *= SQR(quark.charge);

    } // Quark flavor loop

    // We have factorized out \int d^2 b - note the normalization convention!
    // Define sigma_0 = 2 \int d^2 b
    
    // Correspondingly I need to include 1/2
    result *= 1./2.;

        
    // 2pi from overall angular integral
    return fac*result*2.0*M_PI;
}


/*
 * Integrand wrapper for Cuba
 * Evaluate different contributions for the sigma_dip part
 * Integration variables are
 * x[0] = z1
 * x[1] = [0,1] mapped to r = x[1]*maxr
 * x[2] = xi (integration variable xi in (114))
 * x[3] = x (integration variable x in (114)) [when computing the "cd" contribution]
 * 
 * Note: 2pi from the overall angular integration of x01 in NLODIS::Sigma_dip
 */
int integrand_dip_massive(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    auto* const p = static_cast<IntegrationParams*>(userdata);

    if (p->pol == L)
    {
        if (!( (*ndim ==4 and p->contribution=="cd") or (*ndim == 3 and p->contribution =="ab") 
        or (*ndim ==2 and p->contribution=="Omega_L_const") ) )
        {
            cerr << "integrand_dip_massive: ndim " << *ndim << " and contribution " << p->contribution << " do not match" << endl;
            exit(1);
        }   
    }
    else
    {
        if (!( (*ndim ==4 and p->contribution=="T2") or (*ndim == 3 and p->contribution =="T1") 
        or (*ndim ==2 and p->contribution=="T0") ) )
        {
            cerr << "integrand_dip_massive: ndim " << *ndim << " and contribution " << p->contribution << " do not match" << endl;
            exit(1);
        }   
    }

   
    
    double Q2=p->Q2;
    double xbj=p->xbj;
    double mf=p->quark.mass;

    double z1=x[0];
    double x01=p->nlodis->GetMaxR()*x[1];
    double x01sq=SQR(x01);
    
    

    double alphabar=p->nlodis->Alphas(x01)*CF/M_PI;

    // TODO: add more user control for evolution rapidity
    double evolution_rapidity = std::log(1/xbj); 
    double dipole = p->nlodis->GetDipole().DipoleAmplitude(x01,evolution_rapidity);
    double res;
    /////////////// Longitudinal part ///////////////
    if (p->contribution=="ab" and p->pol == L) {
        // "ab" contribution does not have the x integration variable
        double xi=x[2]; 
        res = dipole*(ILdip_massive_Iab(Q2,z1,x01,mf,xi));
    } else if (p->contribution=="cd" and p->pol == L) {
        double xi=x[2]; 
        double intx=x[3];
        res = dipole*(ILdip_massive_Icd(Q2,z1,x01,mf,xi,intx));
    }
    else if (p->contribution=="Omega_L_const" and p->pol == L)
    {
        // only z and r integration
        res = dipole*ILdip_massive_Omega_L_Const(Q2, z1, x01, mf);
    }
    /////////////// Transverse part ///////////////
    else if (p->contribution=="T0" and p->pol == T)
    {
        res = dipole*ITdip_massive_0(Q2, z1, x01, mf);
    }
    else if (p->contribution=="T1" and p->pol == T)
    {
        double xi=x[2]; 
        res = dipole*ITdip_massive_1(Q2, z1, x01, mf, xi);
    }
    else if (p->contribution=="T2" and p->pol == T)
    {
        double xi=x[2]; 
        double intx=x[3];
        res = dipole*ITdip_massive_2(Q2, z1, x01, mf, xi, intx);
    }
    else 
    {
        cerr << "integrand_dip_massive: unknown contribution " << p->contribution << " pol " << p->pol << endl;
        exit(1);
    }

    double jacobian = x01 * p->nlodis->GetMaxR() * 2.0 * M_PI; // Jacobian from d^2r and r = u*MAXR
    res *= jacobian*alphabar; // Jacobian from d^2r

    if(std::isfinite(res)){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


/*
 * \sigma_qg
 * qg part of the NLO cross section
 * 
 * Longitudinal reference: (167) but instead of q^+, k^+ we integrate over z_i
 * Explicit expressoin is docs/NLO_DIS_cross_section_with_massive_quarks.pdf (13) 
*/
double NLODIS::Sigma_qg(double Q2, double xbj, Polarization pol)
{
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^3 below (from d^2x_{01} d^2x_{02} d^2b)
    double fac=4.0*NC*ALPHA_EM/std::pow(2.0*M_PI,3.0);
    

    double result=0;
    for (const auto& quark : quarks) {
        IntegrationParams intparams;
        intparams.nlodis=this;
        intparams.Q2=Q2;
        intparams.xbj=xbj;
        intparams.pol=pol;

        intparams.quark=quark;
        if (pol == L)
        {
            // Note (21): this contribution is split into 3 parts I_1, I_2 and I_3
            double I1,I1err,I1prob;
            intparams.contribution="I1";
            Cuba(cubamethod, 5, integrand_ILqgunsub_massive, &intparams, &I1, &I1err, &I1prob);
            result += I1;

            double I2,I2err,I2prob;
            intparams.contribution="I2";
            Cuba(cubamethod, 6, integrand_ILqgunsub_massive, &intparams, &I2, &I2err, &I2prob);
            result += I2;

            double I3,I3err,I3prob;
            intparams.contribution="I3";
            Cuba(cubamethod, 7, integrand_ILqgunsub_massive, &intparams, &I3, &I3err, &I3prob);
            result += I3;
        }
        result *= SQR(quark.charge);
    }

    


    // We have facotorized out (not performed) \int d^2 b - note the normalization convention!
    // Define sigma_0 = 2 \int d^2 b
    
    // Correspondingly I need to include 1/2 to this result
    result *= 1./2.;

    // 2pi: overall integral over one angle
    return  fac*2.0*M_PI*result;

}



 NLODIS::NLODIS(string bkdata) : dipole(bkdata)
 {
    // Initialize quark flavors and masses
    Quark u; u.type = Quark::U; u.mass = 0.14; u.charge = 2.0/3.0;
    Quark d; d.type = Quark::D; d.mass = 0.14; d.charge = -1.0/3.0;
    Quark s; s.type = Quark::S; s.mass = 0.14;  s.charge = -1.0/3.0;
    Quark c; c.type = Quark::C; c.mass = 1.27;   c.charge = 2.0/3.0;
    //Quark b; b.type = Quark::B; b.mass = 4.18;   b.charge = -1.0/3.0;

    quarks.push_back(u);
    quarks.push_back(d);
    quarks.push_back(s);
    quarks.push_back(c);
    //quarks.push_back(b);

    dipole.SetOutOfRangeErrors(false);

 }
 
 /*
  * Coordinate space coupling
  * TODO: implement flavor thresholds and C^2 scale
  */
 double NLODIS::Alphas(double r) const
 {
    const double LambdaQCD = 0.241; // GeV
    const double b0 = (11.0*NC - 2.0*quarks.size())/(12.0*M_PI);
    double mu2 = 4.0*C2_alpha/(r*r);
    double as = 1.0/(b0*log(mu2/(LambdaQCD*LambdaQCD)));
    if (as > 0.7 or log(mu2/(LambdaQCD*LambdaQCD))<0)
    {
        as = 0.7; // Freeze coupling
    }
    return as;
 }

 /*
  * Lower bound for the z2 integral
  * 
  * Ref https://arxiv.org/pdf/2007.01645 eq (18)
  * 
  */
double NLODIS::z2_lower_bound(double xbj, double Q2)
{
    double W2 = Q2 / xbj;
    return Q0sqr / W2;
}