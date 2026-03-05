
#include <gsl/gsl_errno.h>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cuba.h>

#include "nlodis.hpp"
#include "qcd.hpp"
#include "integration.hpp"
#include "datatypes.hpp"

using namespace std;

// Do not allow z_i < zmin, to avoid numerical instabilities. This is not a physical cutoff, but just a technical one for the numerical integration.
const double ZMIN=1e-8;

double NLODIS::F2(double Q2, double xbj)
{
    double sigmaT = Photon_proton_cross_section_d2b(Q2, xbj, Polarization::T);
    double sigmaL = Photon_proton_cross_section_d2b(Q2, xbj, Polarization::L);

    double res = Q2 / (4.0 * M_PI * M_PI * Constants::AlphaEM) * (sigmaT + sigmaL);

    res *= ProtonTransverseArea(); // Include \int d^2 b
    return res;
}

double NLODIS::FL(double Q2, double xbj)
{
    double sigmaL = Photon_proton_cross_section_d2b(Q2, xbj, Polarization::L);
    double res = Q2 / (4.0 * M_PI * M_PI * Constants::AlphaEM) * sigmaL;
    res *= ProtonTransverseArea(); // Include \int d^2 b
    return res;
}

double NLODIS::FT(double Q2, double xbj)
{
    double sigmaT = Photon_proton_cross_section_d2b(Q2, xbj, Polarization::T);
    double res = Q2 / (4.0 * M_PI * M_PI * Constants::AlphaEM) * sigmaT;
    res *= ProtonTransverseArea(); // Include \int d^2 b
    return res;
}   

double NLODIS::TripoleAmplitude(double x01, double x02, double x21, double Y) 
{
    double S01 = 1-dipole->DipoleAmplitude(x01, Y);
    double S02 = 1-dipole->DipoleAmplitude(x02, Y);
    double S12 = 1-dipole->DipoleAmplitude(x21, Y);

    if (config.nc_scheme == NcScheme::LargeNC)
    {
        cerr << "Warning: Using large-Nc scheme for qqg-target amplitude, but this does not set CF to NC/2 in other parts of the code, i.e. is not a fully consistent large-nc limit!" << endl;
        return 1.0 - S02*S12;
    }
    else if (config.nc_scheme == NcScheme::FiniteNC)
    {
        return 1.0- Constants::NC/(2.0*Constants::CF)*(S02*S12 - 1./SQR(Constants::NC)*S01);
    }
    else
    {
        throw std::runtime_error("NLODIS::TripoleAmplitude: unknown NC scheme");
    }
}

double NLODIS::EvolutionRapidity(double xbj, double Q2, double z2) const
{
    double W2 = Q2 / xbj;
    return std::log(W2*z2/NLODISConfig::Q0sqr);
}

double NLODIS::RunningCouplinScale(double x01, double x02, double x21) const
{
    if (config.rc_scheme == RunningCouplingScheme::SMALLEST)
    {
        return std::min({x01, x02, x21});
    }
    else if (config.rc_scheme == RunningCouplingScheme::PARENT)
    {
        return x01;
    }
    else
    {
        throw std::runtime_error("NLODIS::RunningCouplinScale: unknown running coupling scheme");
    }
}

double NLODIS::Photon_proton_cross_section_d2b(double Q2, double xbj, Polarization pol)
{
    if (config.scheme != SubtractionScheme::UNSUB)
    {
        throw std::runtime_error("Only UNSUB scheme is implemented.");
    }

    if (config.order==Order::LO)
    {
        return Photon_proton_cross_section_LO_d2b(Q2, xbj, pol);
    }
    
    // NLO calculation

    double sigma_LO = Photon_proton_cross_section_LO_d2b(Q2, dipole->X0(), pol);
    double sigma_dip = Sigma_dip_d2b(Q2, xbj, pol);
    double  sigma_qg = Sigma_qg_d2b(Q2,xbj,pol);

    return sigma_LO + sigma_dip + sigma_qg; 
}

/*
 * NLO-dip term
 * References:
 * L polarization: https://arxiv.org/pdf/2103.14549 (166)
 * T polarization:
 */

double NLODIS::Sigma_dip_d2b(double Q2, double xbj, Polarization pol)
{
    double result=0;
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^2 below (from d^2x_{01} d^2b)
    double fac=4.0*Constants::NC*Constants::AlphaEM/SQR(2.0*M_PI); 

    if (config.sigma_dip_scheme != SigmaDipScheme::AnalyticalZ2Int)
    {
        throw std::runtime_error("Only AnalyticalZ2Int scheme for sigma_dip is implemented in Sigma_dip_d2b.");
    }
    

    for (const auto& quark : quarks) {
        double tmpresult=0; // This quark flavor contribution 
    
        IntegrationParams intparams;
        intparams.nlodis=this;
        intparams.Q2=Q2;
        intparams.xbj=xbj;
        intparams.pol=pol;

        intparams.quark=quark;
        if (pol == Polarization::L)
        {
            // 1st line
            double I, Ierr, Iprob;
            intparams.contribution="Omega_L_const";
            Cuba("cuhre", 2, integrand_dip_massive, &intparams, &I, &Ierr, &Iprob, cuba_config);
            tmpresult += I;

            // 2nd line of 2103.14549 (166)
            intparams.contribution="ab";
            double Iab,Iaberr,Iabprob;
            Cuba(cuba_config.method, 3, integrand_dip_massive, &intparams, &Iab, &Iaberr, &Iabprob, cuba_config);
            
            intparams.contribution="cd";
            double Icd,Icderr,Icdprob;
            Cuba(cuba_config.method, 4, integrand_dip_massive, &intparams, &Icd, &Icderr, &Icdprob, cuba_config);
            tmpresult += Iab + Icd;

        }
        else if (pol == Polarization::T)
        {
            // T0 contribution
            intparams.contribution="T0";
            double IT0, IT0err, IT0prob;
            Cuba("cuhre", 2, integrand_dip_massive, &intparams, &IT0, &IT0err, &IT0prob, cuba_config);
            tmpresult += IT0;

            // T1 contribution
            intparams.contribution="T1";
            double IT1, IT1err, IT1prob;
            Cuba(cuba_config.method, 3, integrand_dip_massive, &intparams, &IT1, &IT1err, &IT1prob, cuba_config);
            tmpresult += IT1;

            // T2 contribution
            intparams.contribution="T2";
            double IT2, IT2err, IT2prob;
            Cuba(cuba_config.method, 4, integrand_dip_massive, &intparams, &IT2, &IT2err, &IT2prob, cuba_config);
            tmpresult += IT2;            
        }
        else
        {
            throw std::runtime_error("NLODIS::Sigma_dip: unknown polarization");
        }
        //cout << "result " << result << " Quark " << quark.String() << " tmpres " << tmpresult << " chargefact " << SQR(quark.Charge()) << endl;
        result += tmpresult*SQR(quark.Charge()); // Include electric charge factor for this flavor
        
    } // Quark flavor loop

    // We have factorized out \int d^2 b
        
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
    *f=0; // Default result is 0, in case of invalid integration variables

    if (p->pol == Polarization::L)
    {
        if (!( (*ndim ==4 and p->contribution=="cd") or (*ndim == 3 and p->contribution =="ab") 
        or (*ndim ==2 and p->contribution=="Omega_L_const") ) )
        {
            throw std::runtime_error("integrand_dip_massive: ndim " + std::to_string(*ndim) + " and contribution " + p->contribution + " do not match");
        }   
    }
    else
    {
        if (!( (*ndim ==4 and p->contribution=="T2") or (*ndim == 3 and p->contribution =="T1") 
        or (*ndim ==2 and p->contribution=="T0") ) )
        {
            throw std::runtime_error("integrand_dip_massive: ndim " + std::to_string(*ndim) + " and contribution " + p->contribution + " do not match");
        }   
    }

    // Validate integration variables to avoid some Cuba crashes
    // Todo: would be better to understand why Cuba can sometimes sample integration variables to be NaN
    for (int i = 0; i < *ndim; ++i) {
        if (!std::isfinite(x[i]) or x[i] < 0. or x[i] > 1.) {
            //#ifdef DEBUG
            std::cerr << "Warning: integrand_dip_massive: x[" << i << "] = " << x[i] << std::endl;
            //#endif
            return 0;
        }
    }
 
    
    double Q2=p->Q2;
    double xbj=p->xbj;
    double mf=p->quark.mass;

    double z1=x[0];
    double x01=p->nlodis->GetMaxR()*x[1];

    if (x01 < p->nlodis->GetDipole().MinR() or x01 > p->nlodis->GetDipole().MaxR())
    {
        *f=0;
        return 0;
    }

    double x01sq=SQR(x01);
    
    

    double alphabar=p->nlodis->Alphas(x01)*Constants::CF/M_PI;


    double evolution_rapidity = std::log(1/xbj); 
    double dipole = p->nlodis->GetDipole().DipoleAmplitude(x01,evolution_rapidity);
    double res=0;
    /////////////// Longitudinal part ///////////////
    if (p->contribution=="ab" and p->pol == Polarization::L) {
        // "ab" contribution does not have the x integration variable
        double xi=x[2]; 
        res = dipole*(ILdip_massive_Iab(Q2,z1,x01,mf,xi));
    } else if (p->contribution=="cd" and p->pol == Polarization::L) {
        double xi=x[2]; 
        double intx=x[3];
        res = dipole*(ILdip_massive_Icd(Q2,z1,x01,mf,xi,intx));
    }
    else if (p->contribution=="Omega_L_const" and p->pol == Polarization::L)
    {
        // only z and r integration
        res = dipole*ILdip_massive_Omega_L_Const(Q2, z1, x01, mf);
    }
    /////////////// Transverse part ///////////////
    else if (p->contribution=="T0" and p->pol == Polarization::T)
    {
        res = dipole*ITdip_massive_0(Q2, z1, SQR(x01)  , mf);
    }
    else if (p->contribution=="T1" and p->pol == Polarization::T)
    {
        double xi=x[2]; 
        res = dipole*ITdip_massive_1(Q2, z1, SQR(x01), mf, xi);
    }
    else if (p->contribution=="T2" and p->pol == Polarization::T)
    {
        double xi=x[2]; 
        double intx=x[3];
        res = dipole*ITdip_massive_2(Q2, z1, SQR(x01), mf, xi, intx);
    }
    else 
    {
        throw std::runtime_error("integrand_dip_massive: unknown contribution " + p->contribution + " pol " + PolarizationString(p->pol));
    }

    double jacobian = x01 * p->nlodis->GetMaxR(); // Jacobian from d^2r and r = u*MAXR
    res *= jacobian*alphabar; 

    if(std::isfinite(res)){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

/*
 * qg contribution
 * Longitudinal reference: (167) but instead of q^+, k^+ we integrate over z_i
 * Explicit expressoin is docs/NLO_DIS_cross_section_with_massive_quarks.pdf (13) 
*/
double NLODIS::Sigma_qg_d2b(double Q2, double xbj, Polarization pol)
{
    // Note on factors: the transverse integration measures are defined with 1/(2pi), see
    // 2103.14549. This measure is not visible in the note, but should be there. Therefore
    // we have 1/(2pi)^3 below (from d^2x_{01} d^2x_{02} d^2b)
    double fac=4.0*Constants::NC*Constants::AlphaEM/std::pow(2.0*M_PI,3.0);
    

    double result=0;
    for (const auto& quark : quarks) {
        double tmpresult=0; // This quark flavor contribution
        IntegrationParams intparams;
        intparams.nlodis=this;
        intparams.Q2=Q2;
        intparams.xbj=xbj;
        intparams.pol=pol;

        intparams.quark=quark;

        // Sum different contributions, labeling is the same for T and L polarizations,
        // proper integrand is selcted in integrand_qgunsub_massive according to intarams.pol and intparams.contribution

        // Note (21): this contribution is split into 3 parts I_1, I_2 and I_3
        double I1,I1err,I1prob;
        intparams.contribution="I1";
        Cuba(cuba_config.method, 5, integrand_qgunsub_massive, &intparams, &I1, &I1err, &I1prob, cuba_config);
        tmpresult += I1;

        double I2,I2err,I2prob;
        intparams.contribution="I2";
        Cuba(cuba_config.method, 6, integrand_qgunsub_massive, &intparams, &I2, &I2err, &I2prob, cuba_config);
        tmpresult += I2;

        double I3,I3err,I3prob;
        intparams.contribution="I3";
        Cuba(cuba_config.method, 7, integrand_qgunsub_massive, &intparams, &I3, &I3err, &I3prob, cuba_config);
           
        tmpresult += I3;

        result += tmpresult*SQR(quark.Charge());
    }

    


    // We have factorized out (not performed) \int d^2 b 


    // 2pi: overall integral over one angle
    return  fac*2.0*M_PI*result;

}


/*
* Integrand wrapper for Cuba
* To be used to evaluate I2 part of NLO qg unsubtracted contribution
* Integration variables are
* x[0] = z1
* x[1] = z2
* x[2] = x01 
* x[3] = x02 
* x[4] = phi_x0102 (angle between x01 and x02)
* x[5] = y_t (contribution I2 and I3)
* x[6] = y_t2 (contribution I3)
* 
* Note: overall 2pi integral is included in NLODIS::Sigma_qg
*/
 int integrand_qgunsub_massive(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    auto* const p = static_cast<IntegrationParams*>(userdata);

    if (!( // Note: same contribution labels and ndim's for T and L polarization 
        (*ndim == 5 and p->contribution == "I1") or
        (*ndim ==6 and p->contribution=="I2")   or
        (*ndim ==7 and p->contribution=="I3") )
        )
    {
        throw std::runtime_error("integrand_qgunsub_massive: ndim " + std::to_string(*ndim) + " and contribution " + p->contribution + " do not match");
        exit(1);
    }

    // Validate integration variables to avoid some Cuba crashes
    // Todo: would be better to understand why Cuba can sometimes sample integration variables to be NaN
    for (int i = 0; i < *ndim; ++i) {
        if (!std::isfinite(x[i]) or x[1]< 0. or x[1] > 1.) {
            //#ifdef DEBUG
            std::cerr << "Warning: integrand_qgunsub_massive: x[" << i << "] = " << x[i] << std::endl;
            //#endif
            return 0;
        }
    }
    


    double Q2=p->Q2;
    double xbj=p->xbj;
    double mf=p->quark.mass;

   
    double z2min = p->nlodis->z2_lower_bound(xbj,Q2); 
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;

    double x01=p->nlodis->GetMaxR()*x[2];
    double x02=p->nlodis->GetMaxR()*x[3];
    double phix0102=2.0*M_PI*x[4];
    
    double x01sq=SQR(x01);
    double x02sq=SQR(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);

    double alphabar = p->nlodis->Alphas(p->nlodis->RunningCouplinScale(x01, x02, std::sqrt(x21sq))) * Constants::CF / Constants::NC;

    // Validate
    if (x21sq < 0)
    {
        cerr << "Warning: integrand_qgunsub_massive: x21sq is negative: " << x21sq << ". Setting integrand to 0." << endl;
        *f = 0;
        return 0;
    }



    // Jacobians from Cuba variable changes (z's, 2 distances, 1 angle) and d^2x_01 d^2x_02
    // Note: one overall 2pi from angular integral is included in NLODIS::Sigma_qg
    // All other integrals in I1, I2 and I3 are from 0 to 1
    double jac=(1.0-z2min)*(1.0-z1-z2min)*x01*x02 * p->nlodis->GetMaxR()*p->nlodis->GetMaxR()*2.0*M_PI;

    double evolution_rapidity=p->nlodis->EvolutionRapidity(xbj,Q2,z2);

    if (evolution_rapidity < 0){
        cerr << "Warning: integrand_ILqgunsub_massive: evolution rapidity < 0: " << evolution_rapidity << ", xbj=" << xbj << ", Q2=" << Q2 << ", z2=" << z2 << endl;
        *f=0;
        return 0;
    }

    double SKernel_tripole = p->nlodis->TripoleAmplitude(x01,x02, std::sqrt(x21sq), evolution_rapidity);
    double SKernel_dipole = p->nlodis->GetDipole().DipoleAmplitude(x01, evolution_rapidity);

    double alphafac=p->nlodis->Alphas(p->nlodis->RunningCouplinScale(x01,x02,std::sqrt(x21sq)))*Constants::CF/Constants::NC;

    double res=0;

    if (p->contribution == "I1" and p->pol == Polarization::L)
    {
        double dipole_term  = SKernel_dipole  * ILNLOqg_massive_dipole_uvsub(Q2,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_01 = UV subtraction terms
        double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I1(Q2,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_012

        res =  dipole_term + tripole_term ;
    }
    else if (p->contribution == "I1" and p->pol == Polarization::T)
    {
        double dipole_term = SKernel_dipole * ITNLOqg_massive_dipole_uvsub(Q2,mf,z1,z2,x01sq,x02sq,x21sq);
        double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I1(Q2,mf,z1,z2,x01sq,x02sq,x21sq);

        res = dipole_term + tripole_term;
    }
    else  if (p->contribution == "I2" and p->pol == Polarization::L) {
        double y_t1 = x[5];
       res = SKernel_tripole * ILNLOqg_massive_tripole_part_I2(Q2,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); 
    }
    else if (p->contribution == "I2" and p->pol == Polarization::T) {
        double y_t1 = x[5];
        res = SKernel_tripole * ITNLOqg_massive_tripole_part_I2_fast(Q2,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); 
    }
    else if (p->contribution == "I3" and p->pol == Polarization::L)
    {
        double y_t1 = x[5];
        double y_t2 = x[6];
        res = SKernel_dipole * ILNLOqg_massive_tripole_part_I3(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_t2);
    }
    else if (p->contribution == "I3" and p->pol == Polarization::T)
    {
        double y_t1 = x[5];
        double y_t2 = x[6];
        res = SKernel_dipole * ITNLOqg_massive_tripole_part_I3_fast(Q2, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_t2);
    }
    else
    {
        cerr << "integrand_qgunsub_massive: unknown contribution " + p->contribution << " polarization " << PolarizationString(p->pol) << endl;
        exit(1);
    }

    res *=  jac*alphafac/z2;

    if(std::isfinite(res)){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}





 NLODIS::NLODIS()
 {
    // Initialize quarks with default masses
    quarks = std::vector<Quark> { Quark(Quark::Type::LIGHT), Quark(Quark::Type::C) } ;

 }
 
 double NLODIS::Alphas(double r) const
 {
    // TODO: Implement possibility to have r dependent Nf
    int nf=0;
    for (const auto& quark : quarks) {
        if (quark.type == Quark::Type::LIGHT) 
            nf += 3;
        else
            nf += 1;
    }

    const double b0 = (11.0*Constants::NC - 2.0*nf)/(12.0*M_PI);
    const double scalefactor = 4.0*config.C2_alpha; // Convention: 4C^2/r^2 is the scale at which alpha_s is evaluated in coordinate space

    switch (config.rc_ir_scheme)
    {
        case RunningCouplingIRScheme::FREEZE: {
            double mu2 = scalefactor/(r*r);
            double log_arg = mu2/(Constants::LambdaQCD*Constants::LambdaQCD);

            double as = 1.0/(b0*log(log_arg));
            if (as > 0.7 or log_arg < 1.0)
            {
                as = 0.7; // Freeze coupling
            }
            return as;
            break;
        }
        case RunningCouplingIRScheme::SMOOTH: {
            double scalefactor = 4.0*config.C2_alpha;
            const double alphas_mu0=2.5;    // mu0/lqcd
            const double alphas_freeze_c=0.2;

            double as = 1. / ( b0 * std::log(
            std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/SQR(r*Constants::LambdaQCD), 1.0/alphas_freeze_c), alphas_freeze_c)
            ) );
            return as;
            break;
        }
        default:
            throw std::runtime_error("NLODIS::Alphas: unknown IR scheme");
    }
    
   
 }

 double NLODIS::z2_lower_bound(double xbj, double Q2) const
{
    double W2 = Q2 / xbj;
    return NLODISConfig::Q0sqr / W2;
}



void NLODIS::SetQuarkMass(Quark::Type type, double mass)
{
    if (mass <= 0)
    {
        throw std::runtime_error("Quark mass must be positive. Routines for exactly m=0 quarks have not been implemented.");
    }
    bool found = false;
    for (auto& quark : quarks) {
        if (quark.type == type) {
            quark.mass = mass;
            found=true;
        }
    }
    if (!found) {
        throw std::runtime_error("Trying to set mass for the quark flavor " + Quark(type).String() + ", which is not in the quark list. Note: quarks included = " + Quark::QuarkListToString(quarks));
    }
}
void NLODIS::SetProtonTransverseArea(double transverse_area_, Unit unit)
{
    if (transverse_area_ <= 0)
    {
        throw std::runtime_error("Transverse area must be positive");
    }
    if (unit == Unit::MB)
    {
        // Convert from mb to GeV^-2
        transverse_area = transverse_area_ * 2.56819; 
    }
    else if (unit == Unit::GEVm2)
    {
        transverse_area = transverse_area_;
    }
    else
    {
        throw std::runtime_error("Unknown unit in SetProtonTransverseArea");
    }
}

void NLODIS::SetDipole(std::unique_ptr<Dipole> dipole_)
{
    dipole = std::move(dipole_);
}

void NLODIS::PrintConfiguration(const std::string& lineprefix) const
{
    std::cout << lineprefix << "=== NLODIS Configuration Summary ===" << std::endl;
    
    // Order
    std::cout << lineprefix << "Order: ";
    if (config.order == Order::LO) {
        std::cout << "LO" << std::endl;
    } else if (config.order == Order::NLO) {
        std::cout << "NLO" << std::endl;
    } else {
        std::cout << "Unknown" << std::endl;
    }
    
    // Subtraction Scheme
    std::cout << lineprefix << "Subtraction Scheme: ";
    if (config.scheme == SubtractionScheme::UNSUB) {
        std::cout << "Unsubstracted (UNSUB)" << std::endl;
    } else {
        std::cout << "Unknown" << std::endl;
    }
    
    // Nc Scheme
    std::cout << lineprefix << "Nc Scheme: ";
    if (config.nc_scheme == NcScheme::LargeNC) {
        std::cout << "Large Nc" << std::endl;
    } else if (config.nc_scheme == NcScheme::FiniteNC) {
        std::cout << "Finite Nc" << std::endl;
    } else {
        std::cout << "Unknown" << std::endl;
    }
    
    // Running Coupling Scheme
    std::cout << lineprefix << "Running Coupling Scale: ";
    if (config.rc_scheme == RunningCouplingScheme::SMALLEST) {
        std::cout << "Smallest dipole" << std::endl;
    } else if (config.rc_scheme == RunningCouplingScheme::PARENT) {
        std::cout << "Parent dipole" << std::endl;
    } else {
        std::cout << "Unknown" << std::endl;
    }
    
    // IR Freezing Scheme
    std::cout << lineprefix << "IR Freezing Scheme: ";
    if (config.rc_ir_scheme == RunningCouplingIRScheme::FREEZE) {
        std::cout << "Freeze" << std::endl;
    } else if (config.rc_ir_scheme == RunningCouplingIRScheme::SMOOTH) {
        std::cout << "Smooth" << std::endl;
    } else {
        std::cout << "Unknown" << std::endl;
    }
    
    // Numerical parameters
    std::cout << lineprefix << "Maximum dipole size (maxr): " << config.maxr << " GeV^-1" << std::endl;
    std::cout << lineprefix << "Running coupling C^2 factor: " << config.C2_alpha << std::endl;
    std::cout << lineprefix << "Non-perturbative scale (Q0^2): " << NLODISConfig::Q0sqr << " GeV^2" << std::endl;
    
    // Proton parameters
    std::cout << lineprefix << "Proton transverse area: " << transverse_area << " GeV^-2 = " << transverse_area/2.568 << " mb" << std::endl;
    
    // Quark flavors
    std::cout << lineprefix << "Quark flavors and masses:" << std::endl;
    for (const auto& q : quarks) {
        std::cout << lineprefix << "  " << q.String() << ", m=" << q.mass << " GeV (e=" << q.Charge() << ")" << std::endl;
    }
    std::cout << lineprefix << "Dipole: " << (dipole ? dipole->GetString() : "None") << std::endl;
    const char* cores_env = std::getenv("CUBACORES");
    std::string cores = cores_env ? cores_env : "not set (using default)";
    cout << lineprefix << "Cuba integration method: " << cuba_config.method << ", maxeval " << cuba_config.maxeval <<", relaccuracy " << cuba_config.epsrel << " cores " << cores << std::endl;
    std::cout << lineprefix << "===================================\n" << std::endl;
}