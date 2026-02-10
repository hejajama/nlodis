/*
 * NLO DIS code
 *
 * Compute inclusive structure functions at NLO in the dipole picture, with massive quarks.
 * Currently the code supports only targets with impact-parameter independent dipole amplitude,
 * formally corresponding to infinite large
 * 
 * Based on code originally written by H. Hänninen
 * Modified by J. Penttala, H. Mäntysaari
 */

#pragma once

#include "dipole/dipoleamplitude.hpp"
#include "qcd.hpp"
#include "datatypes.hpp"
#include <vector>
#include <string>
#include <memory>
#include <gsl/gsl_integration.h>

class NLODIS
{
    public:
        
        NLODIS();

        /**
         * @brief Structure function F2 
         * @param Q2 Photon virtuality [GeV^2].
         * @param xbj Bjorken-x.
         * 
         * @return F2
         */
        double F2(double Q2, double xbj);
        /**
         * @brief Structure function FL.
         * @param Q2 Photon virtuality [GeV^2].
         * @param xbj Bjorken-x.
         */
        double FL(double Q2, double xbj);
        /**
         * @brief Structure function FT.
         * @param Q2 Photon virtuality [GeV^2].
         * @param xbj Bjorken-x.
         */
        double FT(double Q2, double xbj);
        /** 
         * @brief Photon-proton cross section (without the impact parameter integral)
         * 
         * Calculate the photon-proton cross section without the impact parameter integral,
         * at LO or NLO depending on the order set by SetOrder().
         * 
         * The result is _not_ integrated over the impact parameter, but factor 2 from
         * the optical theorem is included.
         * 
         * In the case of proton target, the cross section can be obtained by multiplying 
         * this by \sigma_0/2 = \int d^2 b = ProtonTransverseArea()
         * 
         * @param Q2 photon virtuality [GeV^2]
         * @param xbj Bjorken-x
         * @param pol photon polarization (T or L)
         * 
         * @return d\sigma(\gamma+A)/d^2b. 
         */
        double Photon_proton_cross_section_d2b(double Q2, double xbj, Polarization pol);

        /** 
         * @brief Photon-proton cross section at LO
         * 
         * Calculate the photon-proton cross section without the impact parameter integral,
         * at LO
         * 
         * The result is _not_ integrated over the impact parameter, but factor 2 from
         * the optical theorem is included.
         * 
         * In the case of proton target, the cross section can be obtained by multiplying 
         * this by \sigma_0/2 = \int d^2 b = ProtonTransverseArea()
         * 
         * @param Q2 photon virtuality [GeV^2]
         * @param xbj Bjorken-x
         * @param pol photon polarization (T or L)
         * 
         * @return d\sigma(\gamma+A)/d^2b. 
         */
        double Photon_proton_cross_section_LO_d2b(double Q2, double xbj, Polarization pol);

        /**
         * @brief \sigma_dip: qq part of the NLO cross section.
         * 
         * NLO contribution from the diagrams where qq interacts with the shockwave
         * Result is not integrated over the impact parameter, but factor 2 from the 
         * optical theorem is included.
         * 
         * @param Q2 Photon virtuality [GeV^2].
         * @param xbj Bjorken-x.
         * @param pol Photon polarization.
         * 
         * @return d\sigma_dip/d^2b (dimensionless)
         */
        double Sigma_dip_d2b(double Q2, double xbj, Polarization pol);
        /**
         * @brief \sigma_qg: qg part of the NLO cross section.
         * 
         * NLO contribution from the diagrams where qqg interacts with the shockwave
         * Result is not integrated over the impact parameter, but factor 2 from the
         *  optical theorem is included.
         * 
         * @param Q2 Photon virtuality [GeV^2].
         * @param xbj Bjorken-x.
         * @param pol Photon polarization.
         * 
         * @return d\sigma_qg/d^2b (dimensionless)
         */
        double Sigma_qg_d2b(double Q2, double xbj, Polarization pol);

        void SetDipole(std::unique_ptr<Dipole> dipole_);
        
        void SetOrder(Order o) { order = o; }
        double GetMaxR() const { return maxr; }
        
        /**
         * @brief Coordinate space coupling.
         * 
         * TODO: implement n_f dependence
         * 
         * @param r Dipole size [GeV^-1].
         */
        double Alphas(double r) const;

        /**
         * @brief Set proton transverse area = \sigma_0/2
         * 
         * @param transverse_area Proton transverse area \int d^2 b
         * @param unit Unit of transverse_area, default is GeV^-2. If unit is MB, the value will be converted to GeV^-2 internally.
         * 
         */
        void SetProtonTransverseArea(double transverse_area_, Unit unit=Unit::GEVm2);

        /**
         * @brief Proton transverse area in GeV^-2
         */
        double ProtonTransverseArea() const { return transverse_area; }

        Dipole& GetDipole() const { return *dipole; }

        /**
         * @brief Lower bound for the z2 integral.
         * @param xbj Bjorken-x.
         * @param Q2 Photon virtuality [GeV^2].
         * @return Lower bound for z2.
         *
         * Ref https://arxiv.org/pdf/2007.01645 eq (18)
         */
        double z2_lower_bound(double xbj, double Q2) const;

        /**
         * @brief qqg-target scattering amplitude.
         *
         * Ref. https://arxiv.org/pdf/2211.03504 (4). In that notation, this is 1-S_{012}.
         *
         * @param x01 Dipole size [GeV^-1].
         * @param x02 Dipole size [GeV^-1].
         * @param x21 Dipole size [GeV^-1].
         * @param Y Evolution rapidity.
         */
        double TripoleAmplitude(double x01, double x02, double x21, double Y);

        /**
         * @brief Evolution rapidity (in the qqg contribution).
         * @param xbj Bjorken-x.
         * @param Q2 Photon virtuality [GeV^2].
         * @param z2 Gluon longitudinal momentum fraction.
         *
         * Ref https://arxiv.org/pdf/2007.01645 eq (19)
         */
        double EvolutionRapidity(double xbj, double Q2, double z2) const;

        void SetNcScheme(NcScheme scheme_) { nc_scheme = scheme_; }
        void SetRunningCouplingScheme(RunningCouplingScheme rc_) { rc_scheme = rc_; }
        /**
         * @brief Running coupling scale depending on the RC scheme used.
         * @param x01 Dipole size [GeV^-1].
         * @param x02 Dipole size [GeV^-1].
         * @param x21 Dipole size [GeV^-1].
         */
        double RunningCouplinScale(double x01, double x02, double x21) const;

        /**
         * @brief Set scale factor C^2 in the coordinate space running coupling
         * 
         * \alpha_s(r) = \alpha_s(4C^2/r^2 Lambda_QCD^2)
         * 
         * @param c2 Scale factor C^2
         */
        void SetRunningCouplingC2Alpha(double c2) { C2_alpha = c2; }

        void SetQuarks(const std::vector<Quark>& quark_list) { quarks = quark_list; }

        /**
         * @brief Set mass of the given quark flavor.
         * @param type Quark flavor.
         * @param mass Quark mass [GeV].
         */
        void SetQuarkMass(Quark::Type type, double mass);

        /**
         * @brief Control alpha_s in IR
         * 
         * In the FREEZE scheme, the coupling is frozen to a constant value when the scale is below a certain value.
         * In the SMOOTH scheme, the coupling smoothly freezes to a constant value in the IR, without a sharp cutoff.
         */
        void SetRunningCouplingIRScheme(RunningCouplingIRScheme rc_ir_scheme_) { rc_ir_scheme = rc_ir_scheme_; }
    private:
    
        /* Integrand for LO cross section
        * |\Psi|^2 N(r): integrand for the LO cross section, to be integrated over r and z
        * Note: Jacobian of r is _NOT_ included here
        * The dipole-target cross section is sigma_0 \int d^2r dz Integrand_photon_target_LO
        * The factor \sigma_0 = 2 \int d^2 b is NOT included here
        */
        double Integrand_photon_target_LO(double r, double z, double x, double Q2, Polarization pol );
    
        double transverse_area=1; // \sigma_0/2 = \int d^2b in GeV^-2 (proton transverse area)
        //AmplitudeLib dipole;
        std::unique_ptr<Dipole> dipole;
        SubtractionScheme scheme = SubtractionScheme::UNSUB;
        std::vector<Quark> quarks;
        Order order = Order::LO;
        double maxr = 99;
        NcScheme nc_scheme = NcScheme::LargeNC;
        RunningCouplingScheme rc_scheme = RunningCouplingScheme::SMALLEST;
        const double Q0sqr = 1; // Non-perturbative target scale, should match the one used in the NLO DIS fit!
        double C2_alpha = 1.0; // Scale factor in the coordinate space running coupling
        RunningCouplingIRScheme rc_ir_scheme = RunningCouplingIRScheme::FREEZE;
      
};

// Data structure used to store integration parameters
struct IntegrationParams {
        NLODIS* nlodis;
        double Q2;
        double xbj;
        double z;
        Polarization pol;
        gsl_integration_workspace* w_r;
        Quark quark;
        std::string contribution;
    };

constexpr double SQR(double x) noexcept { return x*x; }

//// TODO NOTE: Inconsistency: some functions take r^2, some functions take r as an argumetn. This should be fixed at some point to avoid confusion.

// Helper functions in nlodishelper.cpp that need to be accessed outside
// e.g. for unit tests
int integrand_dip_massive(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata);
double ILdip_massive_Icd(double Q2, double z1, double r, double mf, double xi, double x); 
double ILdip_massive_Iab(double Q2, double z1, double r, double mf, double xi);
double ILdip_massive_Omega_L_Const(double Q2, double z1, double r, double mf);

// Sigma_dip transverse
double ITdip_massive_0(double Q2, double z1, double x01sq, double mf);
double ITdip_massive_1(double Q2, double z1, double x01sq, double mf, double y_chi);
double ITdip_massive_2(double Q2, double z1, double x01sq, double mf, double y_chi, double y_u);

// Sigma_qg longitudinal part helper
int integrand_qgunsub_massive(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata);
double ILNLOqg_massive_tripole_part_I2_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double G_integrand_simplified(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda, double y);
double ILNLOqg_massive_dipole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I3_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);

// Sigma qg T
double ITNLOqg_massive_dipole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ITNLOqg_massive_tripole_part_I1(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_dipole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ITNLOqg_massive_tripole_part_I2_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double ITNLOqg_massive_tripole_part_I3_fast(double Q2, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);

// Helper functions from nlodishelper_transverse.cpp
double OmegaT_V_unsymmetric(double Q, double z, double mf);
double OmegaT_N_unsymmetric(double Q, double z, double mf);
double IT_V1_unsymmetric(double Q, double z, double mf, double r, double xi);
double IT_VMS1_unsymmetric(double Q, double z, double mf, double r, double xi);
double IT_V2_unsymmetric(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_VMS2_unsymmetric(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_N_unsymmetric(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_V1(double Q, double z, double mf, double r, double xi);
double IT_VMS1(double Q, double z, double mf, double r, double xi);
double IT_V2(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_VMS2(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_N(double Q, double z, double mf, double r, double y_chi, double y_u);
double IT_tripole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_dipole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_F_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_Fm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_jk_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double IT_tripole_jkm_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double IT_tripole_F_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double IT_tripole_Fm_I2_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t);
double IT_tripole_jkm_I3_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);
double IT_tripole_jk_I3_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);
double IT_tripole_F_I3_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2);
double IT_tripole_Fm_I3_fast(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_t2); 
