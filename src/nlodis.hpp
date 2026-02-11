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



/**
 * @brief Configuration parameters for NLODIS calculations
 */
struct NLODISConfig {
    Order order = Order::LO;                                                ///< LO or NLO
    SubtractionScheme scheme = SubtractionScheme::UNSUB;                   ///< Subtraction scheme
    NcScheme nc_scheme = NcScheme::FiniteNC;                                 ///< Nc scheme
    RunningCouplingScheme rc_scheme = RunningCouplingScheme::SMALLEST;     ///< Running coupling scale choice
    RunningCouplingIRScheme rc_ir_scheme = RunningCouplingIRScheme::FREEZE; ///< IR freezing scheme for coupling
    double maxr = 99.0;                                                     ///< Maximum dipole size [GeV^-1]
    double C2_alpha = 1.0;                                                  ///< Scale factor C^2 in coordinate space running coupling
    static constexpr double Q0sqr = 1.0;                                   ///< Non-perturbative target scale [GeV^2]
};

/**
 * @class NLODIS
 * @brief Next-to-Leading Order Deep Inelastic Scattering in Dipole Picture.
 * 
 * This class provides methods for calculating structure functions and cross sections
 * in the color glass condensate framework for deep inelastic scattering processes.
 * It supports both leading order (LO) and next-to-leading order (NLO) calculations.
 * 
 */
class NLODIS
{
    public:
        
        NLODIS();

        /**
         * @brief Structure function F2 
         * 
         * Calculates the total F2 structure function, F2 = FL + FT, by combining
         * longitudinal and transverse photon contributions.
         * 
         * @param Q2 Photon virtuality [GeV^2]
         * @param xbj Bjorken-x
         * 
         * @return F2 structure function (dimensionless)
         */
        double F2(double Q2, double xbj);
        
        /**
         * @brief Longitudinal structure function FL
         * 
         * Calculates the contribution from longitudinally polarized virtual photons.
         * 
         * @param Q2 Photon virtuality [GeV^2]
         * @param xbj Bjorken-x
         * 
         * @return FL structure function (dimensionless)
         */
        double FL(double Q2, double xbj);
        
        /**
         * @brief Transverse structure function FT
         * 
         * Calculates the contribution from transversely polarized virtual photons.
         * Note: F2 = FL + FT
         * 
         * @param Q2 Photon virtuality [GeV^2]
         * @param xbj Bjorken-x
         * 
         * @return FT structure function (dimensionless)
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

        /**
         * @brief Set the dipole amplitude model to use for calculations
         * 
         * This method sets the dipole amplitude object that will be used to calculate
         * the dipole-target scattering amplitude N
         * The NLODIS object takes ownership of the passed pointer via 
         * std::unique_ptr, which ensures automatic cleanup when the NLODIS object is
         * destroyed or when SetDipole is called again.
         * 
         * @param dipole_ Unique pointer to a Dipole object. The lifetime of the Dipole object
         *               is managed by NLODIS after this call.
         * 
         * @note This method must be called before performing any structure function
         *       calculations (F2, FL, FT) or cross section calculations. Attempting
         *       calculations without setting a dipole model will result in dereferencing
         *       a null pointer and undefined behavior.
         * 
         * @see Dipole, BKDipole
         */
        void SetDipole(std::unique_ptr<Dipole> dipole_);
        
        /**
         * @brief Set the perturbative order for calculations
         * 
         * @param o Order::LO for leading order or Order::NLO for next-to-leading order
         */
        void SetOrder(Order o) { config.order = o; }
        
        /**
         * @brief Get maximum dipole size used in integrations
         * 
         * @return Maximum dipole size in GeV^-1
         */
        double GetMaxR() const noexcept { return config.maxr; }
        
        /**
         * @brief Coordinate space running coupling constant
         * 
         * Calculates the running coupling α_s evaluated at the scale 4C²/r².
         * The number of active flavors n_f is determined by the number of quarks
         * in the quark list (set via SetQuarks).
         * 
         * The coupling can freeze in the IR according to the scheme set by
         * SetRunningCouplingIRScheme().
         * 
         * @param r Dipole size [GeV^-1]
         * 
         * @return α_s(4C²/r²) (dimensionless)
         */
        double Alphas(double r) const;

        /**
         * @brief Set proton transverse area (σ₀/2)
         * 
         * Sets the proton transverse area ∫d²b, which is used to convert the
         * differential cross sections d²σ/d²b to total cross sections.
         * 
         * @param transverse_area Proton transverse area ∫d²b
         * @param unit Unit of transverse_area. Default is GeV^-2. If Unit::MB is specified,
         *             the value will be converted to GeV^-2 internally using the conversion
         *             factor 1 mb = 2.5684624 GeV^-2
         */
        void SetProtonTransverseArea(double transverse_area_, Unit unit=Unit::GEVm2);

        /**
         * @brief Get proton transverse area
         * 
         * @return Proton transverse area σ₀/2 = ∫d²b in GeV^-2
         */
        double ProtonTransverseArea() const noexcept { return transverse_area; }


        /**
         * @brief Lower bound for the gluon momentum fraction z₂ integral
         * 
         * In NLO calculations, the gluon longitudinal momentum fraction z₂ must
         * satisfy z₂ > Q₀²/W², where W² is the photon-proton center-of-mass energy
         * squared. This ensures the dipole amplitude is evaluated at positive rapidity.
         * 
         * @param xbj Bjorken-x
         * @param Q2 Photon virtuality [GeV^2]
         * 
         * @return Lower bound for z₂ = Q₀²/W² (dimensionless)
         *
         * @see EvolutionRapidity()
         * 
         * Ref: https://arxiv.org/pdf/2007.01645 eq (18)
         */
        double z2_lower_bound(double xbj, double Q2) const;

        /**
         * @brief qqg-target scattering amplitude (tripole)
         *
         * Calculates the scattering amplitude for a quark-antiquark-gluon (qqg) state
         * interacting with the target. This corresponds to 1-S₀₁₂ in the notation of
         * the reference, where S₀₁₂ is the S-matrix element for the three-parton state.
         * 
         * The implementation accounts for different Nc schemes (large Nc vs finite Nc).
         *
         * @param x01 Dipole size between quark and antiquark [GeV^-1]
         * @param x02 Dipole size between quark and gluon [GeV^-1]
         * @param x21 Dipole size between antiquark and gluon [GeV^-1]
         * @param Y Evolution rapidity
         * 
         * @return Tripole amplitude 1-S₀₁₂ (dimensionless)
         * 
         * Ref: https://arxiv.org/pdf/2211.03504 eq (4)
         */
        double TripoleAmplitude(double x01, double x02, double x21, double Y);

        /**
         * @brief Evolution rapidity for NLO qqg contribution
         * 
         * Calculates the rapidity Y at which the dipole/tripole amplitude should be
         * evaluated in the NLO qqg contribution. The rapidity depends on the photon-proton
         * center-of-mass energy W² = Q²/xbj and the gluon momentum fraction z₂.
         * 
         * Y = ln(W²z₂/Q₀²), where Q₀² is the non-perturbative scale.
         * 
         * @param xbj Bjorken-x
         * @param Q2 Photon virtuality [GeV^2]
         * @param z2 Gluon longitudinal momentum fraction
         *
         * @return Evolution rapidity Y = ln(W²z₂/Q₀²) (dimensionless)
         * 
         * @see z2_lower_bound()
         * 
         * Ref: https://arxiv.org/pdf/2007.01645 eq (19)
         */
        double EvolutionRapidity(double xbj, double Q2, double z2) const;

        /**
         * @brief Set Nc counting scheme for tripole amplitude
         * 
         * @param scheme_ NcScheme::LargeNC for large Nc approximation,
         *                NcScheme::FiniteNC for finite Nc (Nc=3) treatment
         */
        void SetNcScheme(NcScheme scheme_) { config.nc_scheme = scheme_; }
        
        /**
         * @brief Set running coupling scale choice for tripole/multi-dipole configurations
         * 
         * @param rc_ RunningCouplingScheme::SMALLEST to use the smallest dipole size,
         *            RunningCouplingScheme::PARENT to use the parent dipole size
         */
        void SetRunningCouplingScheme(RunningCouplingScheme rc_) { config.rc_scheme = rc_; }
        
        /**
         * @brief Get running coupling scale for tripole configuration
         * 
         * Returns the scale at which α_s should be evaluated based on the running
         * coupling scheme (smallest dipole or parent dipole).
         * 
         * @param x01 Dipole size between quark and antiquark [GeV^-1]
         * @param x02 Dipole size between quark and gluon [GeV^-1]
         * @param x21 Dipole size between antiquark and gluon [GeV^-1]
         * 
         * @return Dipole size to use for α_s evaluation [GeV^-1]
         */
        double RunningCouplinScale(double x01, double x02, double x21) const;

        /**
         * @brief Set scale factor C² in coordinate space running coupling
         * 
         * The running coupling in coordinate space is evaluated at the scale
         * μ² = 4C²/r², where r is the dipole size. This method sets the C² factor.
         * 
         * α_s(r) = α_s(4C²/r²)
         * 
         * @param c2 Scale factor C² (dimensionless, parametrically ~1)
         */
        void SetRunningCouplingC2(double c2) { config.C2_alpha = c2; }

        /**
         * @brief Set list of active quark flavors and their masses
         * 
         * @param quark_list Vector of Quark objects containing flavor type, mass, and charge
         * 
         * @see SetQuarkMass()
         */
        void SetQuarks(const std::vector<Quark>& quark_list) { quarks = quark_list; }

        /**
         * @brief Set mass for a specific quark flavor
         * 
         * Updates the mass of an existing quark in the quark list. The quark must
         * already exist in the list (added via SetQuarks).
         * 
         * @param type Quark flavor (Quark::U, Quark::D, Quark::S, Quark::C, Quark::B, Quark::T)
         * @param mass Quark mass [GeV]
         * 
         * @see SetQuarks()
         */
        void SetQuarkMass(Quark::Type type, double mass);

        /**
         * @brief Set infrared behavior of running coupling α_s
         * 
         * Controls how the running coupling behaves at small scales (large dipole sizes).
         * 
         * - FREEZE: Coupling is frozen to a constant value (typically 0.7) when the
         *           argument of the logarithm becomes too small
         * - SMOOTH: Coupling smoothly freezes using a prescription without sharp cutoff
         * 
         * @param rc_ir_scheme_ RunningCouplingIRScheme::FREEZE or RunningCouplingIRScheme::SMOOTH
         * 
         * @see Alphas()
         */
        void SetRunningCouplingIRScheme(RunningCouplingIRScheme rc_ir_scheme_) { config.rc_ir_scheme = rc_ir_scheme_; }

        /**
         * @brief Print detailed configuration summary to stdout
         * 
         * Displays all configuration parameters including:
         * - Calculation order (LO/NLO)
         * - Subtraction scheme
         * - Nc counting scheme
         * - Running coupling settings
         * - Numerical parameters (maxr, C², Q₀²)
         * - Proton transverse area
         * - Active quark flavors and masses
         */
        void PrintConfiguration() const;

        /**
         * @brief Get reference to dipole amplitude object (non-const)
         * 
         * @return Non-const reference to Dipole object
         */
         Dipole& GetDipole() { return *dipole; }
         
        /**
         * @brief Get const reference to dipole amplitude object
         * 
         * @return Const reference to Dipole object
         */
        const Dipole& GetDipole() const { return *dipole; }
    private:
    
        /* Integrand for LO cross section
        * |\Psi|^2 N(r): integrand for the LO cross section, to be integrated over r and z
        * Note: Jacobian of r is _NOT_ included here
        * The dipole-target cross section is sigma_0 \int d^2r dz Integrand_photon_target_LO
        * The factor \sigma_0 = 2 \int d^2 b is NOT included here
        */
        double Integrand_photon_target_LO(double r, double z, double x, double Q2, Polarization pol );
    
        double transverse_area = 1.0;                  ///< \sigma_0/2 = \int d^2b in GeV^-2 (proton transverse area)
        std::unique_ptr<Dipole> dipole;                ///< Dipole amplitude model
        std::vector<Quark> quarks;                     ///< Quark flavors and masses
        NLODISConfig config;                           ///< Configuration parameters
      
};

/**
 * @brief Custom deleter for GSL integration workspace
 */
struct GslIntegrationWorkspaceDeleter {
    void operator()(gsl_integration_workspace* w) const noexcept {
        if (w) gsl_integration_workspace_free(w);
    }
};

/**
 * @brief Parameters passed to integration routines
 */
struct IntegrationParams {
    NLODIS* nlodis;                                                  ///< Pointer to NLODIS instance
    double Q2;                                                        ///< Photon virtuality [GeV^2]
    double xbj;                                                       ///< Bjorken-x
    double z;                                                         ///< Longitudinal momentum fraction
    Polarization pol;                                                 ///< Photon polarization
    std::unique_ptr<gsl_integration_workspace, GslIntegrationWorkspaceDeleter> w_r; ///< Integration workspace
    Quark quark;                                                      ///< Current quark flavor
    std::string contribution;                                         ///< Integration contribution name
};

constexpr double SQR(double x) noexcept { return x*x; }

//// Note:
// Internal helper functions typically take Q [GeV] or squared distances as arguments
// Everything visible to outside from NLODIS always uses Q^2

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
