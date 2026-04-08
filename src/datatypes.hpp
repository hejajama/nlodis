/**
 * @file datatypes.hpp
 * @brief Defines enumerations and constants used throughout the NLO DIS calculation framework.
 */
#pragma once

#include <string>

/**
 * @enum Polarization
 * @brief Represents the polarization state of the virtual photon.
 * @var Polarization::T
 *      Transverse polarization.
 * @var Polarization::L
 *      Longitudinal polarization.
 */
enum class Polarization
{
    T,
    L
};

/**
 * @brief Converts a Polarization enum value to its string representation.
 * @param pol The polarization enum value to convert.
 * @return A string representation of the polarization ("T" or "L").
 */
std::string PolarizationString(Polarization pol);

/**
 * @enum Unit
 * @brief Unit for transverse area
 * @var Unit::MB
 *      Millibarn.
 * @var Unit::GEVm2
 *      GeV^-2 (inverse GeV squared).
 */
enum class Unit
{
    MB,     // millibarn
    GEVm2 // GeV^-2
};


/**
 * @enum SubtractionScheme
 * @brief Specifies the subtraction scheme used in NLO calculations.
 * 
 * Currently only the unsubtracted scheme (UNSUB) is implemented
 * 
 * @var SubtractionScheme::UNSUB
 *      Unsubtracted scheme.
 */
enum class SubtractionScheme
{
    UNSUB
};

/**
 * @enum SigmaDipScheme
 * @brief Specifies the scheme used for the dipole sigma (dip) contributions.
 * @var SigmaDipScheme::AnalyticalZ2Int
 *      The z2 integration is performed analytically in the dip contributions.
 * @var SigmaDipScheme::ExplicitZ2int
 *      The z2 integration is performed explicitly, allowing for a z2-dependent
 *      evolution rapidity. Implements the scheme described in ref. 2112.08818 sec. 3.3.3
 *      for light quarks. Not currently implemented in this program.
 */
enum class SigmaDipScheme
{
    AnalyticalZ2Int,            
    ExplicitZ2int
};

/**
 * @enum HeavyQuarkX
 * @brief Controls whether the momentum fraction x used to determine the evolution
 *        rapidity depends on the quark mass.
 * @var HeavyQuarkX::MassDependentX
 *      The momentum fraction x is defined as x = (Q^2 + 4m^2) / (W^2 + Q^2),
 *      where m is the quark mass.
 * @var HeavyQuarkX::MassIndependentX
 *      The momentum fraction x is defined as x = x_Bj = Q^2 / (W^2 + Q^2),
 *      independent of the quark mass.
 */
enum class HeavyQuarkX
{
    MassDependentX,
    MassIndependentX
};


/**
 * @enum Order
 * @brief Specifies the perturbative order of the calculation.
 * @var Order::LO
 *      Leading order.
 * @var Order::NLO
 *      Next-to-leading order.
 */
enum class Order
{
    LO,
    NLO
};

/**
 * @enum NcScheme
 * @brief Specifies the treatment of the number of colors Nc.
 * @var NcScheme::FiniteNC
 *      Finite Nc scheme, keeping exact Nc dependence.
 * @var NcScheme::LargeNC
 *      Large Nc limit approximation.
 */
enum class NcScheme
{
    FiniteNC,
    LargeNC
};

/**
 * @enum RunningCouplingScheme
 * @brief Specifies the scheme for the running coupling scale setting.
 * @var RunningCouplingScheme::SMALLEST
 *      The scale is set by the smallest dipole size in the process.
 * @var RunningCouplingScheme::PARENT
 *      The scale is set by the parent dipole size.
 */
enum class RunningCouplingScheme
{
    SMALLEST,
    PARENT
};

/**
 * @enum RunningCouplingIRScheme
 * @brief Specifies the infrared (IR) regularization scheme for the running coupling.
 * @var RunningCouplingIRScheme::FREEZE
 *      The running coupling is frozen below a certain IR scale.
 * @var RunningCouplingIRScheme::SMOOTH
 *      The running coupling is smoothly regulated in the IR.
 */
enum class RunningCouplingIRScheme
{
    FREEZE,
    SMOOTH
};

/**
 * @namespace Constants
 * @brief Contains fundamental physical constants used in the calculations.
 * @var Constants::LambdaQCD
 *      QCD scale parameter Lambda_QCD in GeV (0.241 GeV).
 * @var Constants::AlphaEM
 *      Fine-structure constant (electromagnetic coupling), alpha_EM = 1/137.035999
 * @var Constants::NC
 *      Number of colors, Nc = 3.
 * @var Constants::CF
 *      Casimir factor of the fundamental representation, C_F = (Nc^2 - 1) / (2 * Nc).
 */
namespace Constants
{
    constexpr double LambdaQCD = 0.241; // GeV
    constexpr double AlphaEM = 1.0/137.035999; 
    constexpr int NC = 3;
    constexpr double CF = (NC*NC - 1.0)/(2.0*NC);
};
