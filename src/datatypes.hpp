#pragma once

enum class Polarization
{
    T,
    L
};

enum class Unit
{
    MB,     // millibarn
    GEVm2 // GeV^-2
};


/*
 * Convert polarization enum to string for error messages
 */
std::string PolarizationString(Polarization pol);


enum class SubtractionScheme
{
    UNSUB
};

enum class Order
{
    LO,
    NLO
};

enum class NcScheme
{
    FiniteNC,
    LargeNC
};

enum class RunningCouplingScheme
{
    SMALLEST,
    PARENT
};

enum class RunningCouplingIRScheme
{
    FREEZE,
    SMOOTH
};

namespace Constants
{
    constexpr double LambdaQCD = 0.241; // GeV
    constexpr double AlphaEM = 1.0/137.035999084; // Fine-structure constant
    constexpr int NC = 3;
    constexpr double CF = (NC*NC - 1.0)/(2.0*NC);
};