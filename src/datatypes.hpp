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