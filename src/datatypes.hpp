#ifndef _DATATYPES_H
#define _DATATYPES_H

enum Polarization
{
    T,
    L
};

enum Unit
{
    MB,     // millibarn
    GEVm2 // GeV^-2
};


/*
 * Convert polarization enum to string for error messages
 */
std::string PolarizationString(Polarization pol);


enum Scheme
{
    UNSUB
};

enum Order
{
    LO,
    NLO
};

enum NcScheme
{
    FiniteNC,
    LargeNC
};

enum RunningCouplingScheme
{
    SMALLEST,
    PARENT
};

#endif