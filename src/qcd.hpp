#pragma once

#include <string>
#include <vector>

// Basic constants and structures

struct Quark
{
    enum Type
    {
        U,
        D,
        S,
        C,
        B,
        T,
        LIGHT // U+D+S, with an effective electric charge sqrt((2/3)^2 + 2*(-1/3)^2) = sqrt(2/3)
    } type;
    double mass;   // in GeV

    // Constructors
    Quark();                              // Default constructor
    explicit Quark(Type t);               // Auto-set mass/charge based on type
    Quark(Type t, double m);    // Custom mass/charge

    // Member functions
    double Charge() const;
    std::string String() const;
    static std::string QuarkListToString(const std::vector<Quark>& quarks);

};



