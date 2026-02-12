#pragma once

#include <string>

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
        T
    } type;
    double mass; // in GeV
    double charge;
};

double QuarkCharge(Quark q);
std::string QuarkString(Quark q);


