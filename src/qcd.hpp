#pragma once

#include <string>

// Basic constants and structures

const double ALPHA_EM = 1.0/137.035999084; // Fine-structure constant

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


const int NC = 3; 
const double CF = (NC*NC - 1.0)/(2.0*NC); 
