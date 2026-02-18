#include "qcd.hpp"
#include <cmath>
#include <sstream>

// Default constructor
Quark::Quark() : type(LIGHT), mass(0.001) {}

// Constructor that automatically sets mass and charge based on type
Quark::Quark(Type t) : type(t) {
    switch(type) {
        case U:
        case D:
        case S:
        case LIGHT:
            mass = 0.002; 
            break;
        case C:
            mass = 1.4;  
            break;
        case B:
            mass = 4.18; 
            break;
        case T:
            mass = 173.0;
            break;
    }
}

// Constructor for custom mass/charge (backward compatibility)
Quark::Quark(Type t, double m) : type(t), mass(m) {}

double Quark::Charge() const
{
    switch (type)
    {
        case U: return 2.0/3.0;
        case D: return -1.0/3.0;
        case S: return -1.0/3.0;
        case C: return 2.0/3.0;
        case B: return -1.0/3.0;
        case LIGHT: return std::sqrt(std::pow(2.0/3.0, 2) + 2*std::pow(-1.0/3.0, 2)); // Effective charge for u+d+s
        default: 
             throw std::runtime_error("Unknown quark flavor ");
    }
}



std::string Quark::String() const
{
    switch(type) {
        case U: return "u";
        case D: return "d";
        case S: return "s";
        case C: return "c";
        case B: return "b";
        case T: return "t";
        case LIGHT: return "light";
        default: return "unknown";
    }
}

std::string Quark::QuarkListToString(const std::vector<Quark>& quarks)
{
    std::stringstream ss;
    for (size_t i = 0; i < quarks.size(); ++i) {
        ss << quarks[i].String();
        if (i < quarks.size() - 1) {
            ss << ", ";
        }
    }
    return ss.str();
}