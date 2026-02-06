#include "qcd.hpp"
#include <iostream>
#include <string>

double QuarkChage(Quark q)
{
    switch (q.type)
    {
        case Quark::Type::U: return 2.0/3.0;
        case Quark::Type::D: return -1.0/3.0;
        case Quark::Type::S: return -1.0/3.0;
        case Quark::Type::C: return 2.0/3.0;
        case Quark::Type::B: return -1.0/3.0;
        default: 
             throw std::runtime_error("Unknown quark flavor ");
    }
}


std::string QuarkString(Quark q)
{
    switch (q.type)
    {
        case Quark::Type::U: return "u";
        case Quark::Type::D: return "d";
        case Quark::Type::S: return "s";
        case Quark::Type::C: return "c";
        case Quark::Type::B: return "b";
        default: 
             throw std::runtime_error("Unknown quark flavor ");
    }
}