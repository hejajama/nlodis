#include <string>
#include <stdexcept>
#include "datatypes.hpp"

std::string PolarizationString(Polarization pol)
{
    if (pol == Polarization::T)
    {
        return "T";
    }
    else if (pol == Polarization::L)
    {
        return "L";
    }
    else
    {
        throw std::runtime_error("Unknown polarization in PolarizationString.");
    }
}