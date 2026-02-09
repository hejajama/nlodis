#include <string>
#include "datatypes.hpp"

std::string PolarizationString(Polarization pol)
{
    if (pol == T)
    {
        return "T";
    }
    else if (pol == L)
    {
        return "L";
    }
    else
    {
        throw std::runtime_error("Unknown polarization in PolarizationString.");
    }
}