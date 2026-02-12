/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2020
 */

#pragma once

#include <vector>
#include <string>

/**
 * Read BK equation solution from a given file.
 */
class DataFile
{
    public:
        /**
         * Read BK equation solution from the given file
         */
        explicit DataFile(const std::string& fname);
        
        double MinR() const;
        double RMultiplier() const;
        int RPoints() const;
        
        double MaxY() const;
        double Y0() const;

        /**
         * Return dipole amplitudes and evolution rapidities
         *
         * @return pair of (amplitudes, rapidities)
         */
        std::pair<std::vector<std::vector<double>>, std::vector<double>> GetData() const;

    private:
        std::string filename;
        std::vector<std::vector<double>> data;
        std::vector<double> yvals;
        double minr = 0.0;
        double r_multiplier = 0.0;
        int rpoints = 0;
        double y0 = 0.0;
};
