/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */

#include "datafile.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#define LINEINFO __FILE__ << ":" << __LINE__

using namespace std;


DataFile::DataFile(const std::string& fname)
{
    filename=fname;
    ifstream file(fname.c_str());
    if (!file.is_open())
    {
        throw runtime_error("ERROR! Coudn't read BK solution file " + fname);
        exit(1);
    }
    int confid=0;
    double tmp;
    while(!file.eof() and confid < 4)
    {
        string line;
        getline(file, line);
        if (line.substr(0, 3)=="###")
        {                    
            switch (confid)
            {
                case 0:
                    minr = std::stod(line.substr(3));
                    break;
                case 1:
                    r_multiplier = std::stod(line.substr(3,line.length()-3));
                    break;
                case 2:
                    rpoints = std::stoi(line.substr(3,line.length()-3));
                    break;
                case 3:
                    tmp = std::stod(line.substr(3,line.length()-3));
                    y0 = std::log(1/tmp);
                    if (y0 < 0 and y0 > -1e-10)
                        y0=0;
                    
                    if (y0 < 0 )
                    {
                        cerr <<"Warning: Invalid Y0 = " << y0 <<" " << LINEINFO << ". Using Y0=0!" << endl;
                        y0=0;
                    }
                    break;
                default:
                    throw runtime_error("File " + fname + " is formatted incorrectly!");
                    break;
            }
            confid++; 
        }
    }
    // Ok, configurations are read, then read all yvals
    double y=-1;
    std::vector<double> tmpvec;
    while (!file.eof())
    {
        string line;
        getline(file, line);

        if (line.empty())
            continue;

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != rpoints)
            {
                throw runtime_error("File " + fname + ": read " + std::to_string(tmpvec.size()) + " points, but "
                + "there should have been " + std::to_string(rpoints) + " points, y=" 
                + std::to_string(y) + ". ");
            }
            
            y = y0 + std::stod(line.substr(3));

            yvals.push_back(y);
            tmpvec.clear();
            continue;
        }
        else if (line.substr(0,1)=="#")
            continue;   // Comment
        else // new amplitude value
            tmpvec.push_back(std::stod(line));
    }

    // Add last entry
    data.push_back(tmpvec);
    
    file.close();

    if (data[0].size() != rpoints)
    {
        throw runtime_error("File " + fname + ": read " + std::to_string(data.size()) + " rpoints, but "
        + "there should have been " + std::to_string(rpoints) + " points! ");
    }
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> DataFile::GetData() const
{
    return {data, yvals};
}

double DataFile::MinR() const
{
    return minr;
}

double DataFile::RMultiplier() const
{
    return r_multiplier;
}

int DataFile::RPoints() const
{
    return rpoints;
}

double DataFile::MaxY() const
{
    return yvals.empty() ? 0.0 : yvals.back();
}

double DataFile::Y0() const
{
    return y0;
}
