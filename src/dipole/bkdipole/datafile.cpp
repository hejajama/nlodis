/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */


#include "datafile.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
using std::ifstream;
using std::getline;
using std::stringstream;
using std::string;
using std::cerr;
using std::endl;

#define LINEINFO __FILE__ << ":" << __LINE__

DataFile::DataFile(string fname)
{
    filename=fname;
    ifstream file(fname.c_str());
    if (!file.is_open())
    {
        cerr <<  "ERROR! Coudn't read BK solution file " << fname << " " << LINEINFO << endl ;
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
                    minr = std::stod(line.substr(3,line.length()-3));
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
                    cerr << "File " << fname << " is formatted incorrectly!" << endl;
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

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != rpoints)
            {
                cerr << "File " << fname << ": read " << tmpvec.size() << " points, but "
                << "there should have been " << rpoints << " points, y=" 
                << y << ". " << LINEINFO << endl;
            }

            y = y0 + std::stod(line.substr(3,line.length()-3));
            yvals.push_back(y);
            tmpvec.clear();
            continue;   // Next line is probably new amplitude value
        }
        else if (line.substr(0,1)=="#")
            continue;   // Comment

        // Ok, so this a new amplitude value
        if (line.length()==0)
            continue;   // Skip empty lines
        tmpvec.push_back(std::stod(line));
    }

    // Add last entry
    data.push_back(tmpvec);
    
    file.close();

    if (data[0].size() != rpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " rpoints, but "
        << "there should have been " << rpoints << " points???" << LINEINFO << endl;
    }
}

void DataFile::GetData(std::vector< std::vector<double> > &n,
                        std::vector<double> &rapidities)
{
    n.clear();
    rapidities.clear();
    // Return vector where indexes are vec[y][r] containing amplitude

    for (uint yind=0; yind < data.size(); yind++)
    {
        std::vector<double> tmpvec;
        for (uint rind=0; rind<rpoints; rind++)
        {
            tmpvec.push_back(data[yind][rind]);
        }
        n.push_back(tmpvec);

        rapidities.push_back(yvals[yind]);
        
    }
    
}

double DataFile::MinR()
{
    return minr;
}

double DataFile::RMultiplier()
{
    return r_multiplier;
}

int DataFile::RPoints()
{
        return rpoints;
}

double DataFile::MaxY()
{
    return yvals[yvals.size()-1];
}

double DataFile::Y0()
{
    return y0;
}
