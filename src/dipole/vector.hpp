#pragma once

/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 

#include <iostream>
#include <vector>

typedef double REAL;

class Vec
{
    public:
        Vec();
        Vec(REAL x_, REAL y_);
        Vec(REAL x_, REAL y_, REAL z_);
        Vec (const Vec &v);
        REAL GetX() const;
        REAL GetY() const ;
        REAL GetZ() const;
        void SetX(REAL x_);
        void SetY(REAL y_);
        void SetZ(REAL z_);
        Vec& operator+=(const Vec& v);
        Vec& operator-=(const Vec& v);
        Vec& operator=(const Vec& v);
        Vec operator+(const Vec& v) const;
        Vec operator-(const Vec& v) const;
        Vec& operator*=(REAL c);
        Vec operator*(REAL c) const;
    
        void Rotate2D(double angle);  // rotate x,y coords counterclockwise
    
        double operator*(const Vec &v) const;   // dot product
    
        REAL Len() const;
        REAL LenSqr() const;
    
    private:
        REAL x,y,z;
};

std::ostream& operator<<(std::ostream& os, const Vec& ic);

