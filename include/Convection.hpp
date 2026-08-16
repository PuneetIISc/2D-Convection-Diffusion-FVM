#ifndef CONVECTION_HPP
#define CONVECTION_HPP

#include "Vec2.hpp"

class ConvectionScheme
{
public:
    virtual ~ConvectionScheme() {}
    // phiCell : cell value; grad : cell gradient; d : face_centre - cell_centroid
    virtual double state(double phiCell, const Vec2 &grad, const Vec2 &d) const = 0;
};

class FirstOrder : public ConvectionScheme
{
public:
    double state(double phiCell, const Vec2 &, const Vec2 &) const
    { return phiCell; }
};

class SecondOrder : public ConvectionScheme
{
public:
    double state(double phiCell, const Vec2 &grad, const Vec2 &d) const
    { return phiCell + grad.x * d.x + grad.y * d.y; }
};

#endif
