#ifndef VEC2_HPP
#define VEC2_HPP

#include <cmath>

struct Vec2
{
    double x = 0.0;
    double y = 0.0;

    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2 &b) const { return Vec2(x + b.x, y + b.y); }
    Vec2 operator-(const Vec2 &b) const { return Vec2(x - b.x, y - b.y); }
    Vec2 operator*(double s)      const { return Vec2(x * s, y * s); }

    double dot(const Vec2 &b) const { return x * b.x + y * b.y; }
    double length()           const { return std::sqrt(x * x + y * y); }
};

inline Vec2 operator*(double s, const Vec2 &v) { return v * s; }

#endif
