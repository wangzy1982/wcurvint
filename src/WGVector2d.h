#ifndef _WG_VECTOR_2D_
#define _WG_VECTOR_2D_

#include <vector>

class WGVector2d {
public:
    WGVector2d();
    WGVector2d(double x, double y);
    double Length() const;
    double SqrLength() const;
    double Normalize(double epsilon);
    bool IsNormalized() const;
    double Dot(const WGVector2d& other) const;
    double Cross(const WGVector2d& other) const;
    double UnitAngle() const;
    double AngleTo(const WGVector2d& other) const;
    double AngleBetween(const WGVector2d& other) const;
public:
    double X;
    double Y;
};

WGVector2d operator+(const WGVector2d& a, const WGVector2d& b);
WGVector2d operator-(const WGVector2d& a, const WGVector2d& b);
WGVector2d operator-(const WGVector2d& a);
WGVector2d operator*(const WGVector2d& a, double b);
WGVector2d operator*(double a, const WGVector2d& b);
WGVector2d operator/(const WGVector2d& a, double b);

int find_vector2d(const std::vector<WGVector2d>& vts, const WGVector2d& vt, double epsilon);

#endif
