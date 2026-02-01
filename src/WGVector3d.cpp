#include "WGVector3d.h"
#include <assert.h>
#include "WGUtil.h"

WGVector3d::WGVector3d() :
    X(0),
    Y(0),
    Z(0) {
}

WGVector3d::WGVector3d(double x, double y, double z) :
    X(x),
    Y(y),
    Z(z) {
}

double WGVector3d::Length() const {
    return sqrt(X * X + Y * Y + Z * Z);
}

double WGVector3d::SqrLength() const {
    return X * X + Y * Y + Z * Z;
}

double WGVector3d::Normalize(double epsilon) {
    double d = Length();
    if (d > epsilon) {
        X /= d;
        Y /= d;
        Z /= d;
    }
    return d;
}

bool WGVector3d::IsNormalized() const {
    return abs(X * X + Y * Y + Z * Z - 1) <= g_double_epsilon;
}

double WGVector3d::Dot(const WGVector3d& other) const {
    return X * other.X + Y * other.Y + Z * other.Z;
}

WGVector3d WGVector3d::Cross(const WGVector3d& other) const {
    return WGVector3d(Y * other.Z - Z * other.Y, Z * other.X - X * other.Z, X * other.Y - Y * other.X);
}

WGVector3d operator+(const WGVector3d& a, const WGVector3d& b) {
    return WGVector3d(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
}

WGVector3d operator-(const WGVector3d& a, const WGVector3d& b) {
    return WGVector3d(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
}

WGVector3d operator-(const WGVector3d& a) {
    return WGVector3d(-a.X, -a.Y, -a.Z);
}

WGVector3d operator*(const WGVector3d& a, double b) {
    return WGVector3d(a.X * b, a.Y * b, a.Z * b);
}

WGVector3d operator*(double a, const WGVector3d& b) {
    return WGVector3d(a * b.X, a * b.Y, a * b.Z);
}

WGVector3d operator/(const WGVector3d& a, double b) {
    return WGVector3d(a.X / b, a.Y / b, a.Z / b);
}
