#include "WGVector2d.h"
#include <assert.h>
#include "WGUtil.h"

WGVector2d::WGVector2d() :
    X(0), 
    Y(0) {
}

WGVector2d::WGVector2d(double x, double y) :
    X(x),
    Y(y) {
}

double WGVector2d::Length() const {
    return sqrt(X * X + Y * Y);
}

double WGVector2d::SqrLength() const {
    return X * X + Y * Y;
}

double WGVector2d::Normalize(double epsilon) {
    double d = Length();
    if (d > epsilon) {
        X /= d;
        Y /= d;
    }
    return d;
}

bool WGVector2d::IsNormalized() const {
    return abs(X * X + Y * Y - 1) <= g_double_epsilon;
}

double WGVector2d::Dot(const WGVector2d& other) const {
    return X * other.X + Y * other.Y;
}

double WGVector2d::Cross(const WGVector2d& other) const {
    return X * other.Y - Y * other.X;
}

double WGVector2d::UnitAngle() const {
    assert(is_zero(SqrLength() - 1, g_double_epsilon));
    double a = acos_safe(X);
    if (Y < 0) {
        a = g_pi * 2 - a;
    }
    return a;
}

double WGVector2d::AngleTo(const WGVector2d& other) const {
    assert(is_zero(SqrLength() - 1, g_double_epsilon));
    assert(is_zero(other.SqrLength() - 1, g_double_epsilon));
    double a = acos_safe(Dot(other));
    if (Cross(other) < 0) {
        a = g_pi * 2 - a;
    }
    return a;
}

double WGVector2d::AngleBetween(const WGVector2d& other) const {
    double d2 = SqrLength() * other.SqrLength();
    if (d2 == 1) {
        return acos_safe(Dot(other));
    }
    double d = sqrt(d2);
    if (d <= g_double_epsilon) {
        return 0;
    }
    return acos_safe(Dot(other) / d);
}

WGVector2d operator+(const WGVector2d& a, const WGVector2d& b) {
    return WGVector2d(a.X + b.X, a.Y + b.Y);
}

WGVector2d operator-(const WGVector2d& a, const WGVector2d& b) {
    return WGVector2d(a.X - b.X, a.Y - b.Y);
}

WGVector2d operator-(const WGVector2d& a) {
    return WGVector2d(-a.X, -a.Y);
}

WGVector2d operator*(const WGVector2d& a, double b) {
    return WGVector2d(a.X * b, a.Y * b);
}

WGVector2d operator*(double a, const WGVector2d& b) {
    return WGVector2d(a * b.X, a * b.Y);
}

WGVector2d operator/(const WGVector2d& a, double b) {
    return WGVector2d(a.X / b, a.Y / b);
}

int find_vector2d(const std::vector<WGVector2d>& vts, const WGVector2d& vt, double epsilon) {
    for (int i = 0; i < (int)vts.size(); ++i) {
        const WGVector2d& vt2 = vts.at(i);
        if (is_zero(vt2.X - vt.X, epsilon) && is_zero(vt2.Y - vt.Y, epsilon)) {
            return i;
        }
    }
    return -1;
}
