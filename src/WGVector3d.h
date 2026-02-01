#ifndef _WG_VECTOR_3D_
#define _WG_VECTOR_3D_

class WGVector3d {
public:
    WGVector3d();
    WGVector3d(double x, double y, double z);
    double Length() const;
    double SqrLength() const;
    double Normalize(double epsilon);
    bool IsNormalized() const;
    double Dot(const WGVector3d& other) const;
    WGVector3d Cross(const WGVector3d& other) const;
public:
    double X;
    double Y;
    double Z;
};

WGVector3d operator+(const WGVector3d& a, const WGVector3d& b);
WGVector3d operator-(const WGVector3d& a, const WGVector3d& b);
WGVector3d operator-(const WGVector3d& a);
WGVector3d operator*(const WGVector3d& a, double b);
WGVector3d operator*(double a, const WGVector3d& b);
WGVector3d operator/(const WGVector3d& a, double b);

#endif
