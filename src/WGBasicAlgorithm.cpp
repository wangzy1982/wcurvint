#include "WGBasicAlgorithm.h"
#include "WGUtil.h"

bool intersect_beeline_plane(const WGVector3d& beeline_origin, const WGVector3d& beeline_direction,
    const WGVector3d& plane_origin, const WGVector3d& plane_normal, WGVector3d& intersection_point) {
    double d = beeline_direction.Dot(plane_normal);
    if (is_zero(d, g_unit_epsilon)) {
        return false;
    }
    double t = (plane_origin - beeline_origin).Dot(plane_normal) / d;
    intersection_point = beeline_origin + beeline_direction * t;
    return true;
}

double get_point_beeline_nearest(const WGVector3d& point, const WGVector3d& beeline_origin, const WGVector3d& beeline_direction) {
    return (point - beeline_origin).Dot(beeline_direction) / beeline_direction.Dot(beeline_direction);
}

int intersect_beeline_circle(const WGVector2d& beeline_origin, const WGVector2d& beeline_direction,
    const WGVector2d& circle_center, double circle_radius, double distance_epsilon, WGVector2d* result_buffer) {
    double d = (beeline_origin - circle_center).Cross(beeline_direction);
    double c = abs(d);
    if (c < circle_radius) {
        double a = sqrt(circle_radius * circle_radius - c * c);
        WGVector2d p = circle_center + WGVector2d(beeline_direction.Y, -beeline_direction.X) * d;
        WGVector2d vt = beeline_direction * a;
        result_buffer[0] = p - vt;
        result_buffer[1] = p + vt;
        return 2;
    }
    if (c <= circle_radius + distance_epsilon) {
        result_buffer[0] = WGVector2d(beeline_direction.Y, -beeline_direction.X) * d;
        return 1;
    }
    return 0;
}

double calculate_tangent_cone_angle(const WSInterval& x, const WSInterval& y) {
    if (x.Min > 0) {
        if (y.Min > 0) {
            WGVector2d vt0(x.Max, y.Min);
            WGVector2d vt1(x.Min, y.Max);
            return vt0.AngleBetween(vt1);
        }
        else if (y.Max < 0) {
            WGVector2d vt0(x.Min, y.Min);
            WGVector2d vt1(x.Max, y.Max);
            return vt0.AngleBetween(vt1);
        }
        else {
            WGVector2d vt0(x.Min, y.Min);
            WGVector2d vt1(x.Min, y.Max);
            return vt0.AngleBetween(vt1);
        }
    }
    else if (x.Max < 0) {
        if (y.Min > 0) {
            WGVector2d vt0(x.Max, y.Max);
            WGVector2d vt1(x.Min, y.Min);
            return vt0.AngleBetween(vt1);
        }
        else if (y.Max < 0) {
            WGVector2d vt0(x.Min, y.Max);
            WGVector2d vt1(x.Max, y.Min);
            return vt0.AngleBetween(vt1);
        }
        else {
            WGVector2d vt0(x.Max, y.Max);
            WGVector2d vt1(x.Max, y.Min);
            return vt0.AngleBetween(vt1);
        }
    }
    else {
        if (y.Min > 0) {
            WGVector2d vt0(x.Max, y.Min);
            WGVector2d vt1(x.Min, y.Min);
            return vt0.AngleBetween(vt1);
        }
        else if (y.Max < 0) {
            WGVector2d vt0(x.Min, y.Max);
            WGVector2d vt1(x.Max, y.Max);
            return vt0.AngleBetween(vt1);
        }
        else {
            return g_pi * 2;
        }
    }
}