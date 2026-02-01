#ifndef _WG_ARC_2D_
#define _WG_ARC_2D_

#include <assert.h>

class WGArc2d {
public:
    static double CalculateT(double start_angle, double delta_angle, double angle) {
        assert(angle >= 0 && angle < g_pi * 2);
        double end_angle = start_angle + delta_angle;
        if (delta_angle > 0) {
            double d = angle;
            if (d < start_angle) {
                d += g_pi * 2;
            }
            if (d <= end_angle) {
                return (d - start_angle) / delta_angle;
            }
            else {
                if (d <= start_angle + g_pi) {
                    return (d - start_angle) / delta_angle;
                }
                else {
                    return (d - g_pi * 2 - start_angle) / delta_angle;
                }
            }
        }
        else {
            assert(delta_angle < 0);
            double d = angle;
            if (d > start_angle) {
                d -= g_pi * 2;
            }
            if (d >= end_angle) {
                return (d - start_angle) / delta_angle;
            }
            else {
                if (d >= start_angle - g_pi) {
                    return (d - start_angle) / delta_angle;
                }
                else {
                    return (d + g_pi * 2 - start_angle) / delta_angle;
                }
            }
        }
    }

    static WGVector2d CalculatePoint(const WGVector2d& center, double radius, double angle) {
        return WGVector2d(center.X + radius * cos(angle), center.Y + radius * sin(angle));
    }

};

#endif
