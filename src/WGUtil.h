#ifndef _WG_UTIL_
#define _WG_UTIL_

#include <cmath>

const double g_double_epsilon = 1E-12;
const double g_unit_epsilon = 1E-6;
const double g_pi = 3.141592653589793;

inline bool is_zero(double d, double epsilon) {
    return d <= epsilon && d >= -epsilon;
}

inline double acos_safe(double d) {
    if (d >= 1) {
        return 0;
    }
    if (d <= -1) {
        return g_pi;
    }
    return acos(d);
}

#endif
