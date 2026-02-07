#ifndef _WCURVINT_
#define _WCURVINT_

#if defined(_WINDOWS)
    #if defined(WCURVINT_EXPORTS)
        #define WCURVINT_API __declspec(dllexport)
    #elif defined(WCURVINT_STATIC)
        #define WCURVINT_API
    #else
        #define WCURVINT_API __declspec(dllimport)
    #endif
#else
    #define WCURVINT_API
#endif

#define WCURVINT_API_C extern "C" WCURVINT_API

WCURVINT_API_C void* new_line2d();
WCURVINT_API_C void set_line2d_start_point(void* curve, double x, double y);
WCURVINT_API_C void set_line2d_end_point(void* curve, double x, double y);

WCURVINT_API_C void* new_arc2d();
WCURVINT_API_C void set_arc2d_center(void* curve, double x, double y);
WCURVINT_API_C void set_arc2d_radius(void* curve, double radius);
WCURVINT_API_C void set_arc2d_start_angle(void* curve, double start_angle);
WCURVINT_API_C void set_arc2d_delta_angle(void* curve, double delta_angle);

WCURVINT_API_C void* new_bezier_curve2d(int degree);
WCURVINT_API_C void set_bezier_curve2d_control_point(void* curve, int index, double x, double y);

WCURVINT_API_C void* new_rational_bezier_curve2d(int degree);
WCURVINT_API_C void set_rational_bezier_curve2d_control_point(void* curve, int index, double x, double y);
WCURVINT_API_C void set_rational_bezier_curve2d_weight(void* curve, int index, double w);

WCURVINT_API_C void calculate_curve_point(void* curve, double t, double& x, double& y);

WCURVINT_API_C void free_curve2d(void* curve);

WCURVINT_API_C void* new_curve_curve_intersections(void* curve0, void* curve1, double distance_epsilon);

WCURVINT_API_C void free_curve_curve_intersections(void* intersections);

WCURVINT_API_C void is_curve_curve_intersections_singularity(void* intersections, bool& is_singularity0, bool& is_singularity1);

WCURVINT_API_C int get_curve_curve_intersection_count(void* intersections);

WCURVINT_API_C void get_curve_curve_intersection(void* intersections, int index, bool& is_overlap,
    double& start_t0, double& start_t1, bool& is_start_sample, double& end_t0, double& end_t1, bool& is_end_sample);

#endif