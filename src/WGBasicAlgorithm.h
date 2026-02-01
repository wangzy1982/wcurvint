#ifndef _WG_BASIC_ALGORITHM_
#define _WG_BASIC_ALGORITHM_

#include "WGVector2d.h"
#include "WGVector3d.h"
#include "wsolver.h"

bool intersect_beeline_plane(const WGVector3d& beeline_origin, const WGVector3d& beeline_direction,
    const WGVector3d& plane_origin, const WGVector3d& plane_normal, WGVector3d& intersection_point);

double get_point_beeline_nearest(const WGVector3d& point, const WGVector3d& beeline_origin, const WGVector3d& beeline_direction);

int intersect_beeline_circle(const WGVector2d& beeline_origin, const WGVector2d& beeline_direction, 
    const WGVector2d& circle_center, double circle_radius, double distance_epsilon, WGVector2d* result_buffer);

double calculate_tangent_cone_angle(const WSInterval& x, const WSInterval& y);

#endif
