#include "wcurvint.h"
#include "WGIntersectHelper2d.h"
#include <vector>

class WCurvInt {
public:
    enum class CurveType {
        Line2d,
        Arc2d,
        BezierCurve2d,
        RationalBezierCurve2d
    };

    struct Curve2d {
        CurveType Type;
    };

    struct Line2d : public Curve2d {
        WGVector2d StartPoint;
        WGVector2d EndPoint;
    };

    struct Arc2d : public Curve2d {
        WGVector2d Center;
        double Radius;
        double StartAngle;
        double DeltaAngle;
    };

    struct BezierCurve2d : public Curve2d {
        int Degree;
        WGVector2d ControlPoints[16];
    };

    struct RationalBezierCurve2d : public Curve2d {
        int Degree;
        WGVector2d ControlPoints[16];
        double Weights[16];
    };

    struct Intersection2ds {
        bool IsSingularity0;
        bool IsSingularity1;
        std::vector<WGIntersectHelper2d::CurveCurveIntersection> Intersections;
    };

    static void ExchangeCurve(Intersection2ds* intersections) {
        bool b = intersections->IsSingularity0;
        intersections->IsSingularity0 = intersections->IsSingularity1;
        intersections->IsSingularity1 = b;
        for (auto& intersection : intersections->Intersections) {
            if (intersection.PointCount == 1) {
                double t = intersection.Ts0[0];
                intersection.Ts0[0] = intersection.Ts1[0];
                intersection.Ts1[0] = t;
            }
            else {
                double t = intersection.Ts0[0];
                intersection.Ts0[0] = intersection.Ts1[0];
                intersection.Ts1[0] = t;
                t = intersection.Ts0[1];
                intersection.Ts0[1] = intersection.Ts1[1];
                intersection.Ts1[1] = t;
                if (intersection.Ts0[0] > intersection.Ts0[1]) {
                    t = intersection.Ts0[0];
                    intersection.Ts0[0] = intersection.Ts0[1];
                    intersection.Ts0[1] = t;
                    t = intersection.Ts1[0];
                    intersection.Ts1[0] = intersection.Ts1[1];
                    intersection.Ts1[1] = t;
                    b = intersection.IsSamples[0];
                    intersection.IsSamples[0] = intersection.IsSamples[1];
                    intersection.IsSamples[1] = b;
                }
            }
        }
    }
};

void* new_line2d() {
    WCurvInt::Line2d* curve = new WCurvInt::Line2d();
    curve->Type = WCurvInt::CurveType::Line2d;
    return curve;
}

void set_line2d_start_point(void* curve, double x, double y) {
    ((WCurvInt::Line2d*)curve)->StartPoint = WGVector2d(x, y);
}

void set_line2d_end_point(void* curve, double x, double y) {
    ((WCurvInt::Line2d*)curve)->EndPoint = WGVector2d(x, y);
}

void* new_arc2d() {
    WCurvInt::Arc2d* curve = new WCurvInt::Arc2d();
    curve->Type = WCurvInt::CurveType::Arc2d;
    return curve;
}

void set_arc2d_center(void* curve, double x, double y) {
    ((WCurvInt::Arc2d*)curve)->Center = WGVector2d(x, y);
}

void set_arc2d_radius(void* curve, double radius) {
    ((WCurvInt::Arc2d*)curve)->Radius = radius;
}

void set_arc2d_start_angle(void* curve, double start_angle) {
    ((WCurvInt::Arc2d*)curve)->StartAngle = start_angle;
}

void set_arc2d_delta_angle(void* curve, double delta_angle) {
    ((WCurvInt::Arc2d*)curve)->DeltaAngle = delta_angle;
}

void* new_bezier_curve2d(int degree) {
    if (degree >= 16) {
        return nullptr;
    }
    WCurvInt::BezierCurve2d* curve = new WCurvInt::BezierCurve2d();
    curve->Type = WCurvInt::CurveType::BezierCurve2d;
    curve->Degree = degree;
    return curve;
}

void set_bezier_curve2d_control_point(void* curve, int index, double x, double y) {
    ((WCurvInt::BezierCurve2d*)curve)->ControlPoints[index] = WGVector2d(x, y);
}

void* new_rational_bezier_curve2d(int degree) {
    if (degree >= 16) {
        return nullptr;
    }
    WCurvInt::RationalBezierCurve2d* curve = new WCurvInt::RationalBezierCurve2d();
    curve->Type = WCurvInt::CurveType::RationalBezierCurve2d;
    curve->Degree = degree;
    return curve;
}

void set_rational_bezier_curve2d_control_point(void* curve, int index, double x, double y) {
    ((WCurvInt::RationalBezierCurve2d*)curve)->ControlPoints[index] = WGVector2d(x, y);
}

void set_rational_bezier_curve2d_weight(void* curve, int index, double w) {
    ((WCurvInt::RationalBezierCurve2d*)curve)->Weights[index] = w;
}

void free_curve2d(void* curve) {
    delete (WCurvInt::Curve2d*)curve;
}

void* new_curve_curve_intersections(void* curve0, void* curve1, double distance_epsilon) {
    switch (((WCurvInt::Curve2d*)curve0)->Type) {
    case WCurvInt::CurveType::Line2d: {
            WCurvInt::Line2d* line0 = (WCurvInt::Line2d*)curve0;
            switch (((WCurvInt::Curve2d*)curve1)->Type) {
            case WCurvInt::CurveType::Line2d: {
                    WCurvInt::Line2d* line1 = (WCurvInt::Line2d*)curve1;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::CurveCurveIntersection intersections[2];
                    int n = WGIntersectHelper2d::LineLineIntersect(line0->StartPoint, line0->EndPoint,
                        line1->StartPoint, line1->EndPoint, distance_epsilon, 
                        result->IsSingularity0, result->IsSingularity1, intersections);
                    result->Intersections.reserve(n);
                    for (int i = 0; i < n; ++i) {
                        result->Intersections.push_back(intersections[i]);
                    }
                    return result;
                }
            case WCurvInt::CurveType::Arc2d: {
                    WCurvInt::Arc2d* arc1 = (WCurvInt::Arc2d*)curve1;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::CurveCurveIntersection intersections[10];
                    int n = WGIntersectHelper2d::LineArcIntersect(line0->StartPoint, line0->EndPoint,
                        arc1->Center, arc1->Radius, arc1->StartAngle, arc1->DeltaAngle, distance_epsilon, 
                        result->IsSingularity0, result->IsSingularity1, intersections);
                    result->Intersections.reserve(n);
                    for (int i = 0; i < n; ++i) {
                        result->Intersections.push_back(intersections[i]);
                    }
                    return result;
                }
            case WCurvInt::CurveType::BezierCurve2d: {
                    WCurvInt::BezierCurve2d* bezier_curve1 = (WCurvInt::BezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::LineBezierCurveIntersect(line0->StartPoint, line0->EndPoint,
                        bezier_curve1->Degree, bezier_curve1->ControlPoints, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            case WCurvInt::CurveType::RationalBezierCurve2d: {
                    WCurvInt::RationalBezierCurve2d* rational_bezier_curve1 = (WCurvInt::RationalBezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::LineRationalBezierCurveIntersect(line0->StartPoint, line0->EndPoint,
                        rational_bezier_curve1->Degree, rational_bezier_curve1->ControlPoints, rational_bezier_curve1->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            default: {
                    return nullptr;
                }
            }
        }
    case WCurvInt::CurveType::Arc2d: {
            WCurvInt::Arc2d* arc0 = (WCurvInt::Arc2d*)curve0;
            switch (((WCurvInt::Curve2d*)curve1)->Type) {
            case WCurvInt::CurveType::Line2d: {
                    WCurvInt::Line2d* line1 = (WCurvInt::Line2d*)curve1;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::CurveCurveIntersection intersections[10];
                    int n = WGIntersectHelper2d::LineArcIntersect(line1->StartPoint, line1->EndPoint,
                        arc0->Center, arc0->Radius, arc0->StartAngle, arc0->DeltaAngle, distance_epsilon,
                        result->IsSingularity0, result->IsSingularity1, intersections);
                    result->Intersections.reserve(n);
                    for (int i = 0; i < n; ++i) {
                        result->Intersections.push_back(intersections[i]);
                    }
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::Arc2d: {
                    WCurvInt::Arc2d* arc1 = (WCurvInt::Arc2d*)curve1;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::CurveCurveIntersection intersections[10];
                    int n = WGIntersectHelper2d::ArcArcIntersect(arc0->Center, arc0->Radius, arc0->StartAngle, arc0->DeltaAngle,
                        arc1->Center, arc1->Radius, arc1->StartAngle, arc1->DeltaAngle, distance_epsilon,
                        result->IsSingularity0, result->IsSingularity1, intersections);
                    result->Intersections.reserve(n);
                    for (int i = 0; i < n; ++i) {
                        result->Intersections.push_back(intersections[i]);
                    }
                    return result;
                }
            case WCurvInt::CurveType::BezierCurve2d: {
                    WCurvInt::BezierCurve2d* bezier_curve1 = (WCurvInt::BezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache; 
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::ArcBezierCurveIntersect(arc0->Center, arc0->Radius, arc0->StartAngle, arc0->DeltaAngle,
                        bezier_curve1->Degree, bezier_curve1->ControlPoints, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            case WCurvInt::CurveType::RationalBezierCurve2d: {
                    WCurvInt::RationalBezierCurve2d* rational_bezier_curve1 = (WCurvInt::RationalBezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::ArcRationalBezierCurveIntersect(arc0->Center, arc0->Radius, arc0->StartAngle, arc0->DeltaAngle,
                        rational_bezier_curve1->Degree, rational_bezier_curve1->ControlPoints, rational_bezier_curve1->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            default: {
                    return nullptr;
                }
            }
        }
    case WCurvInt::CurveType::BezierCurve2d: {
            WCurvInt::BezierCurve2d* bezier_curve0 = (WCurvInt::BezierCurve2d*)curve0;
            switch (((WCurvInt::Curve2d*)curve1)->Type) {
            case WCurvInt::CurveType::Line2d: {
                    WCurvInt::Line2d* line1 = (WCurvInt::Line2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::LineBezierCurveIntersect(line1->StartPoint, line1->EndPoint,
                        bezier_curve0->Degree, bezier_curve0->ControlPoints, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::Arc2d: {
                    WCurvInt::Arc2d* arc1 = (WCurvInt::Arc2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::ArcBezierCurveIntersect(arc1->Center, arc1->Radius, arc1->StartAngle, arc1->DeltaAngle,
                        bezier_curve0->Degree, bezier_curve0->ControlPoints, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::BezierCurve2d: {
                    WCurvInt::BezierCurve2d* bezier_curve1 = (WCurvInt::BezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::BezierCurveBezierCurveIntersect(bezier_curve0->Degree, bezier_curve0->ControlPoints,
                        bezier_curve1->Degree, bezier_curve1->ControlPoints, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            case WCurvInt::CurveType::RationalBezierCurve2d: {
                    WCurvInt::RationalBezierCurve2d* rational_bezier_curve1 = (WCurvInt::RationalBezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::BezierCurveRationalBezierCurveIntersect(bezier_curve0->Degree, bezier_curve0->ControlPoints,
                        rational_bezier_curve1->Degree, rational_bezier_curve1->ControlPoints, rational_bezier_curve1->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            default: {
                    return nullptr;
                }
            }
        }
    case WCurvInt::CurveType::RationalBezierCurve2d: {
            WCurvInt::RationalBezierCurve2d* rational_bezier_curve0 = (WCurvInt::RationalBezierCurve2d*)curve0;
            switch (((WCurvInt::Curve2d*)curve1)->Type) {
            case WCurvInt::CurveType::Line2d: {
                    WCurvInt::Line2d* line1 = (WCurvInt::Line2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::LineRationalBezierCurveIntersect(line1->StartPoint, line1->EndPoint,
                        rational_bezier_curve0->Degree, rational_bezier_curve0->ControlPoints, rational_bezier_curve0->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::Arc2d: {
                    WCurvInt::Arc2d* arc1 = (WCurvInt::Arc2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::ArcRationalBezierCurveIntersect(arc1->Center, arc1->Radius, arc1->StartAngle, arc1->DeltaAngle,
                        rational_bezier_curve0->Degree, rational_bezier_curve0->ControlPoints, rational_bezier_curve0->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::BezierCurve2d: {
                    WCurvInt::BezierCurve2d* bezier_curve1 = (WCurvInt::BezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::BezierCurveRationalBezierCurveIntersect(bezier_curve1->Degree, bezier_curve1->ControlPoints,
                        rational_bezier_curve0->Degree, rational_bezier_curve0->ControlPoints, rational_bezier_curve0->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    WCurvInt::ExchangeCurve(result);
                    return result;
                }
            case WCurvInt::CurveType::RationalBezierCurve2d: {
                    WCurvInt::RationalBezierCurve2d* rational_bezier_curve1 = (WCurvInt::RationalBezierCurve2d*)curve1;
                    WGIntersectHelper2d::IntersectCache cache;
                    WCurvInt::Intersection2ds* result = new WCurvInt::Intersection2ds;
                    WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveIntersect(
                        rational_bezier_curve0->Degree, rational_bezier_curve0->ControlPoints, rational_bezier_curve0->Weights,
                        rational_bezier_curve1->Degree, rational_bezier_curve1->ControlPoints, rational_bezier_curve1->Weights, distance_epsilon, &cache,
                        result->IsSingularity0, result->IsSingularity1, result->Intersections);
                    return result;
                }
            default: {
                    return nullptr;
                }
            }
        }
    default: {
            return nullptr;
        }
    }
}

void free_curve_curve_intersections(void* intersections) {
    delete (WCurvInt::Intersection2ds*)intersections;
}

void is_curve_curve_intersections_singularity(void* intersections, bool& is_singularity0, bool& is_singularity1) {
    is_singularity0 = ((WCurvInt::Intersection2ds*)intersections)->IsSingularity0;
    is_singularity1 = ((WCurvInt::Intersection2ds*)intersections)->IsSingularity1;
}

int get_curve_curve_intersection_count(void* intersections) {
    return (int)((WCurvInt::Intersection2ds*)intersections)->Intersections.size();
}

void get_curve_curve_intersection(void* intersections, int index, bool& is_overlap,
    double& start_t0, double& start_t1, bool& is_start_sample, double& end_t0, double& end_t1, bool& is_end_sample) {
    const WGIntersectHelper2d::CurveCurveIntersection& intersection = ((WCurvInt::Intersection2ds*)intersections)->Intersections[index];
    is_overlap = intersection.PointCount == 2;
    start_t0 = intersection.Ts0[0];
    start_t1 = intersection.Ts1[0];
    is_start_sample = intersection.IsSamples[0];
    if (is_overlap) {
        end_t0 = intersection.Ts0[1];
        end_t1 = intersection.Ts1[1];
        is_end_sample = intersection.IsSamples[1];
    }
    else {
        end_t0 = intersection.Ts0[0];
        end_t1 = intersection.Ts1[0];
        is_end_sample = intersection.IsSamples[0];
    }
}
