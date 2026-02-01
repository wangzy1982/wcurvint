# wcurvint: A Lightweight Open-Source 2D Curve Intersection C/C++ Library

wcurvint is a lightweight, easy-to-use open-source C/C++ library focused on solving the stability issues of intersection calculation for common 2D curves. It provides a concise C language API, making it ideal for personal projects, research validation, and rapid integration into small and medium-sized applications.

## 📋 Features & Capabilities

### Core Intersection Support

- The library supports pairwise intersection calculations for a wide range of basic 2D curves:
  - Line ↔ Line
  - Line ↔ Circular Arc
  - Line ↔ Bézier Curve
  - Line ↔ Rational Bézier Curve
  - Circular Arc ↔ Circular Arc
  - Circular Arc ↔ Bézier Curve
  - Circular Arc ↔ Rational Bézier Curve
  - Bézier Curve ↔ Bézier Curve
  - Bézier Curve ↔ Rational Bézier Curve
  - Rational Bézier Curve ↔ Rational Bézier Curve
- Rich Intersection Result Information:
  - Number of intersection results
  - Curve singularity check (identifies points with zero derivative vectors)
  - Exact intersection points and overlapping curve segments
  - Intersection type flag (cross point or sample point, for upper-layer application strategy selection)

### Key Advantages
- Robust Calculation Logic
  - Basic interval iteration to prevent missing solutions
  - Interval algebraic polynomial elimination to avoid precision degradation from equation simplification
  - Special endpoint processing for guaranteed topological stability
  - Optimized intersection results that align with practical engineering expectations
- Pure C language API for seamless integration with C/C++, Python, C# and other programming languages

## 💡 Commercial Edition Introduction

This open-source version meets general 2D curve intersection requirements.
For scenarios requiring higher computational efficiency or cross-platform support, we offer a commercial edition of the curve intersection library with exclusive enhancements:
- Performance: The commercial edition delivers <b>5x the efficiency</b> of the open-source version with advanced algorithm optimization.
- Cross-Platform Support: The open-source version is limited to the Windows platform, while the <b>commercial edition supports multiple operating systems and architectures</b>.

<b>For more details about the commercial edition, please contact us via email: 1179422870@qq.com</b>

## 🚀 Quick Start

### Environment Prerequisites

1. The open-source version currently only supports the Windows platform. For cross-platform support, refer to the Commercial Edition introduction above.
2. Clone or download the source code from this repository
3. Compile the project by running the build.sh script
4. Include the header file wcurvint.h in your project
5. Link against the compiled wcurvint.lib import library, and place the runtime libraries wcurvint.dll and wsolver.dll in the same directory as your executable file

### Code Examples

<b>Line & Circular Arc Intersection</b>
```cpp
void test_line_line() {
    // Create a 2D line curve
    void* curve0 = new_line2d();
    set_line2d_start_point(curve0, 0, 0);
    set_line2d_end_point(curve0, 100, 100);
    // Create a 2D circular arc curve
    void* curve1 = new_arc2d();
    set_arc2d_center(curve1, 50, 50);
    set_arc2d_radius(curve1, 10);
    set_arc2d_start_angle(curve1, 0);
    set_arc2d_delta_angle(curve1, 6.28);
    // Calculate curve intersections (tolerance: 1E-6)
    void* intersections = new_curve_curve_intersections(curve0, curve1, 1E-6);
    // Print intersection results
    int count = get_curve_curve_intersection_count(intersections);
    std::cout << "Intersection count: " << count << std::endl;
    for (int i = 0; i < count; ++i) {
        bool is_overlap;
        double start_t0;
        double start_t1;
        bool is_start_sample;
        double end_t0;
        double end_t1;
        bool is_end_sample;
        get_curve_curve_intersection(intersections, i, is_overlap,
            start_t0, start_t1, is_start_sample, end_t0, end_t1, is_end_sample);
        std::cout << "Intersection" << i << ": ";
        if (is_overlap) {
            std::cout << "[" << start_t0 << ", " << end_t0 << "]" << "[" << start_t1 << ", " << end_t1 << "]" << std::endl;
        }
        else {
            std::cout << "[" << start_t0 << "]" << "[" << start_t1 << "]" << std::endl;
        }
    }
    // Release allocated memory
    free(intersections);
    free_curve2d(curve0);
    free_curve2d(curve1);
}
```

<b>Bézier Curve & Bézier Curve Intersection</b>
```cpp
void test_bezier_bezier() {
    // Create a 3rd-order 2D Bézier curve
    void* curve0 = new_bezier_curve2d(3);
    set_bezier_curve2d_control_point(curve0, 0, 0, 0);
    set_bezier_curve2d_control_point(curve0, 1, 100, 100);
    set_bezier_curve2d_control_point(curve0, 2, 200, 100);
    set_bezier_curve2d_control_point(curve0, 3, 300, 0);
    // Create another 3rd-order 2D Bézier curve   
    void* curve1 = new_bezier_curve2d(3);
    set_bezier_curve2d_control_point(curve1, 0, 100, 0);
    set_bezier_curve2d_control_point(curve1, 1, 200, 100);
    set_bezier_curve2d_control_point(curve1, 2, 300, 100);
    set_bezier_curve2d_control_point(curve1, 3, 400, 0);
    // Calculate curve intersections (tolerance: 1E-6)
    void* intersections = new_curve_curve_intersections(curve0, curve1, 1E-6);
    // Print intersection results
    int count = get_curve_curve_intersection_count(intersections);
    std::cout << "Intersection count: " << count << std::endl;
    for (int i = 0; i < count; ++i) {
        bool is_overlap;
        double start_t0;
        double start_t1;
        bool is_start_sample;
        double end_t0;
        double end_t1;
        bool is_end_sample;
        get_curve_curve_intersection(intersections, i, is_overlap,
            start_t0, start_t1, is_start_sample, end_t0, end_t1, is_end_sample);
        std::cout << "Intersection" << i << ": ";
        if (is_overlap) {
            std::cout << "[" << start_t0 << ", " << end_t0 << "]" << "[" << start_t1 << ", " << end_t1 << "]" << std::endl;
        }
        else {
            std::cout << "[" << start_t0 << "]" << "[" << start_t1 << "]" << std::endl;
        }
    }
    // Release allocated memory
    free(intersections);
    free_curve2d(curve0);
    free_curve2d(curve1);
}
```

## 📞 Contact Us

- Email: 1179422870@qq.com