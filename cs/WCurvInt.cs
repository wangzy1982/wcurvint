using System;
using System.Runtime.InteropServices;

namespace WCurvInt
{
    /// <summary>
    /// Represents a 2D point with X and Y coordinates
    /// </summary>
    public struct Point2d
    {
        /// <summary>
        /// X coordinate of the point
        /// </summary>
        public double X;

        /// <summary>
        /// Y coordinate of the point
        /// </summary>
        public double Y;

        /// <summary>
        /// Initializes a new instance of the Point2d struct
        /// </summary>
        /// <param name="x">X coordinate value</param>
        /// <param name="y">Y coordinate value</param>
        public Point2d(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// Returns a string representation of the Point2d
        /// </summary>
        /// <returns>Formatted string with X and Y values</returns>
        public override string ToString()
        {
            return $"({X:F3}, {Y:F3})";
        }
    }

    /// <summary>
    /// Abstract base class for all 2D curve types
    /// Manages native curve pointer and resource disposal
    /// </summary>
    public abstract class Curve2D : IDisposable
    {
        /// <summary>
        /// Native pointer to the curve object in C library
        /// </summary>
        protected IntPtr _curvePtr;

        /// <summary>
        /// Flag to track if object has been disposed
        /// </summary>
        private bool _disposed = false;

        /// <summary>
        /// Finalizer to release unmanaged resources
        /// </summary>
        ~Curve2D()
        {
            Dispose(false);
        }

        /// <summary>
        /// Releases all resources used by the Curve2D object
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Releases the unmanaged resources used by the Curve2D and optionally releases managed resources
        /// </summary>
        /// <param name="disposing">true to release both managed and unmanaged resources; false to release only unmanaged resources</param>
        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                // Release unmanaged resources
                if (_curvePtr != IntPtr.Zero)
                {
                    Interop.free_curve2d(_curvePtr);
                    _curvePtr = IntPtr.Zero;
                }

                _disposed = true;
            }
        }

        /// <summary>
        /// Calculates the point on the curve at the specified parameter t
        /// </summary>
        /// <param name="t">Parameter value along the curve (typically 0.0 to 1.0)</param>
        /// <returns>2D point on the curve at parameter t</returns>
        /// <exception cref="ObjectDisposedException">Thrown if curve is already disposed</exception>
        public Point2d CalculatePoint(double t)
        {
            CheckDisposed();
            Interop.calculate_curve_point(_curvePtr, t, out double x, out double y);
            return new Point2d(x, y);
        }

        /// <summary>
        /// Internal access to native curve pointer (for intersection calculations)
        /// </summary>
        internal IntPtr CurvePtr => _curvePtr;

        /// <summary>
        /// Verifies if the curve object is disposed and throws exception if true
        /// </summary>
        /// <exception cref="ObjectDisposedException">Thrown if object is disposed</exception>
        protected void CheckDisposed()
        {
            if (_disposed || _curvePtr == IntPtr.Zero)
            {
                throw new ObjectDisposedException(GetType().Name);
            }
        }
    }

    /// <summary>
    /// Represents a 2D line segment
    /// </summary>
    public class Line2D : Curve2D
    {
        /// <summary>
        /// Initializes a new instance of the Line2D class
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown if native line creation fails</exception>
        public Line2D()
        {
            _curvePtr = Interop.new_line2d();
            if (_curvePtr == IntPtr.Zero)
            {
                throw new InvalidOperationException("Failed to create native Line2D object");
            }
        }

        /// <summary>
        /// Sets the start point of the line segment
        /// </summary>
        /// <param name="startPoint">Start point of the line</param>
        /// <exception cref="ObjectDisposedException">Thrown if line is already disposed</exception>
        public void SetStartPoint(Point2d startPoint)
        {
            CheckDisposed();
            Interop.set_line2d_start_point(_curvePtr, startPoint.X, startPoint.Y);
        }

        /// <summary>
        /// Sets the end point of the line segment
        /// </summary>
        /// <param name="endPoint">End point of the line</param>
        /// <exception cref="ObjectDisposedException">Thrown if line is already disposed</exception>
        public void SetEndPoint(Point2d endPoint)
        {
            CheckDisposed();
            Interop.set_line2d_end_point(_curvePtr, endPoint.X, endPoint.Y);
        }
    }

    /// <summary>
    /// Represents a 2D circular arc
    /// </summary>
    public class Arc2D : Curve2D
    {
        /// <summary>
        /// Initializes a new instance of the Arc2D class
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown if native arc creation fails</exception>
        public Arc2D()
        {
            _curvePtr = Interop.new_arc2d();
            if (_curvePtr == IntPtr.Zero)
            {
                throw new InvalidOperationException("Failed to create native Arc2D object");
            }
        }

        /// <summary>
        /// Sets the center point of the arc
        /// </summary>
        /// <param name="centerPoint">Center point of the arc</param>
        /// <exception cref="ObjectDisposedException">Thrown if arc is already disposed</exception>
        public void SetCenter(Point2d centerPoint)
        {
            CheckDisposed();
            Interop.set_arc2d_center(_curvePtr, centerPoint.X, centerPoint.Y);
        }

        /// <summary>
        /// Sets the radius of the arc
        /// </summary>
        /// <param name="radius">Radius value (must be non-negative)</param>
        /// <exception cref="ArgumentException">Thrown if radius is negative</exception>
        /// <exception cref="ObjectDisposedException">Thrown if arc is already disposed</exception>
        public void SetRadius(double radius)
        {
            CheckDisposed();
            if (radius < 0)
            {
                throw new ArgumentException("Radius cannot be negative", nameof(radius));
            }
            Interop.set_arc2d_radius(_curvePtr, radius);
        }

        /// <summary>
        /// Sets the start angle of the arc (in radians)
        /// </summary>
        /// <param name="startAngle">Start angle in radians</param>
        /// <exception cref="ObjectDisposedException">Thrown if arc is already disposed</exception>
        public void SetStartAngle(double startAngle)
        {
            CheckDisposed();
            Interop.set_arc2d_start_angle(_curvePtr, startAngle);
        }

        /// <summary>
        /// Sets the delta angle (sweep angle) of the arc (in radians)
        /// </summary>
        /// <param name="deltaAngle">Delta angle in radians</param>
        /// <exception cref="ObjectDisposedException">Thrown if arc is already disposed</exception>
        public void SetDeltaAngle(double deltaAngle)
        {
            CheckDisposed();
            Interop.set_arc2d_delta_angle(_curvePtr, deltaAngle);
        }
    }

    /// <summary>
    /// Represents a 2D Bezier curve
    /// </summary>
    public class BezierCurve2D : Curve2D
    {
        private readonly int _degree;

        /// <summary>
        /// Initializes a new instance of the BezierCurve2D class
        /// </summary>
        /// <param name="degree">Degree of the Bezier curve (must be ≥ 1)</param>
        /// <exception cref="ArgumentException">Thrown if degree is less than 1</exception>
        /// <exception cref="InvalidOperationException">Thrown if native Bezier curve creation fails</exception>
        public BezierCurve2D(int degree)
        {
            if (degree < 1)
            {
                throw new ArgumentException("Curve degree must be greater than or equal to 1", nameof(degree));
            }

            _degree = degree;
            _curvePtr = Interop.new_bezier_curve2d(degree);

            if (_curvePtr == IntPtr.Zero)
            {
                throw new InvalidOperationException("Failed to create native BezierCurve2D object");
            }
        }

        /// <summary>
        /// Sets a control point of the Bezier curve
        /// </summary>
        /// <param name="index">Index of the control point (0 to degree)</param>
        /// <param name="controlPoint">Control point coordinates</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if index is out of valid range</exception>
        /// <exception cref="ObjectDisposedException">Thrown if curve is already disposed</exception>
        public void SetControlPoint(int index, Point2d controlPoint)
        {
            CheckDisposed();

            if (index < 0 || index > _degree)
            {
                throw new ArgumentOutOfRangeException(nameof(index), $"Index must be between 0 and {_degree}");
            }

            Interop.set_bezier_curve2d_control_point(_curvePtr, index, controlPoint.X, controlPoint.Y);
        }
    }

    /// <summary>
    /// Represents a 2D rational Bezier curve with weight values
    /// </summary>
    public class RationalBezierCurve2D : Curve2D
    {
        private readonly int _degree;

        /// <summary>
        /// Initializes a new instance of the RationalBezierCurve2D class
        /// </summary>
        /// <param name="degree">Degree of the rational Bezier curve (must be ≥ 1)</param>
        /// <exception cref="ArgumentException">Thrown if degree is less than 1</exception>
        /// <exception cref="InvalidOperationException">Thrown if native rational Bezier curve creation fails</exception>
        public RationalBezierCurve2D(int degree)
        {
            if (degree < 1)
            {
                throw new ArgumentException("Curve degree must be greater than or equal to 1", nameof(degree));
            }

            _degree = degree;
            _curvePtr = Interop.new_rational_bezier_curve2d(degree);

            if (_curvePtr == IntPtr.Zero)
            {
                throw new InvalidOperationException("Failed to create native RationalBezierCurve2D object");
            }
        }

        /// <summary>
        /// Sets a control point of the rational Bezier curve
        /// </summary>
        /// <param name="index">Index of the control point (0 to degree)</param>
        /// <param name="controlPoint">Control point coordinates</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if index is out of valid range</exception>
        /// <exception cref="ObjectDisposedException">Thrown if curve is already disposed</exception>
        public void SetControlPoint(int index, Point2d controlPoint)
        {
            CheckDisposed();

            if (index < 0 || index > _degree)
            {
                throw new ArgumentOutOfRangeException(nameof(index), $"Index must be between 0 and {_degree}");
            }

            Interop.set_rational_bezier_curve2d_control_point(_curvePtr, index, controlPoint.X, controlPoint.Y);
        }

        /// <summary>
        /// Sets the weight value for a control point of the rational Bezier curve
        /// </summary>
        /// <param name="index">Index of the control point (0 to degree)</param>
        /// <param name="weight">Weight value for the control point</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if index is out of valid range</exception>
        /// <exception cref="ObjectDisposedException">Thrown if curve is already disposed</exception>
        public void SetWeight(int index, double weight)
        {
            CheckDisposed();

            if (index < 0 || index > _degree)
            {
                throw new ArgumentOutOfRangeException(nameof(index), $"Index must be between 0 and {_degree}");
            }

            Interop.set_rational_bezier_curve2d_weight(_curvePtr, index, weight);
        }
    }

    /// <summary>
    /// Represents detailed information about a single curve intersection
    /// </summary>
    public struct IntersectionInfo
    {
        /// <summary>
        /// Indicates if the intersection is an overlap segment (true) or a single point (false)
        /// </summary>
        public bool IsOverlap;

        /// <summary>
        /// Start parameter on the first curve
        /// </summary>
        public double StartT0;

        /// <summary>
        /// Start parameter on the second curve
        /// </summary>
        public double StartT1;

        /// <summary>
        /// Indicates if the start point is a sample point
        /// </summary>
        public bool IsStartSample;

        /// <summary>
        /// End parameter on the first curve (for overlap segments)
        /// </summary>
        public double EndT0;

        /// <summary>
        /// End parameter on the second curve (for overlap segments)
        /// </summary>
        public double EndT1;

        /// <summary>
        /// Indicates if the end point is a sample point
        /// </summary>
        public bool IsEndSample;

        /// <summary>
        /// Calculates the start intersection point using the first curve
        /// </summary>
        /// <param name="curve0">First curve used in intersection calculation</param>
        /// <returns>Start point of the intersection</returns>
        public Point2d GetStartPoint(Curve2D curve0)
        {
            return curve0.CalculatePoint(StartT0);
        }

        /// <summary>
        /// Calculates the end intersection point using the first curve (for overlaps)
        /// </summary>
        /// <param name="curve0">First curve used in intersection calculation</param>
        /// <returns>End point of the intersection segment</returns>
        public Point2d GetEndPoint(Curve2D curve0)
        {
            return curve0.CalculatePoint(EndT0);
        }

        /// <summary>
        /// Returns a string representation of the IntersectionInfo
        /// </summary>
        /// <returns>Formatted string with intersection details</returns>
        public override string ToString()
        {
            if (IsOverlap)
            {
                return $"Overlap - Start (T0: {StartT0:F3}, T1: {StartT1:F3}), End (T0: {EndT0:F3}, T1: {EndT1:F3})";
            }
            return $"Point - T0: {StartT0:F3}, T1: {StartT1:F3}";
        }
    }

    /// <summary>
    /// Represents intersection results between two 2D curves
    /// Manages native intersection pointer and resource disposal
    /// </summary>
    public class CurveCurveIntersections : IDisposable
    {
        /// <summary>
        /// Native pointer to the intersection results in C library
        /// </summary>
        private IntPtr _intersectionsPtr;

        /// <summary>
        /// Flag to track if object has been disposed
        /// </summary>
        private bool _disposed = false;

        /// <summary>
        /// Initializes a new instance of the CurveCurveIntersections class
        /// </summary>
        /// <param name="curve0">First curve for intersection calculation</param>
        /// <param name="curve1">Second curve for intersection calculation</param>
        /// <param name="distanceEpsilon">Tolerance for distance calculations</param>
        /// <exception cref="ArgumentNullException">Thrown if either curve is null</exception>
        /// <exception cref="InvalidOperationException">Thrown if native intersection creation fails</exception>
        public CurveCurveIntersections(Curve2D curve0, Curve2D curve1, double distanceEpsilon)
        {
            if (curve0 == null)
            {
                throw new ArgumentNullException(nameof(curve0), "First curve cannot be null");
            }

            if (curve1 == null)
            {
                throw new ArgumentNullException(nameof(curve1), "Second curve cannot be null");
            }

            _intersectionsPtr = Interop.new_curve_curve_intersections(
                curve0.CurvePtr,
                curve1.CurvePtr,
                distanceEpsilon);

            if (_intersectionsPtr == IntPtr.Zero)
            {
                throw new InvalidOperationException("Failed to create native intersection object");
            }
        }

        /// <summary>
        /// Finalizer to release unmanaged resources
        /// </summary>
        ~CurveCurveIntersections()
        {
            Dispose(false);
        }

        /// <summary>
        /// Releases all resources used by the CurveCurveIntersections object
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Releases the unmanaged resources used by the CurveCurveIntersections and optionally releases managed resources
        /// </summary>
        /// <param name="disposing">true to release both managed and unmanaged resources; false to release only unmanaged resources</param>
        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                // Release unmanaged resources
                if (_intersectionsPtr != IntPtr.Zero)
                {
                    Interop.free_curve_curve_intersections(_intersectionsPtr);
                    _intersectionsPtr = IntPtr.Zero;
                }

                _disposed = true;
            }
        }

        /// <summary>
        /// Gets singularity information for the intersection results
        /// </summary>
        /// <param name="isSingularity0">Output flag indicating if first curve has singularity</param>
        /// <param name="isSingularity1">Output flag indicating if second curve has singularity</param>
        /// <exception cref="ObjectDisposedException">Thrown if intersection object is already disposed</exception>
        public void GetSingularityInfo(out bool isSingularity0, out bool isSingularity1)
        {
            CheckDisposed();
            Interop.is_curve_curve_intersections_singularity(
                _intersectionsPtr,
                out isSingularity0,
                out isSingularity1);
        }

        /// <summary>
        /// Gets the number of intersection points/segments
        /// </summary>
        /// <returns>Count of intersections</returns>
        /// <exception cref="ObjectDisposedException">Thrown if intersection object is already disposed</exception>
        public int GetIntersectionCount()
        {
            CheckDisposed();
            return Interop.get_curve_curve_intersection_count(_intersectionsPtr);
        }

        /// <summary>
        /// Gets detailed information about a specific intersection
        /// </summary>
        /// <param name="index">Index of the intersection (0 to count-1)</param>
        /// <returns>IntersectionInfo struct containing all intersection details</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if index is out of valid range</exception>
        /// <exception cref="ObjectDisposedException">Thrown if intersection object is already disposed</exception>
        public IntersectionInfo GetIntersection(int index)
        {
            CheckDisposed();

            int count = GetIntersectionCount();
            if (index < 0 || index >= count)
            {
                throw new ArgumentOutOfRangeException(nameof(index), $"Index must be between 0 and {count - 1}");
            }

            Interop.get_curve_curve_intersection(
                _intersectionsPtr,
                index,
                out bool isOverlap,
                out double startT0,
                out double startT1,
                out bool isStartSample,
                out double endT0,
                out double endT1,
                out bool isEndSample);

            return new IntersectionInfo
            {
                IsOverlap = isOverlap,
                StartT0 = startT0,
                StartT1 = startT1,
                IsStartSample = isStartSample,
                EndT0 = endT0,
                EndT1 = endT1,
                IsEndSample = isEndSample
            };
        }

        /// <summary>
        /// Verifies if the intersection object is disposed and throws exception if true
        /// </summary>
        /// <exception cref="ObjectDisposedException">Thrown if object is disposed</exception>
        private void CheckDisposed()
        {
            if (_disposed || _intersectionsPtr == IntPtr.Zero)
            {
                throw new ObjectDisposedException(nameof(CurveCurveIntersections));
            }
        }
    }

    /// <summary>
    /// Internal interop class for native C API calls
    /// Do not use directly in application code
    /// </summary>
    internal static class Interop
    {
        // Update with your actual DLL name (e.g., "WCURVINT.dll")
        private const string DllName = "WCURVINT.dll";

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr new_line2d();

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_line2d_start_point(IntPtr curve, double x, double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_line2d_end_point(IntPtr curve, double x, double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr new_arc2d();

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_arc2d_center(IntPtr curve, double x, double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_arc2d_radius(IntPtr curve, double radius);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_arc2d_start_angle(IntPtr curve, double start_angle);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_arc2d_delta_angle(IntPtr curve, double delta_angle);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr new_bezier_curve2d(int degree);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_bezier_curve2d_control_point(IntPtr curve, int index, double x, double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr new_rational_bezier_curve2d(int degree);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_rational_bezier_curve2d_control_point(IntPtr curve, int index, double x, double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void set_rational_bezier_curve2d_weight(IntPtr curve, int index, double w);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void calculate_curve_point(IntPtr curve, double t, out double x, out double y);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void free_curve2d(IntPtr curve);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr new_curve_curve_intersections(IntPtr curve0, IntPtr curve1, double distance_epsilon);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void free_curve_curve_intersections(IntPtr intersections);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void is_curve_curve_intersections_singularity(
            IntPtr intersections,
            [MarshalAs(UnmanagedType.U1)] out bool is_singularity0,
            [MarshalAs(UnmanagedType.U1)] out bool is_singularity1);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern int get_curve_curve_intersection_count(IntPtr intersections);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void get_curve_curve_intersection(
            IntPtr intersections,
            int index,
            [MarshalAs(UnmanagedType.U1)] out bool is_overlap,
            out double start_t0,
            out double start_t1,
            [MarshalAs(UnmanagedType.U1)] out bool is_start_sample,
            out double end_t0,
            out double end_t1,
            [MarshalAs(UnmanagedType.U1)] out bool is_end_sample);
    }
}