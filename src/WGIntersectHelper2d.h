#ifndef _WG_INTERSECT_HELPER_2D_
#define _WG_INTERSECT_HELPER_2D_

#include "wsolver.h"
#include "wsbasis.h"
#include "WGVector2d.h"
#include "WGUtil.h"

class WGIntersectHelper2d {
public:
    class PointCurveEquationSystem : public WSEquationSystem {
    public:
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef) = 0;
        virtual void CalculateD0(double t, WGVector2d& d0) const = 0;
        virtual void CalculateD1(double t, WGVector2d& d1) const = 0;
        virtual void CalculateD1(const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const = 0;
        virtual void CalculateD0(WSEquationsCache* cache, WSInterval& x, WSInterval& y) const = 0;
        virtual void CalculateD1(WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const = 0;
    };

    class CurveCurveEquationSystem : public WSEquationSystem {
    public:
        CurveCurveEquationSystem();
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable) = 0;
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef) = 0;
        virtual int GetSplitIndex(const WSIntervalVector* variable) = 0;
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const = 0;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const = 0;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const = 0;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const = 0;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const = 0;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const = 0;
        virtual int GetSampleCount() = 0;
    public:
        void SetAlgebraEnable(const WSIntervalVector* variable_domain, double algebra_epsilon);
    protected:
        void BuildFittingPolynomials(const WSIntervalVector* variable_domain, WSPolynomial**& polynomials, int& polynomial_count) const;
    protected:
        int m_algebra_enable;
        double m_algebra_epsilon;
    };

    class PointCurveSolver {
    public:
        PointCurveSolver();
        void Execute(PointCurveEquationSystem* equations, const WSIntervalVector* variable_domain,
            int max_clear_root_count, int max_fuzzy_domain_count);
    public:
        int GetClearRootCount();
        const WSIntervalVector* GetClearRoot(int index);
        int GetFuzzyDomainCount();
        const WSIntervalVector* GetFuzzyDomain(int index);
        int GetFuzzyTag(int index);
        bool GetIsSingularity();
    private:
        class HeapItem {
        public:
            HeapItem(WSCache* cache, int variable_cache_offset, int terms_value_cache_offset,
                int terms_partial_derivative_cache_offset, int terms_linear_a_cache_offset, int terms_linear_b_cache_offset,
                int terms_dirty_flag_cache_offset, int variable_count, int term_value_count,
                int term_partial_derivative_count, int term_linear_count) :
                m_variable(cache, variable_cache_offset, variable_count),
                m_terms_value(cache, terms_value_cache_offset, term_value_count),
                m_terms_partial_derivative(cache, terms_partial_derivative_cache_offset, term_partial_derivative_count, variable_count),
                m_terms_linear_a(cache, terms_linear_a_cache_offset, term_linear_count, variable_count),
                m_terms_linear_b(cache, terms_linear_b_cache_offset, term_linear_count),
                m_terms_dirty_flag((int*)cache->Data(terms_dirty_flag_cache_offset)),
                m_tag(0),
                m_priority0(0),
                m_priority1(0) {
            }

            bool LessThan(const HeapItem& other) {
                if (m_priority0 == other.m_priority0) {
                    return m_priority1 < other.m_priority1;
                }
                return m_priority0 < other.m_priority0;
            }
        public:
            WSSliceIntervalVector m_variable;
            WSSliceIntervalVector m_terms_value;
            WSSliceIntervalMatrix m_terms_partial_derivative;
            WSSliceMatrix m_terms_linear_a;
            WSSliceIntervalVector m_terms_linear_b;
            int* m_terms_dirty_flag;
            //m_tag: 
            //    bit0--approximately linear  
            //    bit1--very small  
            int m_tag;
            int m_priority0;
            WSReal m_priority1;
        };

        class Iterator : public WSIterator {
        public:
            Iterator(WSEquationSystem* equations, WSCache* cache, int cache_offset) :
                WSIterator(equations, cache, cache_offset) {
            }
        protected:
            virtual bool CheckTerminateEarly(WSEquationSystem* equations, WSEquationsCache* cache, int& tag);
        };
    private:
        void ReserveClearRoots(int clear_root_capacity);
        void ReserveHeap(int heap_capacity);
        WSIterateResult Iterate(WSIterator* iterator, HeapItem* heap_item, double min_rate,
            int max_iterate_count, int& terminate_early_tag);
        void CopyHeapItem(HeapItem* dst, HeapItem* src);
        void CalculateHeapItemTag(HeapItem* heap_item, int prev_tag, double limit_angle);
        void CalculatePriority(HeapItem* heap_item);
    private:
        PointCurveEquationSystem* m_equations;
        int m_max_clear_root_count;
        int m_max_fuzzy_domain_count;
        int m_clear_root_capacity;
        int m_clear_root_object_offset;
        int m_clear_root_data_offset;
        int m_heap_capacity;
        int m_heap_object_offset;
        int m_heap_data_offset;
        HeapItem* m_heap;
        int m_heap_item_count;
        WSSliceIntervalVector* m_clear_roots;
        int m_clear_root_count;
        HeapItem* m_iterate_heap_item;
        bool m_is_singularity;
    private:
        WSCache m_cache;
        WSCache m_iterator_cache;
    public:
        PointCurveSolver* Next;
    };

    class CurveCurveSolver {
    public:
        CurveCurveSolver();
    public:
        void PrepareExecute(CurveCurveEquationSystem* equations, const WSIntervalVector* variable_domain, 
            int max_clear_root_count, int max_fuzzy_domain_count);
        void PrimaryExecute();
        WSIterateResult FinallyIterate(int heap_index, double distance_epsilon);
        WSIterateResult Iterate(WSIntervalVector* variable, double check_coef, double iterate_coef);
    public:
        int GetClearRootCount();
        const WSIntervalVector* GetClearRoot(int index);
        int GetFuzzyDomainCount();
        const WSIntervalVector* GetFuzzyDomain(int index);
        int GetFuzzyTag(int index);
        bool GetIsSingularity0();
        bool GetIsSingularity1();
    private:
        class HeapItem {
        public:
            HeapItem(WSCache* cache, int variable_cache_offset, int terms_value_cache_offset,
                int terms_partial_derivative_cache_offset, int terms_linear_a_cache_offset, int terms_linear_b_cache_offset,
                int terms_dirty_flag_cache_offset, int variable_count, int term_value_count,
                int term_partial_derivative_count, int term_linear_count) :
                m_variable(cache, variable_cache_offset, variable_count),
                m_terms_value(cache, terms_value_cache_offset, term_value_count),
                m_terms_partial_derivative(cache, terms_partial_derivative_cache_offset, term_partial_derivative_count, variable_count),
                m_terms_linear_a(cache, terms_linear_a_cache_offset, term_linear_count, variable_count),
                m_terms_linear_b(cache, terms_linear_b_cache_offset, term_linear_count),
                m_terms_dirty_flag((int*)cache->Data(terms_dirty_flag_cache_offset)),
                m_tag(0),
                m_priority0(0),
                m_priority1(0) {
            }

            bool LessThan(const HeapItem& other) {
                if (m_priority0 == other.m_priority0) {
                    return m_priority1 < other.m_priority1;
                }
                return m_priority0 < other.m_priority0;
            }
        public:
            WSSliceIntervalVector m_variable;
            WSSliceIntervalVector m_terms_value;
            WSSliceIntervalMatrix m_terms_partial_derivative;
            WSSliceMatrix m_terms_linear_a;
            WSSliceIntervalVector m_terms_linear_b;
            int* m_terms_dirty_flag;
            //m_tag: 
            //    bit0--monotonic  
            //    bit1--approximately linear 0  
            //    bit2--approximately linear 1  
            //    bit3--very small 0  
            //    bit4--very small 1
            int m_tag;
            int m_priority0;
            WSReal m_priority1;
        };

        class PrimaryIterator : public WSIterator {
        public:
            PrimaryIterator(WSEquationSystem* equations, WSCache* cache, int cache_offset) :
                WSIterator(equations, cache, cache_offset) {
            }
        protected:
            virtual bool CheckTerminateEarly(WSEquationSystem* equations, WSEquationsCache* cache, int& tag);
        };
    private:
        void ReserveClearRoots(int clear_root_capacity);
        void ReserveHeap(int heap_capacity);
        WSIterateResult PrimaryIterate(WSIterator* iterator, HeapItem* heap_item, double min_rate, 
            int max_iterate_count, int& terminate_early_tag);
        WSIterateResult FinallyIterate(WSIterator* iterator, HeapItem* heap_item, double min_rate, int max_iterate_count);
        void CopyHeapItem(HeapItem* dst, HeapItem* src);
        void CalculateHeapItemTag(HeapItem* heap_item, int prev_tag, double limit_angle);
        void CalculatePriority(HeapItem* heap_item);
    private:
        CurveCurveEquationSystem* m_equations;
        int m_max_clear_root_count;
        int m_max_fuzzy_domain_count;
        int m_clear_root_capacity;
        int m_clear_root_object_offset;
        int m_clear_root_data_offset;
        int m_heap_capacity;
        int m_heap_object_offset;
        int m_heap_data_offset;
        HeapItem* m_heap;
        int m_heap_item_count;
        WSSliceIntervalVector* m_clear_roots;
        int m_clear_root_count;
        HeapItem* m_iterate_heap_item;
        bool m_is_singularity0;
        bool m_is_singularity1;
    private:
        WSCache m_cache;
        WSCache m_iterator_cache;
    public:
        CurveCurveSolver* Next;
    };

    struct CurveCurveSamplePoint {
        double T0;
        double T1;
        double Distance;
    };

    struct CurveCurveFuzzyPair {
        CurveCurveFuzzyPair* Next;
        double Length;
        WSInterval Domain0;
        WSInterval Domain1;
        WGVector2d SampleDirection;
        bool SameDirection;
        CurveCurveSamplePoint SamplePoints[2];
    };

    struct CurveCurveFuzzyPairList {
        CurveCurveFuzzyPairList* Next;
        WSInterval Domain0;
        WSInterval Domain1;
        bool SameDirection;
        CurveCurveFuzzyPair* TailPair;
    };

    struct CurveCurveIntersection {
        double Ts0[2];
        double Ts1[2];
        bool IsSamples[2];
        int PointCount;
    public:
        CurveCurveIntersection() {}
        CurveCurveIntersection(double t0, double t1, bool is_sample) {
            Ts0[0] = t0;
            Ts1[0] = t1;
            IsSamples[0] = is_sample;
            PointCount = 1;
        }
        CurveCurveIntersection(double t00, double t01, bool is_sample0, double t10, double t11, bool is_sample1) {
            Ts0[0] = t00;
            Ts1[0] = t01;
            IsSamples[0] = is_sample0;
            Ts0[1] = t10;
            Ts1[1] = t11;
            IsSamples[1] = is_sample1;
            PointCount = 2;
        }
    };

    struct PointCurveIntersection {
        double T;
    public:
        PointCurveIntersection() {}
        PointCurveIntersection(double t) {
            T = t;
        }
    };

    class PointBezierCurveEquationSystem : public PointCurveEquationSystem {
    public:
        PointBezierCurveEquationSystem();
        virtual ~PointBezierCurveEquationSystem();
        void Update(const WGVector2d& point, int degree, const WGVector2d* control_points, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual void CalculateD0(double t, WGVector2d& d0) const;
        virtual void CalculateD1(double t, WGVector2d& d1) const;
        virtual void CalculateD1(const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[3];
        WSEquation* m_equations[2];
        WSTerm* m_term_x[2];
        WSTerm* m_term_y[2];
    public:
        PointBezierCurveEquationSystem* Next;
    };

    class PointRationalBezierCurveEquationSystem : public PointCurveEquationSystem {
    public:
        PointRationalBezierCurveEquationSystem();
        virtual ~PointRationalBezierCurveEquationSystem();
        void Update(const WGVector2d& point, int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual void CalculateD0(double t, WGVector2d& d0) const;
        virtual void CalculateD1(double t, WGVector2d& d1) const;
        virtual void CalculateD1(const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[3];
        WSEquation* m_equations[2];
        WSTerm* m_term_x[2];
        WSTerm* m_term_y[2];
    public:
        PointRationalBezierCurveEquationSystem* Next;
    };

    class BezierCurveBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        BezierCurveBezierCurveEquationSystem();
        virtual ~BezierCurveBezierCurveEquationSystem();
        void Update(int degree0, const WGVector2d* control_points0,
            int degree1, const WGVector2d* control_points1, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        void UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
            int degree1, const WGVector2d* control_points1, double distance_epsilon);
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        double m_variable_epsilons[2];
        WSEquationBasis* m_basises[4];
        WSEquation* m_equations[2];
        WSTerm* m_term_x[2];
        WSTerm* m_term_y[2];
    public:
        BezierCurveBezierCurveEquationSystem* Next;
    };

    class BezierCurveRationalBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        BezierCurveRationalBezierCurveEquationSystem();
        virtual ~BezierCurveRationalBezierCurveEquationSystem();
        void Update(int degree0, const WGVector2d* control_points0,
            int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        void UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
            int degree1, const WGVector2d* control_points1, double distance_epsilon);
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        double m_variable_epsilons[2];
        WSEquationBasis* m_basises[10];
        WSEquation* m_equations[5];
        WSTerm* m_term_x[2];
        WSTerm* m_term_y[2];
        WSTerm* m_term_w1;
    public:
        BezierCurveRationalBezierCurveEquationSystem* Next;
    };

    class RationalBezierCurveRationalBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        RationalBezierCurveRationalBezierCurveEquationSystem();
        virtual ~RationalBezierCurveRationalBezierCurveEquationSystem();
        void Update(int degree0, const WGVector2d* control_points0, const double* weights0,
            int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        void UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
            int degree1, const WGVector2d* control_points1, double distance_epsilon);
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        double m_variable_epsilons[2];
        WSEquationBasis* m_basises[16];
        WSEquation* m_equations[8];
        WSTerm* m_term_x[2];
        WSTerm* m_term_y[2];
        WSTerm* m_term_w[2];
    public:
        RationalBezierCurveRationalBezierCurveEquationSystem* Next;
    };

    class LineBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        LineBezierCurveEquationSystem();
        virtual ~LineBezierCurveEquationSystem();
        void Update(const WGVector2d& start_point, const WGVector2d& end_point,
            int degree, const WGVector2d* control_points, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[4];
        WSEquation* m_equations[2];
        WGVector2d m_start_point;
        WGVector2d m_dir;
        WSTerm* m_term_bezier_x;
        WSTerm* m_term_bezier_y;
    public:
        LineBezierCurveEquationSystem* Next;
    };

    class LineRationalBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        LineRationalBezierCurveEquationSystem();
        virtual ~LineRationalBezierCurveEquationSystem();
        void Update(const WGVector2d& start_point, const WGVector2d& end_point,
            int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[5];
        WSEquation* m_equations[3];
        WGVector2d m_start_point;
        WGVector2d m_dir;
        WSTerm* m_term_bezier_x;
        WSTerm* m_term_bezier_y;
        WSTerm* m_term_bezier_w;
    public:
        LineRationalBezierCurveEquationSystem* Next;
    };

    class ArcBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        ArcBezierCurveEquationSystem();
        virtual ~ArcBezierCurveEquationSystem();
        void Update(const WGVector2d& center, double radius, double start_angle, double delta_angle,
            int degree, const WGVector2d* control_points, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[7];
        WSEquation* m_equations[3];
        WGVector2d m_center;
        double m_radius;
        double m_start_angle;
        double m_delta_angle;
        WSTerm* m_term_bezier_x;
        WSTerm* m_term_bezier_y;
    public:
        ArcBezierCurveEquationSystem* Next;
    };

    class ArcRationalBezierCurveEquationSystem : public CurveCurveEquationSystem {
    public:
        ArcRationalBezierCurveEquationSystem();
        virtual ~ArcRationalBezierCurveEquationSystem();
        void Update(const WGVector2d& center, double radius, double start_angle, double delta_angle,
            int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon);
    public:
        virtual int GetVariableCount() const;
        virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const;
        virtual int GetBasisCount() const;
        virtual WSEquationBasis* GetBasis(int index) const;
        virtual int GetEquationCount() const;
        virtual WSEquation* GetEquation(int index) const;
        virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const;
        virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const;
        virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
            WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const;
    public:
        virtual void InitializeIntermediateVariables(WSIntervalVector* variable);
        virtual void SetEpsilonCoef(double check_coef, double iterate_coef);
        virtual int GetSplitIndex(const WSIntervalVector* variable);
        virtual void CalculateD0(int curve_index, double t, WGVector2d& d0) const;
        virtual void CalculateD1(int curve_index, double t, WGVector2d& d1) const;
        virtual void CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const;
        virtual void CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const;
        virtual void CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const;
        virtual int GetSampleCount();
    private:
        bool m_initialized;
        double m_check_coef;
        double m_iterate_coef;
        double m_distance_epsilon;
        WSEquationBasis* m_basises[13];
        WSEquation* m_equations[6];
        WGVector2d m_center;
        double m_radius;
        double m_start_angle;
        double m_delta_angle; 
        WSTerm* m_term_bezier_x;
        WSTerm* m_term_bezier_y;
        WSTerm* m_term_bezier_w;
    public:
        ArcRationalBezierCurveEquationSystem* Next;
    };

    class IntersectCache {
    public:
        IntersectCache();
        virtual ~IntersectCache();
        PointBezierCurveEquationSystem* AllocPointBezierCurveEquations();
        void FreePointBezierCurveEquations(PointBezierCurveEquationSystem* equations);
        PointRationalBezierCurveEquationSystem* AllocPointRationalBezierCurveEquations();
        void FreePointRationalBezierCurveEquations(PointRationalBezierCurveEquationSystem* equations);
        LineBezierCurveEquationSystem* AllocLineBezierCurveEquations();
        void FreeLineBezierCurveEquations(LineBezierCurveEquationSystem* equations);
        LineRationalBezierCurveEquationSystem* AllocLineRationalBezierCurveEquations();
        void FreeLineRationalBezierCurveEquations(LineRationalBezierCurveEquationSystem* equations);
        ArcBezierCurveEquationSystem* AllocArcBezierCurveEquations();
        void FreeArcBezierCurveEquations(ArcBezierCurveEquationSystem* equations);
        ArcRationalBezierCurveEquationSystem* AllocArcRationalBezierCurveEquations();
        void FreeArcRationalBezierCurveEquations(ArcRationalBezierCurveEquationSystem* equations);
        BezierCurveBezierCurveEquationSystem* AllocBezierCurveBezierCurveEquations();
        void FreeBezierCurveBezierCurveEquations(BezierCurveBezierCurveEquationSystem* equations);
        BezierCurveRationalBezierCurveEquationSystem* AllocBezierCurveRationalBezierCurveEquations();
        void FreeBezierCurveRationalBezierCurveEquations(BezierCurveRationalBezierCurveEquationSystem* equations);
        RationalBezierCurveRationalBezierCurveEquationSystem* AllocRationalBezierCurveRationalBezierCurveEquations();
        void FreeRationalBezierCurveRationalBezierCurveEquations(RationalBezierCurveRationalBezierCurveEquationSystem* equations);
        PointCurveSolver* AllocPointCurveSolver();
        void FreePointCurveSolver(PointCurveSolver* solver);
        CurveCurveSolver* AllocCurveCurveSolver();
        void FreeCurveCurveSolver(CurveCurveSolver* solver);
        CurveCurveFuzzyPairList* AllocPairList();
        void FreePairList(CurveCurveFuzzyPairList* fuzzy_pair_list);
        CurveCurveFuzzyPair* AllocPair();
        void FreePair(CurveCurveFuzzyPair* fuzzy_pair);
    private:
        PointBezierCurveEquationSystem* m_first_point_bezier_curve_equations;
        PointRationalBezierCurveEquationSystem* m_first_point_rational_bezier_curve_equations;
        LineBezierCurveEquationSystem* m_first_line_bezier_curve_equations;
        LineRationalBezierCurveEquationSystem* m_first_line_rational_bezier_curve_equations;
        ArcBezierCurveEquationSystem* m_first_arc_bezier_curve_equations;
        ArcRationalBezierCurveEquationSystem* m_first_arc_rational_bezier_curve_equations;
        BezierCurveBezierCurveEquationSystem* m_first_bezier_curve_bezier_curve_equations;
        BezierCurveRationalBezierCurveEquationSystem* m_first_bezier_curve_rational_bezier_curve_equations;
        RationalBezierCurveRationalBezierCurveEquationSystem* m_first_rational_bezier_curve_rational_bezier_curve_equations;
        PointCurveSolver* m_point_curve_solver;
        CurveCurveSolver* m_curve_curve_solver;
        CurveCurveFuzzyPair* m_first_pair;
        CurveCurveFuzzyPairList* m_first_pair_list;
    };
public:
    static void PreparePointCurveIntersect(PointCurveEquationSystem* equations, PointCurveSolver* solver,
        const WSIntervalVector* variable_domain, double distance_epsilon);
    static void BuildPointCurveIntersections(PointCurveEquationSystem* equations, PointCurveSolver* solver, 
        const WGVector2d& point, double distance_epsilon, std::vector<PointCurveIntersection>& intersections);
public:
    static void PrepareCurveCurveIntersect(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
        const WSIntervalVector* variable_domain, double distance_epsilon);
    static void BuildCurveCurveIntersections(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
        WSIntervalVector* runtime_variable, double distance_epsilon, IntersectCache* intersect_cache, 
        std::vector<CurveCurveIntersection>& intersections);
    static void MergeNotSampleIntersections(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable, 
        std::vector<CurveCurveIntersection>& intersections);
    static void AdjustIntersectionsWithEndpoints(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable,
        const std::vector<PointCurveIntersection>& intersections00, const std::vector<PointCurveIntersection>& intersections01,
        const std::vector<PointCurveIntersection>& intersections10, const std::vector<PointCurveIntersection>& intersections11,
        double distance_epsilon, std::vector<CurveCurveIntersection>& intersections);
public:
    static void PointBezierCurveIntersect(const WGVector2d& point, int degree, const WGVector2d* control_points, 
        double distance_epsilon, IntersectCache* intersect_cache, bool& is_singularity, std::vector<PointCurveIntersection>& intersections);
    static void PointRationalBezierCurveIntersect(const WGVector2d& point, int degree, const WGVector2d* control_points, const double* weights,
        double distance_epsilon, IntersectCache* intersect_cache, bool& is_singularity, std::vector<PointCurveIntersection>& intersections);
    static int PointLineIntersect(const WGVector2d& point, const WGVector2d& start_point, const WGVector2d& end_point,
        double distance_epsilon, bool& is_singularity, PointCurveIntersection intersections[1]);
    static int PointArcIntersect(const WGVector2d& point, const WGVector2d& center, double radius, double start_angle, double delta_angle, 
        double distance_epsilon, bool& is_singularity, PointCurveIntersection intersections[2]);
public:
    static void LineBezierCurveIntersect(const WGVector2d& start_point, const WGVector2d& end_point,
        int degree, const WGVector2d* control_points, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void LineRationalBezierCurveIntersect(const WGVector2d& start_point, const WGVector2d& end_point,
        int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void ArcBezierCurveIntersect(const WGVector2d& center, double radius, double start_angle, double delta_angle,
        int degree, const WGVector2d* control_points, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void ArcRationalBezierCurveIntersect(const WGVector2d& center, double radius, double start_angle, double delta_angle,
        int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void BezierCurveBezierCurveIntersect(int degree0, const WGVector2d* control_points0,
        int degree1, const WGVector2d* control_points1, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void BezierCurveRationalBezierCurveIntersect(int degree0, const WGVector2d* control_points0,
        int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static void RationalBezierCurveRationalBezierCurveIntersect(int degree0, const WGVector2d* control_points0, const double* weights0,
        int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon, IntersectCache* intersect_cache,
        bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections);
    static int LineLineIntersect(const WGVector2d& start_point0, const WGVector2d& end_point0,
        const WGVector2d& start_point1, const WGVector2d& end_point1, double distance_epsilon, 
        bool& is_singularity0, bool& is_singularity1, CurveCurveIntersection intersections[2]);
    static int LineArcIntersect(const WGVector2d& line_start_point, const WGVector2d& line_end_point, 
        const WGVector2d& arc_center, double arc_radius, double arc_start_angle, double arc_delta_angle, double distance_epsilon,
        bool& is_line_singularity, bool& is_arc_singularity, CurveCurveIntersection intersections[10]);
    static int ArcArcIntersect(const WGVector2d& center0, double radius0, double start_angle0, double delta_angle0,
        const WGVector2d& center1, double radius1, double start_angle1, double delta_angle1, double distance_epsilon,
        bool& is_singularity0, bool& is_singularity1, CurveCurveIntersection intersections[10]);
public:
    static int BuildArcArcDisjointExtremes(const WGVector2d& center_direction, double center_distance, 
        double radius0, double start_angle0, double delta_angle0, double radius1, double start_angle1, double delta_angle1,
        double distance_epsilon, double extreme_ts0[4], double extreme_ts1[4]);
public:
    static bool Sample(PointCurveEquationSystem* equations, const WGVector2d& point, const WGVector2d& sample_direction, 
        const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance);
    static void BuildFuzzyIntersection(PointCurveEquationSystem* equations, const WGVector2d& point, 
        const WSIntervalVector* variable, double distance_epsilon, std::vector<PointCurveIntersection>& intersections);
public:
    static WSInterval CalculateBezierD1(int degree, const double* coefs, const WSInterval& t);
    static WSInterval CalculateBezierD2(int degree, const double* coefs, const WSInterval& t);
    static WSInterval CalculateBezierD3(int degree, const double* coefs, const WSInterval& t);
    static WSInterval CalculateBezierD4(int degree, const double* coefs, const WSInterval& t);
public:
    static bool CheckMonotonic(CurveCurveEquationSystem* equations, const WSIntervalVector* variable);
    static bool CheckMonotonicOrShortLinear(CurveCurveEquationSystem* equations,
        const WSIntervalVector* variable, double limit_angle, double distance_epsilon);
    static void GetFuzzyDomainInfo(CurveCurveEquationSystem* equations, const WSIntervalVector* variable, 
        WGVector2d& sample_direction, bool& same_direction);
    static bool Sample(CurveCurveEquationSystem* equations, const WGVector2d& base_point, int sample_curve_index,
        const WGVector2d& sample_direction, const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance);
    static bool Sample(CurveCurveEquationSystem* equations, int base_curve_index, double base_t, 
        const WGVector2d& sample_direction, const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance);
    static bool MonotonicIterate(CurveCurveEquationSystem* equations, double distance_epsilon,
        const WSInterval& domain0, const WSInterval& domain1, double& t0, double& t1);
    static void BuildMonotonicIntersection(CurveCurveEquationSystem* equations, double distance_epsilon,
        const WSInterval& domain0, const WSInterval& domain1, std::vector<CurveCurveIntersection>& intersections);
    static CurveCurveFuzzyPair* BuildFuzzyPair(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable, 
        double distance_epsilon, const WSInterval& domain0, const WSInterval& domain1, IntersectCache* intersect_cache);
    static CurveCurveFuzzyPair* SplitFuzzyPair(CurveCurveFuzzyPair* fuzzy_pair, 
        double split_t0, double split_t1, double distance, IntersectCache* intersect_cache);
    static CurveCurveFuzzyPairList* NewFuzzyPairList(IntersectCache* intersect_cache, CurveCurveFuzzyPair* fuzzy_pair);
    static bool MergeFuzzyPairList(CurveCurveFuzzyPairList* fuzzy_pair_list, CurveCurveFuzzyPair* fuzzy_pair);
    static bool MergeFuzzyPairList(CurveCurveFuzzyPairList* dst_fuzzy_pair_list, CurveCurveFuzzyPairList* src_fuzzy_pair_list);
    static void BuildFuzzyIntersection(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
        WSIntervalVector* runtime_variable, double distance_epsilon, CurveCurveFuzzyPair* first_fuzzy_pair, 
        IntersectCache* intersect_cache, std::vector<CurveCurveIntersection>& intersections);
    static void SortTs(double* ts0, double* ts1, int* indices, int tc);
    static void RemoveSameTs(double* ts0, double* ts1, bool* is_samples, int& tc);
    static void RepairCross(double* ts0, double* ts1, bool* is_samples, WGVector2d* points0, WGVector2d* points1, int& tc, double distance_epsilon);
    static void MarkNotSample(double* ts0, double* ts1, bool* is_samples, int tc, CurveCurveIntersection* intersections, int intersection_count);
    static int MergeOverlapIntersections(CurveCurveIntersection* intersections, int intersection_count);
    static void RepairByEndTs(double* end_ts0, double* end_ts1, int end_tc, std::vector<CurveCurveIntersection>& intersections);
    static void RepairByEndTs(double* end_ts0, double* end_ts1, int end_tc, CurveCurveIntersection* intersections, int& intersection_count);
};

#endif
