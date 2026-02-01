#include "WGIntersectHelper2d.h"
#include "WGBasicAlgorithm.h"
#include "WGArc2d.h"
#include <assert.h>

const double g_curve_primary_linear_limit_angle = g_pi * 0.05;
const double g_curve_linear_limit_angle = g_pi * 0.001;

WGIntersectHelper2d::CurveCurveEquationSystem::CurveCurveEquationSystem() : 
    m_algebra_enable(0) {
}

void WGIntersectHelper2d::CurveCurveEquationSystem::SetAlgebraEnable(const WSIntervalVector* variable_domain, double algebra_epsilon) {
    if (m_algebra_enable == 0) {
        m_algebra_enable = 1;
        m_algebra_epsilon = algebra_epsilon;
        RebuildAlgebraRuntime(variable_domain);
        if (GetAlgebraPolynomialCount() == 0) {
            m_algebra_enable = 2;
        }
    }
    else if (m_algebra_enable == 1) {
        if (m_algebra_epsilon != algebra_epsilon) {
            m_algebra_epsilon = algebra_epsilon;
            RebuildAlgebraRuntime(variable_domain);
            if (GetAlgebraPolynomialCount() == 0) {
                m_algebra_enable = 2;
            }
        }
    }
    if (m_algebra_enable == 2) {
        m_algebra_enable = 3;
        m_algebra_epsilon = algebra_epsilon;
        RebuildAlgebraRuntime(variable_domain);
    }
    else if (m_algebra_enable == 3) {
        if (m_algebra_epsilon != algebra_epsilon) {
            m_algebra_epsilon = algebra_epsilon;
            RebuildAlgebraRuntime(variable_domain);
        }
    }
}

void WGIntersectHelper2d::CurveCurveEquationSystem::BuildFittingPolynomials(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count) const {
    const WSInterval& t0 = variable_domain->Get(0);
    const WSInterval& t1 = variable_domain->Get(1);
    double h0 = t0.Length();
    double h1 = t1.Length();
    double m0 = t0.Min;
    double m1 = t1.Min;
    WGVector2d p00, p01, p10, p11;
    CalculateD0(0, t0.Min, p00);
    CalculateD0(0, t0.Max, p01);
    CalculateD0(1, t1.Min, p10);
    CalculateD0(1, t1.Max, p11);
    WGVector2d d00, d01, d10, d11;
    CalculateD1(0, t0.Min, d00);
    CalculateD1(0, t0.Max, d01);
    CalculateD1(1, t1.Min, d10);
    CalculateD1(1, t1.Max, d11);
    double c0 = p00.X;
    double c1 = d00.X;
    double c2 = (3 * (p01.X - p00.X) - h0 * (2 * d00.X + d01.X)) / (h0 * h0);
    double c3 = (h0 * (d00.X + d01.X) - 2 * (p01.X - p00.X)) / (h0 * h0 * h0);
    double a03 = c3;
    double a02 = c2 - 3 * c3 * m0;
    double a01 = c1 - 2 * c2 * m0 + 3 * c3 * m0 * m0;
    double a00 = c0 - c1 * m0 + c2 * m0 * m0 - c3 * m0 * m0 * m0;
    c0 = p10.X;
    c1 = d10.X;
    c2 = (3 * (p11.X - p10.X) - h1 * (2 * d10.X + d11.X)) / (h1 * h1);
    c3 = (h1 * (d10.X + d11.X) - 2 * (p11.X - p10.X)) / (h1 * h1 * h1);
    double a13 = c3;
    double a12 = c2 - 3 * c3 * m1;
    double a11 = c1 - 2 * c2 * m1 + 3 * c3 * m1 * m1;
    double a10 = c0 - c1 * m1 + c2 * m1 * m1 - c3 * m1 * m1 * m1;
    c0 = p00.Y;
    c1 = d00.Y;
    c2 = (3 * (p01.Y - p00.Y) - h0 * (2 * d00.Y + d01.Y)) / (h0 * h0);
    c3 = (h0 * (d00.Y + d01.Y) - 2 * (p01.Y - p00.Y)) / (h0 * h0 * h0);
    double b03 = c3;
    double b02 = c2 - 3 * c3 * m0;
    double b01 = c1 - 2 * c2 * m0 + 3 * c3 * m0 * m0;
    double b00 = c0 - c1 * m0 + c2 * m0 * m0 - c3 * m0 * m0 * m0;
    c0 = p10.Y;
    c1 = d10.Y;
    c2 = (3 * (p11.Y - p10.Y) - h1 * (2 * d10.Y + d11.Y)) / (h1 * h1);
    c3 = (h1 * (d10.Y + d11.Y) - 2 * (p11.Y - p10.Y)) / (h1 * h1 * h1);
    double b13 = c3;
    double b12 = c2 - 3 * c3 * m1;
    double b11 = c1 - 2 * c2 * m1 + 3 * c3 * m1 * m1;
    double b10 = c0 - c1 * m1 + c2 * m1 * m1 - c3 * m1 * m1 * m1;

    WSInterval dx40, dy40;
    CalculateD4(0, variable_domain, dx40, dy40);
    double h40 = pow(variable_domain->Get(0).Length(), 4);
    WSInterval rx0 = dx40 / 384 * h40;
    WSInterval ry0 = dy40 / 384 * h40;
    WSInterval dx41, dy41;
    CalculateD4(1, variable_domain, dx41, dy41);
    double h41 = pow(variable_domain->Get(1).Length(), 4);
    WSInterval rx1 = dx41 / 384 * h41;
    WSInterval ry1 = dy41 / 384 * h41;
    
    int variable_count = GetVariableCount();
    polynomials = new WSPolynomial * [2];
    polynomial_count = 2; WSPolynomial* polynomial = new WSPolynomial(variable_count, 8);
    polynomials[0] = polynomial;
    WSPolynomialTerm* term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = a03;
    term->Powers[0] = 3;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = a02;
    term->Powers[0] = 2;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = a01;
    term->Powers[0] = 1;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -a13;
    term->Powers[1] = 3;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -a12;
    term->Powers[1] = 2;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -a11;
    term->Powers[1] = 1;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = a00 - a10 + rx0 - rx1 + WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
    polynomial->AddLast();

    polynomial = new WSPolynomial(variable_count, 8);
    polynomials[1] = polynomial;
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = b03;
    term->Powers[0] = 3;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = b02;
    term->Powers[0] = 2;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = b01;
    term->Powers[0] = 1;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -b13;
    term->Powers[1] = 3;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -b12;
    term->Powers[1] = 2;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = -b11;
    term->Powers[1] = 1;
    polynomial->AddLast();
    term = polynomial->GetTerm(polynomial->NewTerm());
    memset(term->Powers, 0, variable_count * sizeof(int));
    term->Coef = b00 - b10 + ry0 - ry1 + WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
    polynomial->AddLast();
}

WGIntersectHelper2d::PointCurveSolver::PointCurveSolver() :
    m_equations(nullptr),
    m_is_singularity(false),
    m_heap_item_count(0),
    m_clear_root_count(0) {
}

void WGIntersectHelper2d::PointCurveSolver::Execute(PointCurveEquationSystem* equations, const WSIntervalVector* variable_domain,
    int max_clear_root_count, int max_fuzzy_domain_count) {
    m_equations = equations;
    m_max_clear_root_count = max_clear_root_count;
    m_max_fuzzy_domain_count = max_fuzzy_domain_count;
    m_clear_root_capacity = 0;
    m_heap_capacity = 0;
    m_clear_root_count = 0;
    m_heap_item_count = 0;
    m_is_singularity = false;
    assert(max_clear_root_count > 0 && max_fuzzy_domain_count > 0);
    int variable_count = equations->GetVariableCount();
    int term_value_count = equations->GetTermValueCacheCount();
    int term_partial_derivative_count = equations->GetTermPartialDerivativeCacheCount();
    int term_linear_count = equations->GetTermLinearCacheCount();
    int heap_size = max_fuzzy_domain_count + 1;
    int cache_size = max_clear_root_count * (sizeof(WSSliceIntervalVector) + sizeof(WSInterval) * variable_count) +
        heap_size * (sizeof(HeapItem) + sizeof(WSInterval) * variable_count +
            sizeof(WSInterval) * term_value_count + sizeof(WSInterval) * term_partial_derivative_count * variable_count +
            sizeof(WSReal) * term_linear_count * variable_count + sizeof(WSInterval) * term_linear_count +
            sizeof(int) * term_value_count);
    if (cache_size > m_cache.GetSize()) {
        m_cache.Resize(cache_size);
    }
    m_clear_root_object_offset = 0;
    m_clear_roots = (WSSliceIntervalVector*)m_cache.Data(m_clear_root_object_offset);
    m_clear_root_data_offset = sizeof(WSSliceIntervalVector) * max_clear_root_count;
    m_heap_object_offset = m_clear_root_data_offset + sizeof(WSInterval) * variable_count * max_clear_root_count;
    m_heap = (HeapItem*)m_cache.Data(m_heap_object_offset);
    m_heap_data_offset = m_heap_object_offset + sizeof(HeapItem) * heap_size;
    m_iterate_heap_item = m_heap + max_fuzzy_domain_count;
    int cache_offset1 = m_heap_object_offset;
    int cache_offset2 = m_heap_data_offset;
    int cache_offset3 = cache_offset2 + sizeof(WSInterval) * variable_count * heap_size;
    int cache_offset4 = cache_offset3 + sizeof(WSInterval) * term_value_count * heap_size;
    int cache_offset5 = cache_offset4 + sizeof(WSInterval) * term_partial_derivative_count * variable_count * heap_size;
    int cache_offset6 = cache_offset5 + sizeof(WSReal) * term_linear_count * variable_count * heap_size;
    int cache_offset7 = cache_offset6 + sizeof(WSInterval) * term_linear_count * heap_size;
    cache_offset1 += sizeof(HeapItem) * max_fuzzy_domain_count;
    cache_offset2 += sizeof(WSInterval) * variable_count * max_fuzzy_domain_count;
    cache_offset3 += sizeof(WSInterval) * term_value_count * max_fuzzy_domain_count;
    cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count * max_fuzzy_domain_count;
    cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count * max_fuzzy_domain_count;
    cache_offset6 += sizeof(WSInterval) * term_linear_count * max_fuzzy_domain_count;
    cache_offset7 += sizeof(int) * term_value_count * max_fuzzy_domain_count;
    new(m_cache.Data(cache_offset1)) HeapItem(&m_cache, cache_offset2, cache_offset3, cache_offset4,
        cache_offset5, cache_offset6, cache_offset7, variable_count, term_value_count,
        term_partial_derivative_count, term_linear_count);
    ReserveHeap(1);
    m_heap[m_heap_item_count].m_variable.CopyFrom(variable_domain);
    m_heap[m_heap_item_count].m_priority0 = 0;
    m_heap[m_heap_item_count].m_priority1 = 0;
    for (int i = 0, n = m_equations->GetTermValueCacheCount(); i < n; ++i) {
        m_heap[m_heap_item_count].m_terms_dirty_flag[i] = 15;
    }
    if (m_equations->GetDeterminedNoRoot()) {
        return;
    }
    m_equations->SetEpsilonCoef(1, 1);
    int iterator_cache_size = WSIterator::GetCacheSize(m_equations);
    if (iterator_cache_size > m_iterator_cache.GetSize()) {
        m_iterator_cache.Resize(iterator_cache_size);
    }
    Iterator iterator(m_equations, &m_iterator_cache, 0);
    int terminate_early_tag;
    WSIterateResult iterate_result = Iterate(&iterator, &m_heap[0], 0.6, 30, terminate_early_tag);
    if (iterate_result == WSIterateResult::TerminateEarly && (terminate_early_tag & 2) == 2) {
        m_is_singularity = true;
        iterate_result = WSIterateResult::NoRoot;
    }
    switch (iterate_result) {
    case WSIterateResult::NoRoot: {
            return;
        }
    case WSIterateResult::ClearRoot: {
            ReserveClearRoots(m_clear_root_count + 1);
            m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
            return;
        }
    case WSIterateResult::TerminateEarly: {
            CalculateHeapItemTag(&m_heap[0], terminate_early_tag, g_curve_linear_limit_angle);
            CalculatePriority(&m_heap[0]);
            ++m_heap_item_count;
            if (m_heap_item_count == m_max_fuzzy_domain_count) {
                return;
            }
            break;
        }
    default: {
            CalculateHeapItemTag(&m_heap[0], 0, g_curve_linear_limit_angle);
            CalculatePriority(&m_heap[0]);
            ++m_heap_item_count;
            if (m_heap_item_count == m_max_fuzzy_domain_count) {
                return;
            }
        }
    }
    while (m_heap_item_count > 0) {
        int prev_tag = m_heap[0].m_tag;
        if (m_heap[0].m_priority0 == 1) {
            break;
        }
        ReserveHeap(m_heap_item_count + 1);
        WSInterval v = m_heap[0].m_variable.Get(0);
        WSEquationsCache(m_equations, &m_heap[0].m_variable, &m_heap[0].m_terms_value,
            &m_heap[0].m_terms_partial_derivative, &m_heap[0].m_terms_linear_a, &m_heap[0].m_terms_linear_b,
            m_heap[0].m_terms_dirty_flag).SetTermsDirty(0);
        CopyHeapItem(&m_heap[m_heap_item_count], &m_heap[0]);
        WSReal d = v.Middle();
        m_heap[0].m_variable.Set(0, WSInterval(v.Min, d));
        m_heap[m_heap_item_count].m_variable.Set(0, WSInterval(d, v.Max));
        int terminate_early_tag1;
        int terminate_early_tag2;
        WSIterateResult iterate_result1 = Iterate(&iterator, &m_heap[0], 0.6, 30, terminate_early_tag1);
        if (iterate_result1 == WSIterateResult::TerminateEarly && (terminate_early_tag1 & 2) == 2) {
            m_is_singularity = true;
            iterate_result1 = WSIterateResult::NoRoot;
        }
        WSIterateResult iterate_result2 = Iterate(&iterator, &m_heap[m_heap_item_count], 0.6, 30, terminate_early_tag2);
        if (iterate_result2 == WSIterateResult::TerminateEarly && (terminate_early_tag2 & 2) == 2) {
            m_is_singularity = true;
            iterate_result2 = WSIterateResult::NoRoot;
        }
        switch (iterate_result1)
        {
        case WSIterateResult::NoRoot: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        HeapPop(m_heap, m_heap_item_count);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag2, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                default: {
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                    }
                }
                break;
            }
        case WSIterateResult::ClearRoot: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 2);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            HeapPop(m_heap, m_heap_item_count);
                            return;
                        }
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag2, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                    }
                }
                break;
            }
        case WSIterateResult::TerminateEarly: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag | terminate_early_tag2, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                    }
                }
                break;
            }
        default: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag | terminate_early_tag2, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag, g_curve_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                    }
                }
            }
        }
    }
}

int WGIntersectHelper2d::PointCurveSolver::GetClearRootCount() {
    return m_clear_root_count;
}

const WSIntervalVector* WGIntersectHelper2d::PointCurveSolver::GetClearRoot(int index) {
    return m_clear_roots + index;
}

int WGIntersectHelper2d::PointCurveSolver::GetFuzzyDomainCount() {
    return m_heap_item_count;
}

const WSIntervalVector* WGIntersectHelper2d::PointCurveSolver::GetFuzzyDomain(int index) {
    return &m_heap[index].m_variable;
}

int WGIntersectHelper2d::PointCurveSolver::GetFuzzyTag(int index) {
    return m_heap[index].m_tag;
}

bool WGIntersectHelper2d::PointCurveSolver::GetIsSingularity() {
    return m_is_singularity;
}

bool WGIntersectHelper2d::PointCurveSolver::Iterator::CheckTerminateEarly(WSEquationSystem* equations, WSEquationsCache* cache, int& tag) {
    tag = 0;
    WSInterval dx, dy;
    ((PointCurveEquationSystem*)equations)->CalculateD1(cache, dx, dy);
    double a = calculate_tangent_cone_angle(dx, dy);
    if (a <= g_curve_linear_limit_angle) {
        tag |= 1;
    }
    if ((tag & 1) == 1) {
        return true;
    }
    if ((tag & 2) != 2) {
        WSInterval x, y;
        ((PointCurveEquationSystem*)equations)->CalculateD0(cache, x, y);
        double epsilon = ((PointCurveEquationSystem*)equations)->GetEquationCheckEpsilon(cache, 0) * 0.01;
        if (x.Length() <= epsilon && y.Length() <= epsilon) {
            tag |= 2;
        }
    }
    if ((tag & 2) == 2) {
        return true;
    }
    return false;
}

void WGIntersectHelper2d::PointCurveSolver::ReserveClearRoots(int clear_root_capacity) {
    if (clear_root_capacity > m_clear_root_capacity) {
        if (clear_root_capacity > m_max_clear_root_count) {
            clear_root_capacity = m_max_clear_root_count;
        }
        int variable_count = m_equations->GetVariableCount();
        int cache_offset1 = m_clear_root_object_offset + sizeof(WSSliceIntervalVector) * m_clear_root_capacity;
        int cache_offset2 = m_clear_root_data_offset + sizeof(WSInterval) * variable_count * m_clear_root_capacity;
        for (int i = m_clear_root_capacity; i < clear_root_capacity; ++i) {
            new(m_cache.Data(cache_offset1)) WSSliceIntervalVector(&m_cache, cache_offset2, variable_count);
            cache_offset1 += sizeof(WSSliceIntervalVector);
            cache_offset2 += sizeof(WSInterval) * variable_count;
        }
        m_clear_root_capacity = clear_root_capacity;
    }
}

void WGIntersectHelper2d::PointCurveSolver::ReserveHeap(int heap_capacity) {
    if (heap_capacity > m_heap_capacity) {
        if (heap_capacity > m_max_fuzzy_domain_count) {
            heap_capacity = m_max_fuzzy_domain_count;
        }
        int variable_count = m_equations->GetVariableCount();
        int term_value_count = m_equations->GetTermValueCacheCount();
        int term_partial_derivative_count = m_equations->GetTermPartialDerivativeCacheCount();
        int term_linear_count = m_equations->GetTermLinearCacheCount();
        int heap_size = m_max_fuzzy_domain_count + 1;
        int cache_offset1 = m_heap_object_offset;
        int cache_offset2 = m_heap_data_offset;
        int cache_offset3 = cache_offset2 + sizeof(WSInterval) * variable_count * heap_size;
        int cache_offset4 = cache_offset3 + sizeof(WSInterval) * term_value_count * heap_size;
        int cache_offset5 = cache_offset4 + sizeof(WSInterval) * term_partial_derivative_count * variable_count * heap_size;
        int cache_offset6 = cache_offset5 + sizeof(WSReal) * term_linear_count * variable_count * heap_size;
        int cache_offset7 = cache_offset6 + sizeof(WSInterval) * term_linear_count * heap_size;
        cache_offset1 += sizeof(HeapItem) * m_heap_capacity;
        cache_offset2 += sizeof(WSInterval) * variable_count * m_heap_capacity;
        cache_offset3 += sizeof(WSInterval) * term_value_count * m_heap_capacity;
        cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count * m_heap_capacity;
        cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count * m_heap_capacity;
        cache_offset6 += sizeof(WSInterval) * term_linear_count * m_heap_capacity;
        cache_offset7 += sizeof(int) * term_value_count * m_heap_capacity;
        for (int i = m_heap_capacity; i < heap_capacity; ++i) {
            new(m_cache.Data(cache_offset1)) HeapItem(&m_cache, cache_offset2, cache_offset3, cache_offset4,
                cache_offset5, cache_offset6, cache_offset7, variable_count, term_value_count,
                term_partial_derivative_count, term_linear_count);
            cache_offset1 += sizeof(HeapItem);
            cache_offset2 += sizeof(WSInterval) * variable_count;
            cache_offset3 += sizeof(WSInterval) * term_value_count;
            cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count;
            cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count;
            cache_offset6 += sizeof(WSInterval) * term_linear_count;
            cache_offset7 += sizeof(int) * term_value_count;
        }
        m_heap_capacity = heap_capacity;
    }
}

WSIterateResult WGIntersectHelper2d::PointCurveSolver::Iterate(WSIterator* iterator, HeapItem* heap_item, double min_rate,
    int max_iterate_count, int& terminate_early_tag) {
    WSEquationsCache cache(m_equations, &heap_item->m_variable, &heap_item->m_terms_value,
        &heap_item->m_terms_partial_derivative, &heap_item->m_terms_linear_a, &heap_item->m_terms_linear_b,
        heap_item->m_terms_dirty_flag);
    int iterate_count;
    WSIterateResult r = iterator->Execute(&cache, min_rate, max_iterate_count, iterate_count, terminate_early_tag);
    return r;
}

void WGIntersectHelper2d::PointCurveSolver::CopyHeapItem(HeapItem* dst, HeapItem* src) {
    dst->m_variable.CopyFrom(&src->m_variable);
    dst->m_terms_value.CopyFrom(&src->m_terms_value);
    dst->m_terms_partial_derivative.CopyFrom(&src->m_terms_partial_derivative);
    dst->m_terms_linear_a.CopyFrom(&src->m_terms_linear_a);
    dst->m_terms_linear_b.CopyFrom(&src->m_terms_linear_b);
    memcpy(dst->m_terms_dirty_flag, src->m_terms_dirty_flag, m_equations->GetTermValueCacheCount() * sizeof(int));
    dst->m_tag = src->m_tag;
    dst->m_priority0 = src->m_priority0;
    dst->m_priority1 = src->m_priority1;
}

void WGIntersectHelper2d::PointCurveSolver::CalculateHeapItemTag(HeapItem* heap_item, int prev_tag, double limit_angle) {
    heap_item->m_tag = prev_tag;
    if ((prev_tag & 1) == 0) {
        WSInterval dx, dy;
        m_equations->CalculateD1(&heap_item->m_variable, dx, dy);
        double a = calculate_tangent_cone_angle(dx, dy);
        if (a <= limit_angle) {
            heap_item->m_tag |= 1;
        }
    }
}

void WGIntersectHelper2d::PointCurveSolver::CalculatePriority(HeapItem* heap_item) {
    if ((heap_item->m_tag & 1) == 1) {
        heap_item->m_priority0 = 1;
    }
    else {
        heap_item->m_priority0 = 0;
    }
    heap_item->m_priority1 = heap_item->m_variable.Get(0).Length();
}

WGIntersectHelper2d::CurveCurveSolver::CurveCurveSolver() :
    m_equations(nullptr),
    m_is_singularity0(false),
    m_is_singularity1(false),
    m_heap_item_count(0),
    m_clear_root_count(0) {
}

void WGIntersectHelper2d::CurveCurveSolver::PrepareExecute(CurveCurveEquationSystem* equations, const WSIntervalVector* variable_domain,
    int max_clear_root_count, int max_fuzzy_domain_count) {
    m_equations = equations;
    m_max_clear_root_count = max_clear_root_count;
    m_max_fuzzy_domain_count = max_fuzzy_domain_count;
    m_clear_root_capacity = 0;
    m_heap_capacity = 0;
    m_clear_root_count = 0;
    m_heap_item_count = 0;
    m_is_singularity0 = false;
    m_is_singularity1 = false;
    assert(max_clear_root_count > 0 && max_fuzzy_domain_count > 0);
    int variable_count = equations->GetVariableCount();
    int term_value_count = equations->GetTermValueCacheCount();
    int term_partial_derivative_count = equations->GetTermPartialDerivativeCacheCount();
    int term_linear_count = equations->GetTermLinearCacheCount();
    int heap_size = max_fuzzy_domain_count + 1;
    int cache_size = max_clear_root_count * (sizeof(WSSliceIntervalVector) + sizeof(WSInterval) * variable_count) +
        heap_size * (sizeof(HeapItem) + sizeof(WSInterval) * variable_count +
            sizeof(WSInterval) * term_value_count + sizeof(WSInterval) * term_partial_derivative_count * variable_count +
            sizeof(WSReal) * term_linear_count * variable_count + sizeof(WSInterval) * term_linear_count +
            sizeof(int) * term_value_count);
    if (cache_size > m_cache.GetSize()) {
        m_cache.Resize(cache_size);
    }
    m_clear_root_object_offset = 0;
    m_clear_roots = (WSSliceIntervalVector*)m_cache.Data(m_clear_root_object_offset);
    m_clear_root_data_offset = sizeof(WSSliceIntervalVector) * max_clear_root_count;
    m_heap_object_offset = m_clear_root_data_offset + sizeof(WSInterval) * variable_count * max_clear_root_count;
    m_heap = (HeapItem*)m_cache.Data(m_heap_object_offset);
    m_heap_data_offset = m_heap_object_offset + sizeof(HeapItem) * heap_size;
    m_iterate_heap_item = m_heap + max_fuzzy_domain_count;
    int cache_offset1 = m_heap_object_offset;
    int cache_offset2 = m_heap_data_offset;
    int cache_offset3 = cache_offset2 + sizeof(WSInterval) * variable_count * heap_size;
    int cache_offset4 = cache_offset3 + sizeof(WSInterval) * term_value_count * heap_size;
    int cache_offset5 = cache_offset4 + sizeof(WSInterval) * term_partial_derivative_count * variable_count * heap_size;
    int cache_offset6 = cache_offset5 + sizeof(WSReal) * term_linear_count * variable_count * heap_size;
    int cache_offset7 = cache_offset6 + sizeof(WSInterval) * term_linear_count * heap_size;
    cache_offset1 += sizeof(HeapItem) * max_fuzzy_domain_count;
    cache_offset2 += sizeof(WSInterval) * variable_count * max_fuzzy_domain_count;
    cache_offset3 += sizeof(WSInterval) * term_value_count * max_fuzzy_domain_count;
    cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count * max_fuzzy_domain_count;
    cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count * max_fuzzy_domain_count;
    cache_offset6 += sizeof(WSInterval) * term_linear_count * max_fuzzy_domain_count;
    cache_offset7 += sizeof(int) * term_value_count * max_fuzzy_domain_count;
    new(m_cache.Data(cache_offset1)) HeapItem(&m_cache, cache_offset2, cache_offset3, cache_offset4,
        cache_offset5, cache_offset6, cache_offset7, variable_count, term_value_count,
        term_partial_derivative_count, term_linear_count);
    ReserveHeap(1);
    m_heap[m_heap_item_count].m_variable.CopyFrom(variable_domain);
    m_heap[m_heap_item_count].m_priority0 = 0;
    m_heap[m_heap_item_count].m_priority1 = 0;
    for (int i = 0, n = m_equations->GetTermValueCacheCount(); i < n; ++i) {
        m_heap[m_heap_item_count].m_terms_dirty_flag[i] = 15;
    }
}

void WGIntersectHelper2d::CurveCurveSolver::PrimaryExecute() {
    if (m_equations->GetDeterminedNoRoot()) {
        return;
    }
    m_equations->SetEpsilonCoef(1, 1);
    int iterator_cache_size = WSIterator::GetCacheSize(m_equations);
    if (iterator_cache_size > m_iterator_cache.GetSize()) {
        m_iterator_cache.Resize(iterator_cache_size);
    }
    PrimaryIterator iterator(m_equations, &m_iterator_cache, 0);
    int terminate_early_tag;
    WSIterateResult iterate_result = PrimaryIterate(&iterator, &m_heap[0], 0.6, 30, terminate_early_tag);
    if (iterate_result == WSIterateResult::TerminateEarly) {
        if ((terminate_early_tag & 8) == 8) {
            m_is_singularity0 = true;
            iterate_result = WSIterateResult::NoRoot;
        }
        if ((terminate_early_tag & 16) == 16) {
            m_is_singularity1 = true;
            iterate_result = WSIterateResult::NoRoot;
        }
    }
    switch (iterate_result) {
    case WSIterateResult::NoRoot: {
            return;
        }
    case WSIterateResult::ClearRoot: {
            ReserveClearRoots(m_clear_root_count + 1);
            m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
            return;
        }
    case WSIterateResult::TerminateEarly: {
            CalculateHeapItemTag(&m_heap[0], terminate_early_tag, g_curve_primary_linear_limit_angle);
            CalculatePriority(&m_heap[0]);
            ++m_heap_item_count;
            if (m_heap_item_count == m_max_fuzzy_domain_count) {
                return;
            }
            break;
        }
    default: {
            CalculateHeapItemTag(&m_heap[0], 0, g_curve_primary_linear_limit_angle);
            CalculatePriority(&m_heap[0]);
            ++m_heap_item_count;
            if (m_heap_item_count == m_max_fuzzy_domain_count) {
                return;
            }
        }
    }
    while (m_heap_item_count > 0) {
        int prev_tag = m_heap[0].m_tag;
        if (m_heap[0].m_priority0 == 1) {
            break;
        }
        ReserveHeap(m_heap_item_count + 1);
        int split_index = m_equations->GetSplitIndex(&m_heap[0].m_variable);
        WSInterval v = m_heap[0].m_variable.Get(split_index);
        WSEquationsCache(m_equations, &m_heap[0].m_variable, &m_heap[0].m_terms_value,
            &m_heap[0].m_terms_partial_derivative, &m_heap[0].m_terms_linear_a, &m_heap[0].m_terms_linear_b,
            m_heap[0].m_terms_dirty_flag).SetTermsDirty(split_index);
        CopyHeapItem(&m_heap[m_heap_item_count], &m_heap[0]);
        WSReal d = v.Middle();
        m_heap[0].m_variable.Set(split_index, WSInterval(v.Min, d));
        m_heap[m_heap_item_count].m_variable.Set(split_index, WSInterval(d, v.Max));
        int terminate_early_tag1;
        int terminate_early_tag2;
        WSIterateResult iterate_result1 = PrimaryIterate(&iterator, &m_heap[0], 0.6, 30, terminate_early_tag1);
        if (iterate_result1 == WSIterateResult::TerminateEarly) {
            if ((terminate_early_tag1 & 8) == 8) {
                m_is_singularity0 = true;
                iterate_result1 = WSIterateResult::NoRoot;
            }
            if ((terminate_early_tag1 & 16) == 16) {
                m_is_singularity1 = true;
                iterate_result1 = WSIterateResult::NoRoot;
            }
        }
        WSIterateResult iterate_result2 = PrimaryIterate(&iterator, &m_heap[m_heap_item_count], 0.6, 30, terminate_early_tag2);
        if (iterate_result2 == WSIterateResult::TerminateEarly) {
            if ((terminate_early_tag2 & 8) == 8) {
                m_is_singularity0 = true;
                iterate_result2 = WSIterateResult::NoRoot;
            }
            if ((terminate_early_tag2 & 16) == 16) {
                m_is_singularity1 = true;
                iterate_result2 = WSIterateResult::NoRoot;
            }
        }
        switch (iterate_result1)
        {
        case WSIterateResult::NoRoot: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        HeapPop(m_heap, m_heap_item_count);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag2, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                default: {
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                    }
                }
                break;
            }
        case WSIterateResult::ClearRoot: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 2);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            HeapPop(m_heap, m_heap_item_count);
                            return;
                        }
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        HeapPop(m_heap, m_heap_item_count);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag2, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[0].m_variable);
                        HeapItem temp_heap_item = m_heap[0];
                        m_heap[0] = m_heap[m_heap_item_count];
                        m_heap[m_heap_item_count] = temp_heap_item;
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                    }
                }
                break;
            }
        case WSIterateResult::TerminateEarly: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag | terminate_early_tag2, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag | terminate_early_tag1, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                    }
                }
                break;
            }
        default: {
                switch (iterate_result2)
                {
                case WSIterateResult::NoRoot: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        break;
                    }
                case WSIterateResult::ClearRoot: {
                        ReserveClearRoots(m_clear_root_count + 1);
                        m_clear_roots[m_clear_root_count++].CopyFrom(&m_heap[m_heap_item_count].m_variable);
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        if (m_clear_root_count == m_max_clear_root_count) {
                            return;
                        }
                        break;
                    }
                case WSIterateResult::TerminateEarly: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag | terminate_early_tag2, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                        break;
                    }
                default: {
                        CalculateHeapItemTag(&m_heap[0], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[0]);
                        HeapMoveDown(m_heap, m_heap_item_count, 0);
                        CalculateHeapItemTag(&m_heap[m_heap_item_count], prev_tag, g_curve_primary_linear_limit_angle);
                        CalculatePriority(&m_heap[m_heap_item_count]);
                        ++m_heap_item_count;
                        HeapMoveUp(m_heap, m_heap_item_count, m_heap_item_count - 1);
                        if (m_heap_item_count == m_max_fuzzy_domain_count) {
                            return;
                        }
                    }
                }
            }
        }
    }
}

WSIterateResult WGIntersectHelper2d::CurveCurveSolver::FinallyIterate(int heap_index, double distance_epsilon) {
    m_equations->InitializeIntermediateVariables(&m_heap[heap_index].m_variable);
    m_equations->SetAlgebraEnable(&m_heap[heap_index].m_variable, distance_epsilon);
    if (m_equations->GetDeterminedNoRoot()) {
        return WSIterateResult::NoRoot;
    }
    m_equations->SetEpsilonCoef(1, 1);
    int iterator_cache_size = WSIterator::GetCacheSize(m_equations);
    if (iterator_cache_size > m_iterator_cache.GetSize()) {
        m_iterator_cache.Resize(iterator_cache_size);
    }
    WSIterator iterator(m_equations, &m_iterator_cache, 0);
    return FinallyIterate(&iterator, &m_heap[heap_index], 0.8, 30);
}

WSIterateResult WGIntersectHelper2d::CurveCurveSolver::Iterate(WSIntervalVector* variable, double check_coef, double iterate_coef) {
    m_iterate_heap_item->m_variable.CopyFrom(variable);
    for (int i = 0, n = m_equations->GetTermValueCacheCount(); i < n; ++i) {
        m_iterate_heap_item->m_terms_dirty_flag[i] = 15;
    }
    m_equations->SetEpsilonCoef(check_coef, iterate_coef);
    int iterator_cache_size = WSIterator::GetCacheSize(m_equations);
    if (iterator_cache_size > m_iterator_cache.GetSize()) {
        m_iterator_cache.Resize(iterator_cache_size);
    }
    WSIterator iterator(m_equations, &m_iterator_cache, 0);
    WSEquationsCache cache(m_equations, &m_iterate_heap_item->m_variable, &m_iterate_heap_item->m_terms_value,
        &m_iterate_heap_item->m_terms_partial_derivative, &m_iterate_heap_item->m_terms_linear_a, 
        &m_iterate_heap_item->m_terms_linear_b, m_iterate_heap_item->m_terms_dirty_flag);
    int iterate_count;
    int terminate_early_tag;
    WSIterateResult iterate_result = iterator.Execute(&cache, 0.8, 30, iterate_count, terminate_early_tag);
    if (iterate_result != WSIterateResult::NoRoot) {
        variable->Set(0, m_iterate_heap_item->m_variable.Get(0));
        variable->Set(1, m_iterate_heap_item->m_variable.Get(1));
    }
    return iterate_result;
}

int WGIntersectHelper2d::CurveCurveSolver::GetClearRootCount() {
    return m_clear_root_count;
}

const WSIntervalVector* WGIntersectHelper2d::CurveCurveSolver::GetClearRoot(int index) {
    return m_clear_roots + index;
}

int WGIntersectHelper2d::CurveCurveSolver::GetFuzzyDomainCount() {
    return m_heap_item_count;
}

const WSIntervalVector* WGIntersectHelper2d::CurveCurveSolver::GetFuzzyDomain(int index) {
    return &m_heap[index].m_variable;
}

int WGIntersectHelper2d::CurveCurveSolver::GetFuzzyTag(int index) {
    return m_heap[index].m_tag;
}

bool WGIntersectHelper2d::CurveCurveSolver::GetIsSingularity0() {
    return m_is_singularity0;
}

bool WGIntersectHelper2d::CurveCurveSolver::GetIsSingularity1() {
    return m_is_singularity1;
}

bool WGIntersectHelper2d::CurveCurveSolver::PrimaryIterator::CheckTerminateEarly(WSEquationSystem* equations, WSEquationsCache* cache, int& tag) {
    tag = 0;
    WSInterval dx0, dy0, dx1, dy1;
    ((CurveCurveEquationSystem*)equations)->CalculateD1(0, cache, dx0, dy0);
    ((CurveCurveEquationSystem*)equations)->CalculateD1(1, cache, dx1, dy1);
    WSInterval d = dx0 * dy1 - dx1 * dy0;
    if (d.Min > 0 || d.Max < 0) {
        tag |= 1;
    }
    double a0 = calculate_tangent_cone_angle(dx0, dy0);
    if (a0 <= g_curve_primary_linear_limit_angle) {
        tag |= 2;
    }
    double a1 = calculate_tangent_cone_angle(dx1, dy1);
    if (a1 <= g_curve_primary_linear_limit_angle) {
        tag |= 4;
    }
    if ((tag & 1) == 1 || (tag & 6) == 6) {
        return true;
    }
    if ((tag & 8) != 8) {
        WSInterval x0, y0;
        ((CurveCurveEquationSystem*)equations)->CalculateD0(0, cache, x0, y0);
        double epsilon = ((CurveCurveEquationSystem*)equations)->GetEquationCheckEpsilon(cache, 0) * 0.01;
        if (x0.Length() <= epsilon && y0.Length() <= epsilon) {
            tag |= 8;
        }
    }
    if ((tag & 16) != 16) {
        WSInterval x1, y1;
        ((CurveCurveEquationSystem*)equations)->CalculateD0(1, cache, x1, y1);
        double epsilon = ((CurveCurveEquationSystem*)equations)->GetEquationCheckEpsilon(cache, 1) * 0.01;
        if (x1.Length() <= epsilon && y1.Length() <= epsilon) {
            tag |= 16;
        }
    }
    if ((tag & 8) == 8 || (tag & 16) == 16) {
        return true;
    }
    return false;
}

void WGIntersectHelper2d::CurveCurveSolver::ReserveClearRoots(int clear_root_capacity) {
    if (clear_root_capacity > m_clear_root_capacity) {
        if (clear_root_capacity > m_max_clear_root_count) {
            clear_root_capacity = m_max_clear_root_count;
        }
        int variable_count = m_equations->GetVariableCount();
        int cache_offset1 = m_clear_root_object_offset + sizeof(WSSliceIntervalVector) * m_clear_root_capacity;
        int cache_offset2 = m_clear_root_data_offset + sizeof(WSInterval) * variable_count * m_clear_root_capacity;
        for (int i = m_clear_root_capacity; i < clear_root_capacity; ++i) {
            new(m_cache.Data(cache_offset1)) WSSliceIntervalVector(&m_cache, cache_offset2, variable_count);
            cache_offset1 += sizeof(WSSliceIntervalVector);
            cache_offset2 += sizeof(WSInterval) * variable_count;
        }
        m_clear_root_capacity = clear_root_capacity;
    }
}

void WGIntersectHelper2d::CurveCurveSolver::ReserveHeap(int heap_capacity) {
    if (heap_capacity > m_heap_capacity) {
        if (heap_capacity > m_max_fuzzy_domain_count) {
            heap_capacity = m_max_fuzzy_domain_count;
        }
        int variable_count = m_equations->GetVariableCount();
        int term_value_count = m_equations->GetTermValueCacheCount();
        int term_partial_derivative_count = m_equations->GetTermPartialDerivativeCacheCount();
        int term_linear_count = m_equations->GetTermLinearCacheCount();
        int heap_size = m_max_fuzzy_domain_count + 1;
        int cache_offset1 = m_heap_object_offset;
        int cache_offset2 = m_heap_data_offset;
        int cache_offset3 = cache_offset2 + sizeof(WSInterval) * variable_count * heap_size;
        int cache_offset4 = cache_offset3 + sizeof(WSInterval) * term_value_count * heap_size;
        int cache_offset5 = cache_offset4 + sizeof(WSInterval) * term_partial_derivative_count * variable_count * heap_size;
        int cache_offset6 = cache_offset5 + sizeof(WSReal) * term_linear_count * variable_count * heap_size;
        int cache_offset7 = cache_offset6 + sizeof(WSInterval) * term_linear_count * heap_size;
        cache_offset1 += sizeof(HeapItem) * m_heap_capacity;
        cache_offset2 += sizeof(WSInterval) * variable_count * m_heap_capacity;
        cache_offset3 += sizeof(WSInterval) * term_value_count * m_heap_capacity;
        cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count * m_heap_capacity;
        cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count * m_heap_capacity;
        cache_offset6 += sizeof(WSInterval) * term_linear_count * m_heap_capacity;
        cache_offset7 += sizeof(int) * term_value_count * m_heap_capacity;
        for (int i = m_heap_capacity; i < heap_capacity; ++i) {
            new(m_cache.Data(cache_offset1)) HeapItem(&m_cache, cache_offset2, cache_offset3, cache_offset4,
                cache_offset5, cache_offset6, cache_offset7, variable_count, term_value_count,
                term_partial_derivative_count, term_linear_count);
            cache_offset1 += sizeof(HeapItem);
            cache_offset2 += sizeof(WSInterval) * variable_count;
            cache_offset3 += sizeof(WSInterval) * term_value_count;
            cache_offset4 += sizeof(WSInterval) * term_partial_derivative_count * variable_count;
            cache_offset5 += sizeof(WSReal) * term_linear_count * variable_count;
            cache_offset6 += sizeof(WSInterval) * term_linear_count;
            cache_offset7 += sizeof(int) * term_value_count;
        }
        m_heap_capacity = heap_capacity;
    }
}

WSIterateResult WGIntersectHelper2d::CurveCurveSolver::PrimaryIterate(WSIterator* iterator, HeapItem* heap_item, 
    double min_rate, int max_iterate_count, int& terminate_early_tag) {
    WSEquationsCache cache(m_equations, &heap_item->m_variable, &heap_item->m_terms_value,
        &heap_item->m_terms_partial_derivative, &heap_item->m_terms_linear_a, &heap_item->m_terms_linear_b,
        heap_item->m_terms_dirty_flag);
    int iterate_count;
    WSIterateResult r = iterator->Execute(&cache, min_rate, max_iterate_count, iterate_count, terminate_early_tag);
    return r;
}

WSIterateResult WGIntersectHelper2d::CurveCurveSolver::FinallyIterate(WSIterator* iterator, HeapItem* heap_item, double min_rate, int max_iterate_count) {
    WSEquationsCache cache(m_equations, &heap_item->m_variable, &heap_item->m_terms_value,
        &heap_item->m_terms_partial_derivative, &heap_item->m_terms_linear_a, &heap_item->m_terms_linear_b,
        heap_item->m_terms_dirty_flag);
    int iterate_count;
    int terminate_early_tag;
    WSIterateResult r = iterator->Execute(&cache, min_rate, max_iterate_count, iterate_count, terminate_early_tag);
    return r;
}

void WGIntersectHelper2d::CurveCurveSolver::CopyHeapItem(HeapItem* dst, HeapItem* src) {
    dst->m_variable.CopyFrom(&src->m_variable);
    dst->m_terms_value.CopyFrom(&src->m_terms_value);
    dst->m_terms_partial_derivative.CopyFrom(&src->m_terms_partial_derivative);
    dst->m_terms_linear_a.CopyFrom(&src->m_terms_linear_a);
    dst->m_terms_linear_b.CopyFrom(&src->m_terms_linear_b);
    memcpy(dst->m_terms_dirty_flag, src->m_terms_dirty_flag, m_equations->GetTermValueCacheCount() * sizeof(int));
    dst->m_tag = src->m_tag;
    dst->m_priority0 = src->m_priority0;
    dst->m_priority1 = src->m_priority1;
}

void WGIntersectHelper2d::CurveCurveSolver::CalculateHeapItemTag(HeapItem* heap_item, int prev_tag, double limit_angle) {
    heap_item->m_tag = prev_tag;
    if ((prev_tag & 1) == 0) {
        WSInterval dx0, dy0, dx1, dy1;
        m_equations->CalculateD1(0, &heap_item->m_variable, dx0, dy0);
        m_equations->CalculateD1(1, &heap_item->m_variable, dx1, dy1);
        WSInterval d = dx0 * dy1 - dx1 * dy0;
        if (d.Min > 0 || d.Max < 0) {
            heap_item->m_tag |= 1;
        }
        if ((prev_tag & 2) == 0) {
            double a0 = calculate_tangent_cone_angle(dx0, dy0);
            if (a0 <= limit_angle) {
                heap_item->m_tag |= 2;
            }
        }
        if ((prev_tag & 4) == 0) {
            double a1 = calculate_tangent_cone_angle(dx1, dy1);
            if (a1 <= limit_angle) {
                heap_item->m_tag |= 4;
            }
        }
    }
}

void WGIntersectHelper2d::CurveCurveSolver::CalculatePriority(HeapItem* heap_item) {
    if ((heap_item->m_tag & 6) == 6) {
        heap_item->m_priority0 = 1;
    }
    else if ((heap_item->m_tag & 1) == 1) {
        heap_item->m_priority0 = 1;
    }
    else {
        heap_item->m_priority0 = 0;
    }
    heap_item->m_priority1 = heap_item->m_variable.Get(0).Length() * heap_item->m_variable.Get(1).Length();
}

WGIntersectHelper2d::PointBezierCurveEquationSystem::PointBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::PointBezierCurveEquationSystem::~PointBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::Update(const WGVector2d& point, 
    int degree, const WGVector2d* control_points, double distance_epsilon) {
    /*
    v0:     t
    ---------------------
    b0:     x
    b1:     y
    b2:     const
    ---------------------
    e0:     x - point.x = 0
    e1:     y - point.y = 0
    */
    m_distance_epsilon = distance_epsilon;
    if (!m_initialized) {
        m_basises[0] = new WSBernsteinBasis(0, degree);
        m_basises[1] = new WSBernsteinBasis(0, degree);
        m_basises[2] = new WSConstBasis();
        m_equations[0] = (new WSEquation(this, 2, true, true, true))->AddTerm(0, 1)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 2, true, true, true))->AddTerm(1, 1)->AddTerm(2, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_term_x[0] = term;
            }
            else if (term->GetIndex() == 2) {
                m_term_x[1] = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_term_y[0] = term;
            }
            else if (term->GetIndex() == 2) {
                m_term_y[1] = term;
            }
        }
        m_initialized = true;
    }
    else {
        if (((WSBernsteinBasis*)m_basises[0])->GetDegree() != degree) {
            delete m_basises[0];
            delete m_basises[1];
            m_basises[0] = new WSBernsteinBasis(0, degree);
            m_basises[1] = new WSBernsteinBasis(0, degree);
        }
    }
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[0])->SetCoef(i, control_points[i].X);
        ((WSBernsteinBasis*)m_basises[1])->SetCoef(i, control_points[i].Y);
    }
    m_term_x[1]->SetCoef(-point.X);
    m_term_y[1]->SetCoef(-point.Y);
}

int WGIntersectHelper2d::PointBezierCurveEquationSystem::GetVariableCount() const {
    return 1;
}

WSReal WGIntersectHelper2d::PointBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::PointBezierCurveEquationSystem::GetBasisCount() const {
    return 3;
}

WSEquationBasis* WGIntersectHelper2d::PointBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::PointBezierCurveEquationSystem::GetEquationCount() const {
    return 2;
}

WSEquation* WGIntersectHelper2d::PointBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::PointBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    return m_distance_epsilon * m_check_coef;
}

WSReal WGIntersectHelper2d::PointBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    return m_distance_epsilon * m_iterate_coef;
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    polynomials = nullptr;
    polynomial_count = 0;
    leader_variables = nullptr;
    leader_variable_count = 0;
    max_term_count = 0;
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::CalculateD0(double t, WGVector2d& d0) const {
    d0.X = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
    d0.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::CalculateD1(double t, WGVector2d& d1) const {
    d1.X = ((WSBernsteinBasis*)m_basises[0])->CalculateDerivative(t);
    d1.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateDerivative(t);
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::CalculateD1(const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    dx = ((WSBernsteinBasis*)m_basises[0])->CalculatePartialDerivative(variable, 0);
    dy = ((WSBernsteinBasis*)m_basises[1])->CalculatePartialDerivative(variable, 0);
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::CalculateD0(WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    x = cache->GetTermsValue(m_term_x[0]);
    y = cache->GetTermsValue(m_term_y[0]);
}

void WGIntersectHelper2d::PointBezierCurveEquationSystem::CalculateD1(WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    dx = cache->GetTermsPartialDerivative(m_term_x[0], 0);
    dy = cache->GetTermsPartialDerivative(m_term_y[0], 0);
}

WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::PointRationalBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::~PointRationalBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::Update(const WGVector2d& point,
    int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon) {
    /*
    v0:     t
    ---------------------
    b0:     x
    b1:     y
    b2:     w
    ---------------------
    e0:     x - point.x * w = 0
    e1:     y - point.y * w = 0
    */
    m_distance_epsilon = distance_epsilon;
    if (!m_initialized) {
        m_basises[0] = new WSBernsteinBasis(0, degree);
        m_basises[1] = new WSBernsteinBasis(0, degree);
        m_basises[2] = new WSBernsteinBasis(0, degree);
        m_equations[0] = (new WSEquation(this, 2, true, true, true))->AddTerm(0, 1)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 2, true, true, true))->AddTerm(1, 1)->AddTerm(2, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_term_x[0] = term;
            }
            else if (term->GetIndex() == 2) {
                m_term_x[1] = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_term_y[0] = term;
            }
            else if (term->GetIndex() == 2) {
                m_term_y[1] = term;
            }
        }
        m_initialized = true;
    }
    else {
        if (((WSBernsteinBasis*)m_basises[0])->GetDegree() != degree) {
            delete m_basises[0];
            delete m_basises[1];
            delete m_basises[2];
            m_basises[0] = new WSBernsteinBasis(0, degree);
            m_basises[1] = new WSBernsteinBasis(0, degree);
            m_basises[2] = new WSBernsteinBasis(0, degree);
        }
    }
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[0])->SetCoef(i, control_points[i].X * weights[i]);
        ((WSBernsteinBasis*)m_basises[1])->SetCoef(i, control_points[i].Y * weights[i]);
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, weights[i]);
    }
    m_term_x[1]->SetCoef(-point.X);
    m_term_y[1]->SetCoef(-point.Y);
}

int WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetVariableCount() const {
    return 1;
}

WSReal WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetBasisCount() const {
    return 3;
}

WSEquationBasis* WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetEquationCount() const {
    return 2;
}

WSEquation* WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(cache->GetVariable());
    double epsilon = m_distance_epsilon * w.Min * m_check_coef;
    if (epsilon < g_double_epsilon) {
        epsilon = g_double_epsilon;
    }
    return epsilon;
}

WSReal WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(cache->GetVariable());
    double epsilon = m_distance_epsilon * w.Max * m_iterate_coef;
    if (epsilon < g_double_epsilon) {
        epsilon = g_double_epsilon;
    }
    return epsilon;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    polynomials = nullptr;
    polynomial_count = 0;
    leader_variables = nullptr;
    leader_variable_count = 0;
    max_term_count = 0;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::CalculateD0(double t, WGVector2d& d0) const {
    double x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
    double y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
    double w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
    d0.X = x / w;
    d0.Y = y / w;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::CalculateD1(double t, WGVector2d& d1) const {
    double x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
    double y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
    double w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
    double dx = ((WSBernsteinBasis*)m_basises[0])->CalculateDerivative(t);
    double dy = ((WSBernsteinBasis*)m_basises[1])->CalculateDerivative(t);
    double dw = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
    double w2 = w * w;
    d1.X = (dx * w - x * dw) / w2;
    d1.Y = (dy * w - y * dw) / w2;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::CalculateD1(const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    WSInterval x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(variable);
    WSInterval y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(variable);
    WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
    dx = ((WSBernsteinBasis*)m_basises[0])->CalculatePartialDerivative(variable, 0);
    dy = ((WSBernsteinBasis*)m_basises[1])->CalculatePartialDerivative(variable, 0);
    WSInterval dw = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 0);
    WSInterval w2 = pow(w, 2);
    dx = (dx * w - x * dw) / w2;
    dy = (dy * w - y * dw) / w2;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::CalculateD0(WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    x = cache->GetTermsValue(m_term_x[0]);
    y = cache->GetTermsValue(m_term_y[0]);
    WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(cache->GetVariable());
    x = x / w;
    y = y / w;
}

void WGIntersectHelper2d::PointRationalBezierCurveEquationSystem::CalculateD1(WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    WSInterval x = cache->GetTermsValue(m_term_x[0]);
    WSInterval y = cache->GetTermsValue(m_term_y[0]);
    WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(cache->GetVariable());
    dx = cache->GetTermsPartialDerivative(m_term_x[0], 0);
    dy = cache->GetTermsPartialDerivative(m_term_y[0], 0);
    WSInterval dw = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(cache->GetVariable(), 0);
    WSInterval w2 = pow(w, 2);
    dx = (dx * w - x * dw) / w2;
    dy = (dy * w - y * dw) / w2;
}

WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::BezierCurveBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::~BezierCurveBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::Update(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    ---------------------
    b0:     x0
    b1:     y0
    b2:     x1
    b3:     y1
    ---------------------
    e0:     x0 - x1 = 0
    e1:     y0 - y1 = 0
    */
    m_distance_epsilon = distance_epsilon;
    UpdateVariableEpsilons(degree0, control_points0, degree1, control_points1, distance_epsilon);
    if (!m_initialized) {
        m_basises[0] = new WSBernsteinBasis(0, degree0);
        m_basises[1] = new WSBernsteinBasis(0, degree0);
        m_basises[2] = new WSBernsteinBasis(1, degree1);
        m_basises[3] = new WSBernsteinBasis(1, degree1);
        m_equations[0] = (new WSEquation(this, 2, true, true, true))->AddTerm(0, 1)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 2, true, true, true))->AddTerm(1, 1)->AddTerm(3, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_term_x[0] = term;
            }
            else if (term->GetIndex() == 2) {
                m_term_x[1] = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_term_y[0] = term;
            }
            else if (term->GetIndex() == 3) {
                m_term_y[1] = term;
            }
        }
        m_initialized = true;
    }
    else {
        if (((WSBernsteinBasis*)m_basises[0])->GetDegree() != degree0) {
            delete m_basises[0];
            delete m_basises[1];
            m_basises[0] = new WSBernsteinBasis(0, degree0);
            m_basises[1] = new WSBernsteinBasis(0, degree0);
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree1) {
            delete m_basises[2];
            delete m_basises[3];
            m_basises[2] = new WSBernsteinBasis(1, degree1);
            m_basises[3] = new WSBernsteinBasis(1, degree1);
        }
        ClearAlgebraRuntime();
    }
    for (int i = 0; i <= degree0; ++i) {
        ((WSBernsteinBasis*)m_basises[0])->SetCoef(i, control_points0[i].X);
        ((WSBernsteinBasis*)m_basises[1])->SetCoef(i, control_points0[i].Y);
    }
    for (int i = 0; i <= degree1; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points1[i].X);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points1[i].Y);
    }
    m_algebra_enable = 0;
}

int WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetVariableCount() const {
    return 2;
}

WSReal WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetBasisCount() const {
    return 4;
}

WSEquationBasis* WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetEquationCount() const {
    return 2;
}

WSEquation* WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_check_coef;
    }
    return g_double_epsilon;
}

WSReal WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_iterate_coef;
    }
    return g_double_epsilon;
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree0 = ((WSBernsteinBasis*)m_basises[0])->GetDegree();
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [2];
        polynomial_count = 2;

        polynomials[0] = new WSPolynomial(GetVariableCount(), degree0 + degree1 + 2);
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_basises[0]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
            else if (term->GetIndex() == 2) {
                m_basises[2]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
        }
        WSPolynomialTerm* term = polynomials[0]->GetTerm(polynomials[0]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[0]->AddLast();

        polynomials[1] = new WSPolynomial(GetVariableCount(), degree0 + degree1 + 2);
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
            else if (term->GetIndex() == 3) {
                m_basises[3]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
        }
        term = polynomials[1]->GetTerm(polynomials[1]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[1]->AddLast();
        
        leader_variables = new int[2];
        if (degree0 <= degree1) {
            leader_variables[0] = 0;
            leader_variables[1] = 1;
        }
        else {
            leader_variables[0] = 1;
            leader_variables[1] = 0;
        }
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);        
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    return variable->Get(0).Length() / m_variable_epsilons[0] < variable->Get(1).Length() / m_variable_epsilons[1] ? 1 : 0;
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        d0.X = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
        d0.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
    }
    else {
        d0.X = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        d0.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        d1.X = ((WSBernsteinBasis*)m_basises[0])->CalculateDerivative(t);
        d1.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateDerivative(t);
    }
    else {
        d1.X = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        d1.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = ((WSBernsteinBasis*)m_basises[0])->CalculatePartialDerivative(variable, 0);
        dy = ((WSBernsteinBasis*)m_basises[1])->CalculatePartialDerivative(variable, 0);
    }
    else {
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        x = cache->GetTermsValue(m_term_x[curve_index]);
        y = cache->GetTermsValue(m_term_y[curve_index]);
    }
    else {
        x = -cache->GetTermsValue(m_term_x[curve_index]);
        y = -cache->GetTermsValue(m_term_y[curve_index]);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = cache->GetTermsPartialDerivative(m_term_x[curve_index], curve_index);
        dy = cache->GetTermsPartialDerivative(m_term_y[curve_index], curve_index);
    }
    else {
        dx = -cache->GetTermsPartialDerivative(m_term_x[curve_index], curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_y[curve_index], curve_index);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        const WSInterval& t = variable->Get(0);
        dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[0])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[1])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
    }
    else {
        const WSInterval& t = variable->Get(1);
        dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
    }
}

int WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::GetSampleCount() {
    return ((WSBernsteinBasis*)m_basises[0])->GetDegree() * ((WSBernsteinBasis*)m_basises[2])->GetDegree() + 1;
}

void WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem::UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, double distance_epsilon) {
    double d0 = 0;
    for (int i = 0; i < degree0; ++i) {
        d0 += (control_points0[i + 1] - control_points0[i]).Length();
    }
    if (d0 == 0) {
        m_variable_epsilons[0] = 1;
    }
    else {
        m_variable_epsilons[0] = distance_epsilon / d0;
        if (m_variable_epsilons[0] > 1) {
            m_variable_epsilons[0] = 1;
        }
    }
    double d1 = 0;
    for (int i = 0; i < degree1; ++i) {
        d1 += (control_points1[i + 1] - control_points1[i]).Length();
    }
    if (d1 == 0) {
        m_variable_epsilons[1] = 1;
    }
    else {
        m_variable_epsilons[1] = distance_epsilon / d1;
        if (m_variable_epsilons[1] > 1) {
            m_variable_epsilons[1] = 1;
        }
    }
}

WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::BezierCurveRationalBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::~BezierCurveRationalBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::Update(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    v2:     x0
    v3:     y0
    v4:     w1
    ---------------------
    b0:     x0
    b1:     y0
    b2:     x1
    b3:     y1
    b4:     w1
    b5:     v2 * v4
    b6:     v3 * v4
    b7:     v2
    b8:     v3
    b9:     v4
    ---------------------
    e0:     x0 - v2 = 0
    e1:     y0 - v3 = 0
    e2:     w1 - v4 = 0
    e3:     v2 * v4 - x1 = 0
    e4:     v3 * v4 - y1 = 0
    */
    m_distance_epsilon = distance_epsilon;
    UpdateVariableEpsilons(degree0, control_points0, degree1, control_points1, distance_epsilon);
    if (!m_initialized) {
        m_basises[0] = new WSBernsteinBasis(0, degree0);
        m_basises[1] = new WSBernsteinBasis(0, degree0);
        m_basises[2] = new WSBernsteinBasis(1, degree1);
        m_basises[3] = new WSBernsteinBasis(1, degree1);
        m_basises[4] = new WSBernsteinBasis(1, degree1);
        m_basises[5] = new WSMulBasis(2, 4);
        m_basises[6] = new WSMulBasis(3, 4);
        m_basises[7] = new WSPowerBasis(2, 1);
        m_basises[8] = new WSPowerBasis(3, 1);
        m_basises[9] = new WSPowerBasis(4, 1);
        m_equations[0] = (new WSEquation(this, 2, false, true, true))->AddTerm(0, 1)->AddTerm(7, -1);
        m_equations[1] = (new WSEquation(this, 2, false, true, true))->AddTerm(1, 1)->AddTerm(8, -1);
        m_equations[2] = (new WSEquation(this, 2, false, true, true))->AddTerm(4, 1)->AddTerm(9, -1);
        m_equations[3] = (new WSEquation(this, 2, true, true, true))->AddTerm(5, 1)->AddTerm(2, -1);
        m_equations[4] = (new WSEquation(this, 2, true, true, true))->AddTerm(6, 1)->AddTerm(3, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_term_x[0] = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_term_y[0] = term;
            }
        }
        for (int i = 0; i < m_equations[3]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[3]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_x[1] = term;
            }
        }
        for (int i = 0; i < m_equations[4]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[4]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_y[1] = term;
            }
        }
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 4) {
                m_term_w1 = term;
            }
        }
        m_initialized = true;
    }
    else {
        if (((WSBernsteinBasis*)m_basises[0])->GetDegree() != degree0) {
            delete m_basises[0];
            delete m_basises[1];
            m_basises[0] = new WSBernsteinBasis(0, degree0);
            m_basises[1] = new WSBernsteinBasis(0, degree0);
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree1) {
            delete m_basises[2];
            delete m_basises[3];
            delete m_basises[4];
            m_basises[2] = new WSBernsteinBasis(1, degree1);
            m_basises[3] = new WSBernsteinBasis(1, degree1);
            m_basises[4] = new WSBernsteinBasis(1, degree1);
        }
        ClearAlgebraRuntime();
    }
    for (int i = 0; i <= degree0; ++i) {
        ((WSBernsteinBasis*)m_basises[0])->SetCoef(i, control_points0[i].X);
        ((WSBernsteinBasis*)m_basises[1])->SetCoef(i, control_points0[i].Y);
    }
    for (int i = 0; i <= degree1; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points1[i].X * weights1[i]);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points1[i].Y * weights1[i]);
        ((WSBernsteinBasis*)m_basises[4])->SetCoef(i, weights1[i]);
    }
    m_algebra_enable = 0;
}

int WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetVariableCount() const {
    return 5;
}

WSReal WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetBasisCount() const {
    return 10;
}

WSEquationBasis* WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetEquationCount() const {
    return 5;
}

WSEquation* WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            return m_distance_epsilon * m_check_coef;
        }
    case 1: {
            return m_distance_epsilon * m_check_coef;
        }
    case 2: {
            return m_distance_epsilon * m_check_coef;
        }
    case 3: {
            WSInterval w1 = cache->GetTermsValue(m_term_w1);
            double epsilon = w1.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 4: {
            WSInterval w1 = cache->GetTermsValue(m_term_w1);
            double epsilon = w1.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

WSReal WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 1: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 2: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 3: {
            WSInterval w1 = cache->GetTermsValue(m_term_w1);
            double epsilon = w1.Max * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 4: {
            WSInterval w1 = cache->GetTermsValue(m_term_w1);
            double epsilon = w1.Max * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree0 = ((WSBernsteinBasis*)m_basises[0])->GetDegree();
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [5];
        polynomial_count = 5;

        polynomials[0] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_basises[0]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
            else if (term->GetIndex() == 7) {
                m_basises[7]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
        }

        polynomials[1] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
            else if (term->GetIndex() == 8) {
                m_basises[8]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
        }

        polynomials[2] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 4) {
                m_basises[4]->AddPolynomialTerm(polynomials[2], term->GetCoef());
            }
            else if (term->GetIndex() == 9) {
                m_basises[9]->AddPolynomialTerm(polynomials[2], term->GetCoef());
            }
        }

        polynomials[3] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[3]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[3]->GetTerm(i);
            if (term->GetIndex() == 5) {
                m_basises[5]->AddPolynomialTerm(polynomials[3], term->GetCoef());
            }
            else if (term->GetIndex() == 2) {
                m_basises[2]->AddPolynomialTerm(polynomials[3], term->GetCoef());
            }
        }
        WSPolynomialTerm* term = polynomials[3]->GetTerm(polynomials[3]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Powers[4] = 1;
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[3]->AddLast();
        
        polynomials[4] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[4]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[4]->GetTerm(i);
            if (term->GetIndex() == 6) {
                m_basises[6]->AddPolynomialTerm(polynomials[4], term->GetCoef());
            }
            else if (term->GetIndex() == 3) {
                m_basises[3]->AddPolynomialTerm(polynomials[4], term->GetCoef());
            }
        }
        term = polynomials[4]->GetTerm(polynomials[4]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Powers[4] = 1;
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[4]->AddLast();
        
        leader_variables = new int[5];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variables[2] = 2;
        leader_variables[3] = 3;
        leader_variables[4] = 4;
        leader_variable_count = 5;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
    variable->Set(2, ((WSBernsteinBasis*)m_basises[0])->CalculateValue(variable));
    variable->Set(3, ((WSBernsteinBasis*)m_basises[1])->CalculateValue(variable));
    variable->Set(4, ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable));
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    return variable->Get(0).Length() / m_variable_epsilons[0] < variable->Get(1).Length() / m_variable_epsilons[1] ? 1 : 0;
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        d0.X = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
        d0.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t) ;
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        d0.X = x / w;
        d0.Y = y / w;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        d1.X = ((WSBernsteinBasis*)m_basises[0])->CalculateDerivative(t);
        d1.Y = ((WSBernsteinBasis*)m_basises[1])->CalculateDerivative(t);
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        double dx = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        double dy = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
        double dw = ((WSBernsteinBasis*)m_basises[4])->CalculateDerivative(t);
        double w2 = w * w;
        d1.X = (dx * w - x * dw) / w2;
        d1.Y = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = ((WSBernsteinBasis*)m_basises[0])->CalculatePartialDerivative(variable, 0);
        dy = ((WSBernsteinBasis*)m_basises[1])->CalculatePartialDerivative(variable, 0);
    }
    else {
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
        WSInterval dw = ((WSBernsteinBasis*)m_basises[4])->CalculatePartialDerivative(variable, 1);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        x = cache->GetTermsValue(m_term_x[curve_index]);
        y = cache->GetTermsValue(m_term_y[curve_index]);
    }
    else {
        x = -cache->GetTermsValue(m_term_x[curve_index]);
        y = -cache->GetTermsValue(m_term_y[curve_index]);
        WSInterval w = cache->GetTermsValue(m_term_w1);
        x = x / w;
        y = y / w;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = cache->GetTermsPartialDerivative(m_term_x[curve_index], curve_index);
        dy = cache->GetTermsPartialDerivative(m_term_y[curve_index], curve_index);
    }
    else {
        WSInterval x = -cache->GetTermsValue(m_term_x[curve_index]);
        WSInterval y = -cache->GetTermsValue(m_term_y[curve_index]);
        WSInterval w = cache->GetTermsValue(m_term_w1);
        dx = -cache->GetTermsPartialDerivative(m_term_x[curve_index], curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_y[curve_index], curve_index);
        WSInterval dw = cache->GetTermsPartialDerivative(m_term_w1, curve_index);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        const WSInterval& t = variable->Get(0);
        dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[0])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[1])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
    }
    else {
        const WSInterval& t = variable->Get(1);
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        WSInterval dx = CalculateBezierD1(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy = CalculateBezierD1(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw = CalculateBezierD1(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval w2 = pow(w, 2);
        WSInterval w3 = pow(w, 3);
        WSInterval w4 = pow(w, 4);
        WSInterval w5 = pow(w, 5);
        WSInterval a4 = w4;
        WSInterval a3 = -4 * w3 * dw;
        WSInterval a2 = 6 * w2 * pow(dw, 2) - 4 * w3 * dw2;
        WSInterval a1 = 4 * w3 * dw3 - 12 * w2 * dw * dw2 + 6 * w * pow(dw, 3);
        WSInterval a0 = w4 * dw4 - 4 * w3 * dw * dw3 - 3 * pow(w2, 3) + 12 * w2 * pow(dw, 2) * dw2 - 6 * pow(dw, 4);
        dx4 = (a4 * dx4 + a3 * dx3 + a2 * dx2 + a1 * dx + a0 * x) / w5;
        dy4 = (a4 * dy4 + a3 * dy3 + a2 * dy2 + a1 * dy + a0 * y) / w5;
    }
}

int WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::GetSampleCount() {
    int degree0 = ((WSBernsteinBasis*)m_basises[0])->GetDegree();
    int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();
    return (degree0 + degree1) * degree1 + 1;
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem::UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, double distance_epsilon) {
    double d0 = 0;
    for (int i = 0; i < degree0; ++i) {
        d0 += (control_points0[i + 1] - control_points0[i]).Length();
    }
    if (d0 == 0) {
        m_variable_epsilons[0] = 1;
    }
    else {
        m_variable_epsilons[0] = distance_epsilon / d0;
        if (m_variable_epsilons[0] > 1) {
            m_variable_epsilons[0] = 1;
        }
    }
    double d1 = 0;
    for (int i = 0; i < degree1; ++i) {
        d1 += (control_points1[i + 1] - control_points1[i]).Length();
    }
    if (d1 == 0) {
        m_variable_epsilons[1] = 1;
    }
    else {
        m_variable_epsilons[1] = distance_epsilon / d1;
        if (m_variable_epsilons[1] > 1) {
            m_variable_epsilons[1] = 1;
        }
    }
}

WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::RationalBezierCurveRationalBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::~RationalBezierCurveRationalBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::Update(int degree0, const WGVector2d* control_points0, const double* weights0,
    int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    v2:     x0
    v3:     y0
    v4:     w0
    v5:     x1
    v6:     y1
    v7:     w1
    ---------------------
    b0:     x0
    b1:     y0
    b2:     w0
    b3:     x1
    b4:     y1
    b5:     w1
    b6:     v2 * v7
    b7:     v3 * v7
    b8:     v4 * v5
    b9:     v4 * v6
    b10:    v2
    b11:    v3
    b12:    v4
    b13:    v5
    b14:    v6
    b15:    v7
    ---------------------
    e0:     x0 - v2 = 0
    e1:     y0 - v3 = 0
    e2:     w0 - v4 = 0
    e3:     x1 - v5 = 0
    e4:     y1 - v6 = 0
    e5:     w1 - v7 = 0
    e6:     v2 * v7 - v4 * v5 = 0
    e7:     v3 * v7 - v4 * v6 = 0
    */
    m_distance_epsilon = distance_epsilon;
    UpdateVariableEpsilons(degree0, control_points0, degree1, control_points1, distance_epsilon);
    if (!m_initialized) {
        m_basises[0] = new WSBernsteinBasis(0, degree0);
        m_basises[1] = new WSBernsteinBasis(0, degree0);
        m_basises[2] = new WSBernsteinBasis(0, degree0);
        m_basises[3] = new WSBernsteinBasis(1, degree1);
        m_basises[4] = new WSBernsteinBasis(1, degree1);
        m_basises[5] = new WSBernsteinBasis(1, degree1);
        m_basises[6] = new WSMulBasis(2, 7);
        m_basises[7] = new WSMulBasis(3, 7);
        m_basises[8] = new WSMulBasis(4, 5);
        m_basises[9] = new WSMulBasis(4, 6);
        m_basises[10] = new WSPowerBasis(2, 1);
        m_basises[11] = new WSPowerBasis(3, 1);
        m_basises[12] = new WSPowerBasis(4, 1);
        m_basises[13] = new WSPowerBasis(5, 1);
        m_basises[14] = new WSPowerBasis(6, 1);
        m_basises[15] = new WSPowerBasis(7, 1);
        m_equations[0] = (new WSEquation(this, 2, false, true, true))->AddTerm(0, 1)->AddTerm(10, -1);
        m_equations[1] = (new WSEquation(this, 2, false, true, true))->AddTerm(1, 1)->AddTerm(11, -1);
        m_equations[2] = (new WSEquation(this, 2, false, true, true))->AddTerm(2, 1)->AddTerm(12, -1);
        m_equations[3] = (new WSEquation(this, 2, false, true, true))->AddTerm(3, 1)->AddTerm(13, -1);
        m_equations[4] = (new WSEquation(this, 2, false, true, true))->AddTerm(4, 1)->AddTerm(14, -1);
        m_equations[5] = (new WSEquation(this, 2, false, true, true))->AddTerm(5, 1)->AddTerm(15, -1);
        m_equations[6] = (new WSEquation(this, 2, true, true, true))->AddTerm(6, 1)->AddTerm(8, -1);
        m_equations[7] = (new WSEquation(this, 2, true, true, true))->AddTerm(7, 1)->AddTerm(9, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_term_x[0] = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 1) {
                m_term_y[0] = term;
            }
        }
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_w[0] = term;
            }
        }
        for (int i = 0; i < m_equations[3]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[3]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_x[1] = term;
            }
        }
        for (int i = 0; i < m_equations[4]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[4]->GetTerm(i);
            if (term->GetIndex() == 4) {
                m_term_y[1] = term;
            }
        }
        for (int i = 0; i < m_equations[5]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[5]->GetTerm(i);
            if (term->GetIndex() == 5) {
                m_term_w[1] = term;
            }
        }
        m_initialized = true;
    }
    else {
        if (((WSBernsteinBasis*)m_basises[0])->GetDegree() != degree0) {
            delete m_basises[0];
            delete m_basises[1];
            delete m_basises[2];
            m_basises[0] = new WSBernsteinBasis(0, degree0);
            m_basises[1] = new WSBernsteinBasis(0, degree0);
            m_basises[2] = new WSBernsteinBasis(0, degree0);
        }
        if (((WSBernsteinBasis*)m_basises[3])->GetDegree() != degree1) {
            delete m_basises[3];
            delete m_basises[4];
            delete m_basises[5];
            m_basises[3] = new WSBernsteinBasis(1, degree1);
            m_basises[4] = new WSBernsteinBasis(1, degree1);
            m_basises[5] = new WSBernsteinBasis(1, degree1);
        }
        ClearAlgebraRuntime();
    }
    for (int i = 0; i <= degree0; ++i) {
        ((WSBernsteinBasis*)m_basises[0])->SetCoef(i, control_points0[i].X * weights0[i]);
        ((WSBernsteinBasis*)m_basises[1])->SetCoef(i, control_points0[i].Y * weights0[i]);
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, weights0[i]);
    }
    for (int i = 0; i <= degree1; ++i) {
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points1[i].X * weights1[i]);
        ((WSBernsteinBasis*)m_basises[4])->SetCoef(i, control_points1[i].Y * weights1[i]);
        ((WSBernsteinBasis*)m_basises[5])->SetCoef(i, weights1[i]);
    }
    m_algebra_enable = 0;
}

int WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetVariableCount() const {
    return 8;
}

WSReal WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetBasisCount() const {
    return 16;
}
WSEquationBasis* WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetEquationCount() const {
    return 8;
}
WSEquation* WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            return m_distance_epsilon * m_check_coef;
        }
    case 1: {
            return m_distance_epsilon * m_check_coef;
        }
    case 3: {
            return m_distance_epsilon * m_check_coef;
        }
    case 4: {
            return m_distance_epsilon * m_check_coef;
        }
    case 6: {
            WSInterval w0 = cache->GetTermsValue(m_term_w[0]);
            WSInterval w1 = cache->GetTermsValue(m_term_w[1]);
            double epsilon = w0.Min * w1.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 7: {
            WSInterval w0 = cache->GetTermsValue(m_term_w[0]);
            WSInterval w1 = cache->GetTermsValue(m_term_w[1]);
            double epsilon = w0.Min * w1.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

WSReal WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 1: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 3: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 4: {
            return m_distance_epsilon * m_iterate_coef;
        }
    case 6: {
            WSInterval w0 = cache->GetTermsValue(m_term_w[0]);
            WSInterval w1 = cache->GetTermsValue(m_term_w[1]);
            double epsilon = w0.Min * w1.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 7: {
            WSInterval w0 = cache->GetTermsValue(m_term_w[0]);
            WSInterval w1 = cache->GetTermsValue(m_term_w[1]);
            double epsilon = w0.Min * w1.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree0 = ((WSBernsteinBasis*)m_basises[0])->GetDegree();
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [8];
        polynomial_count = 8;
        
        polynomials[0] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[0]->AddPolynomialTerm(polynomials[0], 1);
        m_basises[10]->AddPolynomialTerm(polynomials[0], -1);

        polynomials[1] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[1]->AddPolynomialTerm(polynomials[1], 1);
        m_basises[11]->AddPolynomialTerm(polynomials[1], -1);

        polynomials[2] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[2]->AddPolynomialTerm(polynomials[2], 1);
        m_basises[12]->AddPolynomialTerm(polynomials[2], -1);

        polynomials[3] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[3]->AddPolynomialTerm(polynomials[3], 1);
        m_basises[13]->AddPolynomialTerm(polynomials[3], -1);

        polynomials[4] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[4]->AddPolynomialTerm(polynomials[4], 1);
        m_basises[14]->AddPolynomialTerm(polynomials[4], -1);

        polynomials[5] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[5]->AddPolynomialTerm(polynomials[5], 1);
        m_basises[15]->AddPolynomialTerm(polynomials[5], -1);

        polynomials[6] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[6]->AddPolynomialTerm(polynomials[6], 1);
        m_basises[8]->AddPolynomialTerm(polynomials[6], -1);
        WSPolynomialTerm* term = polynomials[6]->GetTerm(polynomials[6]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Powers[4] = 1;
        term->Powers[7] = 1;
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[6]->AddLast();

        polynomials[7] = new WSPolynomial(GetVariableCount(), degree0 + 2);
        m_basises[7]->AddPolynomialTerm(polynomials[7], 1);
        m_basises[9]->AddPolynomialTerm(polynomials[7], -1);
        term = polynomials[7]->GetTerm(polynomials[7]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Powers[4] = 1;
        term->Powers[7] = 1;
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[7]->AddLast();

        leader_variables = new int[8];
        leader_variables[0] = 2;
        leader_variables[1] = 3;
        leader_variables[2] = 4;
        leader_variables[3] = 5;
        leader_variables[4] = 6;
        leader_variables[5] = 7;
        leader_variables[6] = 0;
        leader_variables[7] = 1;
        leader_variable_count = 8;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
    variable->Set(2, ((WSBernsteinBasis*)m_basises[0])->CalculateValue(variable));
    variable->Set(3, ((WSBernsteinBasis*)m_basises[1])->CalculateValue(variable));
    variable->Set(4, ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable));
    variable->Set(5, ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable));
    variable->Set(6, ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable));
    variable->Set(7, ((WSBernsteinBasis*)m_basises[5])->CalculateValue(variable));
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    return variable->Get(0).Length() / m_variable_epsilons[0] < variable->Get(1).Length() / m_variable_epsilons[1] ? 1 : 0;
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        double x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        d0.X = x / w;
        d0.Y = y / w;
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[5])->CalculateValue(t);
        d0.X = x / w;
        d0.Y = y / w;
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        double x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double dx = ((WSBernsteinBasis*)m_basises[0])->CalculateDerivative(t);
        double dy = ((WSBernsteinBasis*)m_basises[1])->CalculateDerivative(t);
        double dw = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        double w2 = w * w;
        d1.X = (dx * w - x * dw) / w2;
        d1.Y = (dy * w - y * dw) / w2;
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[5])->CalculateValue(t);
        double dx = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
        double dy = ((WSBernsteinBasis*)m_basises[4])->CalculateDerivative(t);
        double dw = ((WSBernsteinBasis*)m_basises[5])->CalculateDerivative(t);
        double w2 = w * w;
        d1.X = (dx * w - x * dw) / w2;
        d1.Y = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        WSInterval x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        dx = ((WSBernsteinBasis*)m_basises[0])->CalculatePartialDerivative(variable, 0);
        dy = ((WSBernsteinBasis*)m_basises[1])->CalculatePartialDerivative(variable, 0);
        WSInterval dw = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 0);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
    else {
        WSInterval x = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[5])->CalculateValue(variable);
        dx = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[4])->CalculatePartialDerivative(variable, 1);
        WSInterval dw = ((WSBernsteinBasis*)m_basises[5])->CalculatePartialDerivative(variable, 1);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    x = cache->GetTermsValue(m_term_x[curve_index]);
    y = cache->GetTermsValue(m_term_y[curve_index]);
    WSInterval w = cache->GetTermsValue(m_term_w[curve_index]);
    x = x / w;
    y = y / w;
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    WSInterval x = cache->GetTermsValue(m_term_x[curve_index]);
    WSInterval y = cache->GetTermsValue(m_term_y[curve_index]);
    WSInterval w = cache->GetTermsValue(m_term_w[curve_index]);
    dx = cache->GetTermsPartialDerivative(m_term_x[curve_index], curve_index);
    dy = cache->GetTermsPartialDerivative(m_term_y[curve_index], curve_index);
    WSInterval dw = cache->GetTermsPartialDerivative(m_term_w[curve_index], curve_index);
    WSInterval w2 = pow(w, 2);
    dx = (dx * w - x * dw) / w2;
    dy = (dy * w - y * dw) / w2;
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        const WSInterval& t = variable->Get(0);
        WSInterval x = ((WSBernsteinBasis*)m_basises[0])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[1])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval dx = CalculateBezierD1(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        WSInterval dy = CalculateBezierD1(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
        WSInterval dw = CalculateBezierD1(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dx2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        WSInterval dy2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
        WSInterval dw2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dx3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        WSInterval dy3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
        WSInterval dw3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[0])->GetCoefs(), t);
        WSInterval dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[1])->GetCoefs(), t);
        WSInterval dw4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval w2 = pow(w, 2);
        WSInterval w3 = pow(w, 3);
        WSInterval w4 = pow(w, 4);
        WSInterval w5 = pow(w, 5);
        WSInterval a4 = w4;
        WSInterval a3 = -4 * w3 * dw;
        WSInterval a2 = 6 * w2 * pow(dw, 2) - 4 * w3 * dw2;
        WSInterval a1 = 4 * w3 * dw3 - 12 * w2 * dw * dw2 + 6 * w * pow(dw, 3);
        WSInterval a0 = w4 * dw4 - 4 * w3 * dw * dw3 - 3 * pow(w2, 3) + 12 * w2 * pow(dw, 2) * dw2 - 6 * pow(dw, 4);
        dx4 = (a4 * dx4 + a3 * dx3 + a2 * dx2 + a1 * dx + a0 * x) / w5;
        dy4 = (a4 * dy4 + a3 * dy3 + a2 * dy2 + a1 * dy + a0 * y) / w5;
    }
    else {
        const WSInterval& t = variable->Get(1);
        WSInterval x = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[5])->CalculateValue(variable);
        WSInterval dx = CalculateBezierD1(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dy = CalculateBezierD1(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dw = CalculateBezierD1(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[5])->GetCoefs(), t);
        WSInterval dx2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dy2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dw2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[5])->GetCoefs(), t);
        WSInterval dx3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dy3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dw3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[5])->GetCoefs(), t);
        WSInterval dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dw4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[5])->GetCoefs(), t);
        WSInterval w2 = pow(w, 2);
        WSInterval w3 = pow(w, 3);
        WSInterval w4 = pow(w, 4);
        WSInterval w5 = pow(w, 5);
        WSInterval a4 = w4;
        WSInterval a3 = -4 * w3 * dw;
        WSInterval a2 = 6 * w2 * pow(dw, 2) - 4 * w3 * dw2;
        WSInterval a1 = 4 * w3 * dw3 - 12 * w2 * dw * dw2 + 6 * w * pow(dw, 3);
        WSInterval a0 = w4 * dw4 - 4 * w3 * dw * dw3 - 3 * pow(w2, 3) + 12 * w2 * pow(dw, 2) * dw2 - 6 * pow(dw, 4);
        dx4 = (a4 * dx4 + a3 * dx3 + a2 * dx2 + a1 * dx + a0 * x) / w5;
        dy4 = (a4 * dy4 + a3 * dy3 + a2 * dy2 + a1 * dy + a0 * y) / w5;
    }
}

int WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::GetSampleCount() {
    int degree0 = ((WSBernsteinBasis*)m_basises[0])->GetDegree();
    int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();
    return (degree0 + degree1) * (degree0 + degree1) + 1;
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem::UpdateVariableEpsilons(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, double distance_epsilon) {
    double d0 = 0;
    for (int i = 0; i < degree0; ++i) {
        d0 += (control_points0[i + 1] - control_points0[i]).Length();
    }
    if (d0 == 0) {
        m_variable_epsilons[0] = 1;
    }
    else {
        m_variable_epsilons[0] = distance_epsilon / d0;
        if (m_variable_epsilons[0] > 1) {
            m_variable_epsilons[0] = 1;
        }
    }
    double d1 = 0;
    for (int i = 0; i < degree1; ++i) {
        d1 += (control_points1[i + 1] - control_points1[i]).Length();
    }
    if (d1 == 0) {
        m_variable_epsilons[1] = 1;
    }
    else {
        m_variable_epsilons[1] = distance_epsilon / d1;
        if (m_variable_epsilons[1] > 1) {
            m_variable_epsilons[1] = 1;
        }
    }
}

WGIntersectHelper2d::LineBezierCurveEquationSystem::LineBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::LineBezierCurveEquationSystem::~LineBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::Update(const WGVector2d& start_point, const WGVector2d& end_point,
    int degree, const WGVector2d* control_points, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    ---------------------
    b0:     const
    b1:     v0
    b2:     x1
    b3:     y1
    ---------------------
    e0:     start_point.x + vt.x * v0 - x1 = 0
    e1:     start_point.y + vt.y * v0 - y1 = 0
    */
    m_distance_epsilon = distance_epsilon;
    WGVector2d vt = end_point - start_point;
    if (!m_initialized) {
        m_basises[0] = new WSConstBasis();
        m_basises[1] = new WSPowerBasis(0, 1);
        m_basises[2] = new WSBernsteinBasis(1, degree);
        m_basises[3] = new WSBernsteinBasis(1, degree);
        m_equations[0] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, start_point.X)->AddTerm(1, vt.X)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, start_point.Y)->AddTerm(1, vt.Y)->AddTerm(3, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_bezier_x = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_bezier_y = term;
            }
        }
        m_initialized = true;
    }
    else {
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(start_point.X);
            }
            else if (term->GetIndex() == 1) {
                term->SetCoef(vt.X);
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(start_point.Y);
            }
            else if (term->GetIndex() == 1) {
                term->SetCoef(vt.Y);
            }
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree) {
            delete m_basises[2];
            delete m_basises[3];
            m_basises[2] = new WSBernsteinBasis(1, degree);
            m_basises[3] = new WSBernsteinBasis(1, degree);
        }
        ClearAlgebraRuntime();
    }
    m_start_point = start_point;
    m_dir = vt;
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points[i].X);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points[i].Y);
    }
    m_algebra_enable = 0;
}

int WGIntersectHelper2d::LineBezierCurveEquationSystem::GetVariableCount() const {
    return 2;
}

WSReal WGIntersectHelper2d::LineBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::LineBezierCurveEquationSystem::GetBasisCount() const {
    return 4;
}

WSEquationBasis* WGIntersectHelper2d::LineBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::LineBezierCurveEquationSystem::GetEquationCount() const {
    return 2;
}

WSEquation* WGIntersectHelper2d::LineBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::LineBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_check_coef;
    }
    return g_double_epsilon;
}

WSReal WGIntersectHelper2d::LineBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_iterate_coef;
    }
    return g_double_epsilon;
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [2];
        polynomial_count = 2;

        double d = 0;
        polynomials[0] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                d = term->GetCoef();
            }
            else if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
            else if (term->GetIndex() == 2) {
                m_basises[2]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
        }
        WSPolynomialTerm* term = polynomials[0]->GetTerm(polynomials[0]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = d + WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[0]->AddLast();

        polynomials[1] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 0) {
                d = term->GetCoef();
            }
            else if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
            else if (term->GetIndex() == 3) {
                m_basises[3]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
        }
        term = polynomials[1]->GetTerm(polynomials[1]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon) + d;
        polynomials[1]->AddLast();
        
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::LineBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    return 1;
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        d0.X = m_start_point.X + m_dir.X * t;
        d0.Y = m_start_point.Y + m_dir.Y * t;
    }
    else {
        d0.X = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        d0.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        d1.X = m_dir.X;
        d1.Y = m_dir.Y;
    }
    else {
        d1.X = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        d1.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = m_dir.X;
        dy = m_dir.Y;
    }
    else {
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        WSInterval t = cache->GetVariable()->Get(0);
        x = m_start_point.X + m_dir.X * t;
        y = m_start_point.Y + m_dir.Y * t;
    }
    else {
        x = -cache->GetTermsValue(m_term_bezier_x);
        y = -cache->GetTermsValue(m_term_bezier_y);
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = m_dir.X;
        dy = m_dir.Y;
    }
    else {
        dx = -cache->GetTermsPartialDerivative(m_term_bezier_x, curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_bezier_y, curve_index);
    }
}

void WGIntersectHelper2d::LineBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        dx4 = 0;
        dy4 = 0;
    }
    else {
        const WSInterval& t = variable->Get(1);
        dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
    }
}

int WGIntersectHelper2d::LineBezierCurveEquationSystem::GetSampleCount() {
    return ((WSBernsteinBasis*)m_basises[2])->GetDegree() + 1;
}

WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::LineRationalBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::~LineRationalBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::Update(const WGVector2d& start_point, const WGVector2d& end_point,
    int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    v2:     w
    ---------------------
    b0:     v2
    b1:     v0 * v2
    b2:     x
    b3:     y
    b4:     w
    ---------------------
    e0:     start_point.x * v2 + vt.x * v0 * v2 - x = 0
    e1:     start_point.y * v2 + vt.y * v0 * v2 - y = 0
    e2:     v2 - w = 0;
    */
    m_distance_epsilon = distance_epsilon;
    WGVector2d vt = end_point - start_point;
    if (!m_initialized) {
        m_basises[0] = new WSPowerBasis(2, 1);
        m_basises[1] = new WSMulBasis(0, 2);
        m_basises[2] = new WSBernsteinBasis(1, degree);
        m_basises[3] = new WSBernsteinBasis(1, degree);
        m_basises[4] = new WSBernsteinBasis(1, degree);
        m_equations[0] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, start_point.X)->AddTerm(1, vt.X)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, start_point.Y)->AddTerm(1, vt.Y)->AddTerm(3, -1);
        m_equations[2] = (new WSEquation(this, 2, false, true, true))->AddTerm(0, 1)->AddTerm(4, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_bezier_x = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_bezier_y = term;
            }
        }
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 4) {
                m_term_bezier_w = term;
            }
        }
        m_initialized = true;
    }
    else {
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(start_point.X);
            }
            else if (term->GetIndex() == 1) {
                term->SetCoef(vt.X);
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(start_point.Y);
            }
            else if (term->GetIndex() == 1) {
                term->SetCoef(vt.Y);
            }
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree) {
            delete m_basises[2];
            delete m_basises[3];
            delete m_basises[4];
            m_basises[2] = new WSBernsteinBasis(1, degree);
            m_basises[3] = new WSBernsteinBasis(1, degree);
            m_basises[4] = new WSBernsteinBasis(1, degree);
        }
        ClearAlgebraRuntime();
    }
    m_start_point = start_point;
    m_dir = vt;
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points[i].X * weights[i]);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points[i].Y * weights[i]);
        ((WSBernsteinBasis*)m_basises[4])->SetCoef(i, weights[i]);
    }
    m_algebra_enable = 0;
}

int WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetVariableCount() const {
    return 3;
}

WSReal WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetBasisCount() const {
    return 5;
}

WSEquationBasis* WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetEquationCount() const {
    return 3;
}

WSEquation* WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 1: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

WSReal WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 1: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [3];
        polynomial_count = 3;

        polynomials[0] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_basises[0]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
            else if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
            else if (term->GetIndex() == 2) {
                m_basises[2]->AddPolynomialTerm(polynomials[0], term->GetCoef());
            }
        }
        WSPolynomialTerm* term = polynomials[0]->GetTerm(polynomials[0]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[0]->AddLast();

        polynomials[1] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_basises[0]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
            else if (term->GetIndex() == 1) {
                m_basises[1]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
            else if (term->GetIndex() == 3) {
                m_basises[3]->AddPolynomialTerm(polynomials[1], term->GetCoef());
            }
        }
        term = polynomials[1]->GetTerm(polynomials[1]->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomials[1]->AddLast();

        polynomials[2] = new WSPolynomial(GetVariableCount(), degree1 + 2);
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 0) {
                m_basises[0]->AddPolynomialTerm(polynomials[2], term->GetCoef());
            }
            else if (term->GetIndex() == 4) {
                m_basises[4]->AddPolynomialTerm(polynomials[2], term->GetCoef());
            }
        }

        leader_variables = new int[3];
        leader_variables[0] = 2;
        leader_variables[1] = 0;
        leader_variables[2] = 1;
        leader_variable_count = 3;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
    variable->Set(2, ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable));
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    return 1;
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        d0.X = m_start_point.X + m_dir.X * t;
        d0.Y = m_start_point.Y + m_dir.Y * t;
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        d0.X = x / w;
        d0.Y = y / w;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        d1.X = m_dir.X;
        d1.Y = m_dir.Y;
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        double dx = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        double dy = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
        double dw = ((WSBernsteinBasis*)m_basises[4])->CalculateDerivative(t);
        double w2 = w * w;
        d1.X = (dx * w - x * dw) / w2;
        d1.Y = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = m_dir.X;
        dy = m_dir.Y;
    }
    else {
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
        WSInterval dw = ((WSBernsteinBasis*)m_basises[4])->CalculatePartialDerivative(variable, 1);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        WSInterval t = cache->GetVariable()->Get(0);
        x = m_start_point.X + m_dir.X * t;
        y = m_start_point.Y + m_dir.Y * t;
    }
    else {
        x = -cache->GetTermsValue(m_term_bezier_x);
        y = -cache->GetTermsValue(m_term_bezier_y);
        WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
        x = x / w;
        y = y / w;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        dx = m_dir.X;
        dy = m_dir.Y;
    }
    else {
        WSInterval x = -cache->GetTermsValue(m_term_bezier_x);
        WSInterval y = -cache->GetTermsValue(m_term_bezier_y);
        WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
        dx = -cache->GetTermsPartialDerivative(m_term_bezier_x, curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_bezier_y, curve_index);
        WSInterval dw = -cache->GetTermsPartialDerivative(m_term_bezier_w, curve_index);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        dx4 = 0;
        dy4 = 0;
    }
    else {
        const WSInterval& t = variable->Get(1);
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        WSInterval dx = CalculateBezierD1(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy = CalculateBezierD1(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw = CalculateBezierD1(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval w2 = pow(w, 2);
        WSInterval w3 = pow(w, 3);
        WSInterval w4 = pow(w, 4);
        WSInterval w5 = pow(w, 5);
        WSInterval a4 = w4;
        WSInterval a3 = -4 * w3 * dw;
        WSInterval a2 = 6 * w2 * pow(dw, 2) - 4 * w3 * dw2;
        WSInterval a1 = 4 * w3 * dw3 - 12 * w2 * dw * dw2 + 6 * w * pow(dw, 3);
        WSInterval a0 = w4 * dw4 - 4 * w3 * dw * dw3 - 3 * pow(w2, 3) + 12 * w2 * pow(dw, 2) * dw2 - 6 * pow(dw, 4);
        dx4 = (a4 * dx4 + a3 * dx3 + a2 * dx2 + a1 * dx + a0 * x) / w5;
        dy4 = (a4 * dy4 + a3 * dy3 + a2 * dy2 + a1 * dy + a0 * y) / w5;
    }
}

int WGIntersectHelper2d::LineRationalBezierCurveEquationSystem::GetSampleCount() {
    return ((WSBernsteinBasis*)m_basises[2])->GetDegree() + 2;
}

WGIntersectHelper2d::ArcBezierCurveEquationSystem::ArcBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::ArcBezierCurveEquationSystem::~ArcBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::Update(const WGVector2d& center, double radius, double start_angle, double delta_angle,
    int degree, const WGVector2d* control_points, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    v2:     angle
    v3:     cos(v2) for polynomial
    v4:     sin(v2) for polynomial
    ---------------------
    b0:     const
    b1:     v0
    b2:     x
    b3:     y
    b4:     v2
    b5:     cos(v2)
    b6:     sin(v2)
    ---------------------
    e0:     center.x + radius * cos(v2) - x = 0
    e1:     center.y + radius * sin(v2) - y = 0
    e2:     start_angle + delta_angle * v0 - v2 = 0
    */
    m_distance_epsilon = distance_epsilon;
    if (!m_initialized) {
        m_basises[0] = new WSConstBasis();
        m_basises[1] = new WSPowerBasis(0, 1);
        m_basises[2] = new WSBernsteinBasis(1, degree);
        m_basises[3] = new WSBernsteinBasis(1, degree);
        m_basises[4] = new WSPowerBasis(2, 1);
        m_basises[5] = new WSCosBasis(2);
        m_basises[6] = new WSSinBasis(2);
        m_equations[0] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, center.X)->AddTerm(5, radius)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 3, true, true, true))->AddTerm(0, center.Y)->AddTerm(6, radius)->AddTerm(3, -1);
        m_equations[2] = (new WSEquation(this, 3, false, true, true))->AddTerm(0, start_angle)->AddTerm(1, delta_angle)->AddTerm(4, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_bezier_x = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_bezier_y = term;
            }
        }
        m_initialized = true;
    }
    else {
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(center.X);
            }
            else if (term->GetIndex() == 5) {
                term->SetCoef(radius);
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(center.Y);
            }
            else if (term->GetIndex() == 6) {
                term->SetCoef(radius);
            }
        }
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 0) {
                term->SetCoef(start_angle);
            }
            else if (term->GetIndex() == 1) {
                term->SetCoef(delta_angle);
            }
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree) {
            delete m_basises[2];
            delete m_basises[3];
            m_basises[2] = new WSBernsteinBasis(1, degree);
            m_basises[3] = new WSBernsteinBasis(1, degree);
        }
        ClearAlgebraRuntime();
    }
    m_center = center;
    m_radius = radius;
    m_start_angle = start_angle;
    m_delta_angle = delta_angle;
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points[i].X);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points[i].Y);
    }
    m_algebra_enable = false;
}

int WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetVariableCount() const {
    return 5;
}

WSReal WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetBasisCount() const {
    return 7;
}

WSEquationBasis* WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetEquationCount() const {
    return 3;
}

WSEquation* WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_check_coef;
    }
    return g_double_epsilon;
}

WSReal WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    if (index == 0 || index == 1) {
        return m_distance_epsilon * m_iterate_coef;
    }
    return g_double_epsilon;
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [3];
        polynomial_count = 3;

        WSPolynomial* polynomial = new WSPolynomial(GetVariableCount(), degree1 + 2);
        polynomials[0] = polynomial;
        m_basises[2]->AddPolynomialTerm(polynomial, -1);
        WSPolynomialTerm* term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = m_radius;
        term->Powers[3] = 1;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon) + m_center.X;
        polynomial->AddLast();
        
        polynomial = new WSPolynomial(GetVariableCount(), degree1 + 2);
        polynomials[1] = polynomial;
        m_basises[3]->AddPolynomialTerm(polynomial, -1);
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = m_radius;
        term->Powers[4] = 1;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon) + m_center.Y;
        polynomial->AddLast();

        polynomial = new WSPolynomial(GetVariableCount(), 3);
        polynomials[2] = polynomial;
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = 1;
        term->Powers[3] = 2;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = 1;
        term->Powers[4] = 2;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = -1;
        polynomial->AddLast();

        leader_variables = new int[2];
        leader_variables[0] = 3;
        leader_variables[1] = 4;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
    WSInterval a = m_start_angle + variable->Get(0) * m_delta_angle;
    variable->Set(2, a);
    variable->Set(3, cos(a));
    variable->Set(4, sin(a));
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    if (variable->Get(4).Length() > g_pi) {
        return 0;
    }
    return 1;
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        double a = m_start_angle + t * m_delta_angle;
        d0.X = m_center.X + m_radius * cos(a);
        d0.Y = m_center.Y + m_radius * sin(a);
    }
    else {
        d0.X = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        d0.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        double a = m_start_angle + t * d;
        d1.X = -m_radius * d * sin(a);
        d1.Y = m_radius * d * cos(a);
    }
    else {
        d1.X = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        d1.Y = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + variable->Get(0) * d;
        dx = -m_radius * d * sin(a);
        dy = m_radius * d * cos(a);
    }
    else {
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        WSInterval t = cache->GetVariable()->Get(0);
        WSInterval a = m_start_angle + t * m_delta_angle;
        x = m_center.X + m_radius * cos(a);
        y = m_center.Y + m_radius * sin(a);
    }
    else {
        x = -cache->GetTermsValue(m_term_bezier_x);
        y = -cache->GetTermsValue(m_term_bezier_y);
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + cache->GetVariable()->Get(0) * d;
        dx = -m_radius * d * sin(a);
        dy = m_radius * d * cos(a);
    }
    else {
        dx = -cache->GetTermsPartialDerivative(m_term_bezier_x, curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_bezier_y, curve_index);
    }
}

void WGIntersectHelper2d::ArcBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + variable->Get(0) * d;
        double d4 = d * d * d * d;
        dx4 = m_radius * d4 * cos(a);
        dy4 = m_radius * d4 * sin(a);
    }
    else {
        const WSInterval& t = variable->Get(1);
        dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
    }
}

int WGIntersectHelper2d::ArcBezierCurveEquationSystem::GetSampleCount() {
    return ((WSBernsteinBasis*)m_basises[2])->GetDegree() * 2 + 1;
}

WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::ArcRationalBezierCurveEquationSystem() :
    m_initialized(false),
    m_check_coef(1),
    m_iterate_coef(1) {
}

WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::~ArcRationalBezierCurveEquationSystem() {
    if (m_initialized) {
        for (int i = GetBasisCount() - 1; i >= 0; --i) {
            delete m_basises[i];
        }
        for (int i = GetEquationCount() - 1; i >= 0; --i) {
            delete m_equations[i];
        }
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::Update(const WGVector2d& center, double radius, double start_angle, double delta_angle,
    int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon) {
    /*
    v0:     t0
    v1:     t1
    v2:     angle
    v3:     cos(v2)
    v4:     sin(v2)
    v5:     w
    ---------------------
    b0:     cos(v2)
    b1:     sin(v2)
    b2:     x
    b3:     y
    b4:     w
    b5:     v3 * v5
    b6:     v4 * v5
    b7:     v0
    b8:     v2
    b9:     const
    b10:    v3
    b11:    v4
    b12:    v5
    ---------------------
    e0:     center.x * w + radius * v3 * v5 - x = 0
    e1:     center.y * w + radius * v4 * v5 - y = 0
    e2:     start_angle + delta_angle * v0 - v2 = 0
    e3:     v3 - cos(v2) = 0
    e4:     v4 - sin(v2) = 0
    e5:     v5 - w = 0
    */
    m_distance_epsilon = distance_epsilon;
    if (!m_initialized) {
        m_basises[0] = new WSCosBasis(2);
        m_basises[1] = new WSSinBasis(2);
        m_basises[2] = new WSBernsteinBasis(1, degree);
        m_basises[3] = new WSBernsteinBasis(1, degree);
        m_basises[4] = new WSBernsteinBasis(1, degree);
        m_basises[5] = new WSMulBasis(3, 5);
        m_basises[6] = new WSMulBasis(4, 5);
        m_basises[7] = new WSPowerBasis(0, 1);
        m_basises[8] = new WSPowerBasis(2, 1);
        m_basises[9] = new WSConstBasis();
        m_basises[10] = new WSPowerBasis(3, 1);
        m_basises[11] = new WSPowerBasis(4, 1);
        m_basises[12] = new WSPowerBasis(5, 1);
        m_equations[0] = (new WSEquation(this, 3, true, true, true))->AddTerm(4, center.X)->AddTerm(5, radius)->AddTerm(2, -1);
        m_equations[1] = (new WSEquation(this, 3, true, true, true))->AddTerm(4, center.Y)->AddTerm(6, radius)->AddTerm(3, -1);
        m_equations[2] = (new WSEquation(this, 3, false, true, true))->AddTerm(9, start_angle)->AddTerm(7, delta_angle)->AddTerm(8, -1);
        m_equations[3] = (new WSEquation(this, 2, false, true, true))->AddTerm(10, 1)->AddTerm(0, -1);
        m_equations[4] = (new WSEquation(this, 2, false, true, true))->AddTerm(11, 1)->AddTerm(1, -1);
        m_equations[5] = (new WSEquation(this, 2, false, true, true))->AddTerm(12, 1)->AddTerm(4, -1);
        BuildRuntime();
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 2) {
                m_term_bezier_x = term;
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 3) {
                m_term_bezier_y = term;
            }
        }
        for (int i = 0; i < m_equations[5]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[5]->GetTerm(i);
            if (term->GetIndex() == 4) {
                m_term_bezier_w = term;
            }
        }
        m_initialized = true;
    }
    else {
        for (int i = 0; i < m_equations[0]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[0]->GetTerm(i);
            if (term->GetIndex() == 4) {
                term->SetCoef(center.X);
            }
            else if (term->GetIndex() == 5) {
                term->SetCoef(radius);
            }
        }
        for (int i = 0; i < m_equations[1]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[1]->GetTerm(i);
            if (term->GetIndex() == 4) {
                term->SetCoef(center.Y);
            }
            else if (term->GetIndex() == 6) {
                term->SetCoef(radius);
            }
        }
        for (int i = 0; i < m_equations[2]->GetTermCount(); ++i) {
            WSTerm* term = m_equations[2]->GetTerm(i);
            if (term->GetIndex() == 9) {
                term->SetCoef(start_angle);
            }
            else if (term->GetIndex() == 7) {
                term->SetCoef(delta_angle);
            }
        }
        if (((WSBernsteinBasis*)m_basises[2])->GetDegree() != degree) {
            delete m_basises[2];
            delete m_basises[3];
            delete m_basises[4];
            m_basises[2] = new WSBernsteinBasis(1, degree);
            m_basises[3] = new WSBernsteinBasis(1, degree);
            m_basises[4] = new WSBernsteinBasis(1, degree);
        }
        ClearAlgebraRuntime();
    }
    m_center = center;
    m_radius = radius;
    m_start_angle = start_angle;
    m_delta_angle = delta_angle;
    for (int i = 0; i <= degree; ++i) {
        ((WSBernsteinBasis*)m_basises[2])->SetCoef(i, control_points[i].X * weights[i]);
        ((WSBernsteinBasis*)m_basises[3])->SetCoef(i, control_points[i].Y * weights[i]);
        ((WSBernsteinBasis*)m_basises[4])->SetCoef(i, weights[i]);
    }
    m_algebra_enable = false;
}

int WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetVariableCount() const {
    return 6;
}

WSReal WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetVariableEpsilon(const WSIntervalVector* variable, int index) const {
    return g_double_epsilon;
}

int WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetBasisCount() const {
    return 13;
}
WSEquationBasis* WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetBasis(int index) const {
    return m_basises[index];
}

int WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetEquationCount() const {
    return 6;
}
WSEquation* WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetEquation(int index) const {
    return m_equations[index];
}

WSReal WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 1: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_check_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

WSReal WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const {
    switch (index) {
    case 0: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    case 1: {
            WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
            double epsilon = w.Min * m_distance_epsilon * m_iterate_coef;
            if (epsilon < g_double_epsilon) {
                epsilon = g_double_epsilon;
            }
            return epsilon;
        }
    default: {
            return g_double_epsilon;
        }
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain,
    WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const {
    if (m_algebra_enable == 1) {
        int degree1 = ((WSBernsteinBasis*)m_basises[2])->GetDegree();

        polynomials = new WSPolynomial * [4];
        polynomial_count = 4;

        WSPolynomial* polynomial = new WSPolynomial(GetVariableCount(), degree1 + 2);
        polynomials[0] = polynomial;
        m_basises[4]->AddPolynomialTerm(polynomial, m_center.X);
        m_basises[5]->AddPolynomialTerm(polynomial, m_radius);
        m_basises[2]->AddPolynomialTerm(polynomial, -1);
        WSPolynomialTerm* term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomial->AddLast();

        polynomial = new WSPolynomial(GetVariableCount(), degree1 + 2);
        polynomials[1] = polynomial;
        m_basises[4]->AddPolynomialTerm(polynomial, m_center.Y);
        m_basises[6]->AddPolynomialTerm(polynomial, m_radius);
        m_basises[3]->AddPolynomialTerm(polynomial, -1);
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = WSInterval(-m_algebra_epsilon, m_algebra_epsilon);
        polynomial->AddLast();

        polynomial = new WSPolynomial(GetVariableCount(), 3);
        polynomials[2] = polynomial;
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = 1;
        term->Powers[3] = 2;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = 1;
        term->Powers[4] = 2;
        polynomial->AddLast();
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = -1;
        polynomial->AddLast();

        polynomial = new WSPolynomial(GetVariableCount(), degree1 + 2);
        polynomials[3] = polynomial;
        m_basises[4]->AddPolynomialTerm(polynomial, -1);
        term = polynomial->GetTerm(polynomial->NewTerm());
        memset(term->Powers, 0, GetVariableCount() * sizeof(int));
        term->Coef = 1;
        term->Powers[5] = 1;
        polynomial->AddLast();

        leader_variables = new int[3];
        leader_variables[0] = 5;
        leader_variables[1] = 3;
        leader_variables[2] = 4;
        leader_variable_count = 3;
        max_term_count = 256;
    }
    else if (m_algebra_enable == 3) {
        BuildFittingPolynomials(variable_domain, polynomials, polynomial_count);
        leader_variables = new int[2];
        leader_variables[0] = 0;
        leader_variables[1] = 1;
        leader_variable_count = 2;
        max_term_count = 256;
    }
    else {
        polynomials = nullptr;
        polynomial_count = 0;
        leader_variables = nullptr;
        leader_variable_count = 0;
        max_term_count = 0;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::InitializeIntermediateVariables(WSIntervalVector* variable) {
    WSInterval a = m_start_angle + variable->Get(0) * m_delta_angle;
    variable->Set(2, a);
    variable->Set(3, cos(a));
    variable->Set(4, sin(a));
    variable->Set(5, ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable));
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::SetEpsilonCoef(double check_coef, double iterate_coef) {
    m_check_coef = check_coef;
    m_iterate_coef = iterate_coef;
}

int WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetSplitIndex(const WSIntervalVector* variable) {
    if (variable->Get(4).Length() > g_pi) {
        return 0;
    }
    return 1;
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD0(int curve_index, double t, WGVector2d& d0) const {
    if (curve_index == 0) {
        double a = m_start_angle + t * m_delta_angle;
        d0.X = m_center.X + m_radius * cos(a);
        d0.Y = m_center.Y + m_radius * sin(a);
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        d0.X = x / w;
        d0.Y = y / w;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD1(int curve_index, double t, WGVector2d& d1) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        double a = m_start_angle + t * d;
        d1.X = -m_radius * d * sin(a);
        d1.Y = m_radius * d * cos(a);
    }
    else {
        double x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(t);
        double y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(t);
        double w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(t);
        double dx = ((WSBernsteinBasis*)m_basises[2])->CalculateDerivative(t);
        double dy = ((WSBernsteinBasis*)m_basises[3])->CalculateDerivative(t);
        double dw = ((WSBernsteinBasis*)m_basises[4])->CalculateDerivative(t);
        double w2 = w * w;
        d1.X = (dx * w - x * dw) / w2;
        d1.Y = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD1(int curve_index, const WSIntervalVector* variable, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + variable->Get(0) * d;
        dx = -m_radius * d * sin(a);
        dy = m_radius * d * cos(a);
    }
    else {
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        dx = ((WSBernsteinBasis*)m_basises[2])->CalculatePartialDerivative(variable, 1);
        dy = ((WSBernsteinBasis*)m_basises[3])->CalculatePartialDerivative(variable, 1);
        WSInterval dw = ((WSBernsteinBasis*)m_basises[4])->CalculatePartialDerivative(variable, 1);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD0(int curve_index, WSEquationsCache* cache, WSInterval& x, WSInterval& y) const {
    if (curve_index == 0) {
        WSInterval t = cache->GetVariable()->Get(0);
        WSInterval a = m_start_angle + t * m_delta_angle;
        x = m_center.X + m_radius * cos(a);
        y = m_center.Y + m_radius * sin(a);
    }
    else {
        x = -cache->GetTermsValue(m_term_bezier_x);
        y = -cache->GetTermsValue(m_term_bezier_y);
        WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
        x = x / w;
        y = y / w;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD1(int curve_index, WSEquationsCache* cache, WSInterval& dx, WSInterval& dy) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + cache->GetVariable()->Get(0) * d;
        dx = -m_radius * d * sin(a);
        dy = m_radius * d * cos(a);
    }
    else {
        WSInterval x = -cache->GetTermsValue(m_term_bezier_x);
        WSInterval y = -cache->GetTermsValue(m_term_bezier_y);
        WSInterval w = -cache->GetTermsValue(m_term_bezier_w);
        dx = -cache->GetTermsPartialDerivative(m_term_bezier_x, curve_index);
        dy = -cache->GetTermsPartialDerivative(m_term_bezier_y, curve_index);
        WSInterval dw = -cache->GetTermsPartialDerivative(m_term_bezier_w, curve_index);
        WSInterval w2 = pow(w, 2);
        dx = (dx * w - x * dw) / w2;
        dy = (dy * w - y * dw) / w2;
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::CalculateD4(int curve_index, const WSIntervalVector* variable, WSInterval& dx4, WSInterval& dy4) const {
    if (curve_index == 0) {
        double d = m_delta_angle;
        WSInterval a = m_start_angle + variable->Get(0) * d;
        double d4 = d * d * d * d;
        dx4 = m_radius * d4 * cos(a);
        dy4 = m_radius * d4 * sin(a);
    }
    else {
        const WSInterval& t = variable->Get(1);
        WSInterval x = ((WSBernsteinBasis*)m_basises[2])->CalculateValue(variable);
        WSInterval y = ((WSBernsteinBasis*)m_basises[3])->CalculateValue(variable);
        WSInterval w = ((WSBernsteinBasis*)m_basises[4])->CalculateValue(variable);
        WSInterval dx = CalculateBezierD1(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy = CalculateBezierD1(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw = CalculateBezierD1(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw2 = CalculateBezierD2(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw3 = CalculateBezierD3(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval dx4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[2])->GetDegree(), ((WSBernsteinBasis*)m_basises[2])->GetCoefs(), t);
        WSInterval dy4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[3])->GetDegree(), ((WSBernsteinBasis*)m_basises[3])->GetCoefs(), t);
        WSInterval dw4 = CalculateBezierD4(((WSBernsteinBasis*)m_basises[4])->GetDegree(), ((WSBernsteinBasis*)m_basises[4])->GetCoefs(), t);
        WSInterval w2 = pow(w, 2);
        WSInterval w3 = pow(w, 3);
        WSInterval w4 = pow(w, 4);
        WSInterval w5 = pow(w, 5);
        WSInterval a4 = w4;
        WSInterval a3 = -4 * w3 * dw;
        WSInterval a2 = 6 * w2 * pow(dw, 2) - 4 * w3 * dw2;
        WSInterval a1 = 4 * w3 * dw3 - 12 * w2 * dw * dw2 + 6 * w * pow(dw, 3);
        WSInterval a0 = w4 * dw4 - 4 * w3 * dw * dw3 - 3 * pow(w2, 3) + 12 * w2 * pow(dw, 2) * dw2 - 6 * pow(dw, 4);
        dx4 = (a4 * dx4 + a3 * dx3 + a2 * dx2 + a1 * dx + a0 * x) / w5;
        dy4 = (a4 * dy4 + a3 * dy3 + a2 * dy2 + a1 * dy + a0 * y) / w5;
    }
}

int WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem::GetSampleCount() {
    return ((WSBernsteinBasis*)m_basises[2])->GetDegree() * (((WSBernsteinBasis*)m_basises[2])->GetDegree() + 1) + 1;
}

int WGIntersectHelper2d::LineLineIntersect(const WGVector2d& start_point0, const WGVector2d& end_point0,
    const WGVector2d& start_point1, const WGVector2d& end_point1, double distance_epsilon,
    bool& is_singularity0, bool& is_singularity1, CurveCurveIntersection intersections[2]) {
    WGVector2d vt0 = end_point0 - start_point0;
    WGVector2d vt1 = end_point1 - start_point1;
    if (is_zero(vt0.X, g_double_epsilon) && is_zero(vt0.Y, g_double_epsilon)) {
        if (is_zero(vt1.X, g_double_epsilon) && is_zero(vt1.Y, g_double_epsilon)) {
            is_singularity0 = true;
            is_singularity1 = true;
            return 0;
        }
        else {
            is_singularity0 = true;
            is_singularity1 = false;
            return 0;
        }
    }
    else {
        if (is_zero(vt1.X, g_double_epsilon) && is_zero(vt1.Y, g_double_epsilon)) {
            is_singularity0 = false;
            is_singularity1 = true;
            return 0;
        }
        else {
            is_singularity0 = false;
            is_singularity1 = false;
            double dot01 = vt0.Dot(vt1);
            if (dot01 == 0) {
                double cross01 = vt0.Cross(vt1);
                WGVector2d vt = start_point1 - start_point0;
                double t0 = vt.Cross(vt1) / cross01;
                if (t0 >= 0 && t0 <= 1) {
                    double t1 = vt.Cross(vt0) / cross01;
                    if (t1 >= 0 && t1 <= 1) {
                        intersections[0] = CurveCurveIntersection(t0, t1, false);
                        return 1;
                    }
                }
                return 0;
            }
            else {
                double epsilon2 = distance_epsilon * distance_epsilon;
                double ts0[5];
                double ts1[5];
                bool is_samples[5];
                int tc = 0;
                double d0 = vt0.Dot(vt0);
                double d1 = vt1.Dot(vt1);
                WGVector2d start_vt01 = start_point1 - start_point0;
                WGVector2d end_vt01 = end_point1 - end_point0;
                ts0[tc] = start_vt01.Dot(vt0) / d0;
                if (ts0[tc] < 0) {
                    ts0[tc] = 0;
                }
                if (ts0[tc] > 1) {
                    ts0[tc] = 1;
                }
                if ((start_vt01 - vt0 * ts0[tc]).SqrLength() <= epsilon2) {
                    ts1[tc] = 0;
                    is_samples[tc] = true;
                    ++tc;
                }
                ts1[tc] = -start_vt01.Dot(vt1) / d1;
                if (ts1[tc] < 0) {
                    ts1[tc] = 0;
                }
                if (ts1[tc] > 1) {
                    ts1[tc] = 1;
                }
                if ((start_vt01 + vt1 * ts1[tc]).SqrLength() <= epsilon2) {
                    ts0[tc] = 0;
                    is_samples[tc] = true;
                    ++tc;
                }                
                ts0[tc] = 1 + end_vt01.Dot(vt0) / d0;
                if (ts0[tc] < 0) {
                    ts0[tc] = 0;
                }
                if (ts0[tc] > 1) {
                    ts0[tc] = 1;
                }
                if ((end_vt01 + vt0 * (1 - ts0[tc])).SqrLength() <= epsilon2) {
                    ts1[tc] = 1;
                    is_samples[tc] = true;
                    ++tc;
                }
                ts1[tc] = 1 - end_vt01.Dot(vt1) / d1;
                if (ts1[tc] < 0) {
                    ts1[tc] = 0;
                }
                if (ts1[tc] > 1) {
                    ts1[tc] = 1;
                }
                if ((end_vt01 - vt1 * (1 - ts1[tc])).SqrLength() <= epsilon2) {
                    ts0[tc] = 1;
                    is_samples[tc] = true;
                    ++tc;
                }
                double cross01 = vt0.Cross(vt1);
                if (cross01 != 0) {
                    ts0[tc] = start_vt01.Cross(vt1) / cross01;
                    ts1[tc] = start_vt01.Cross(vt0) / cross01;
                    if (ts0[tc] >= 0 && ts0[tc] <= 1 && ts1[tc] >= 0 && ts1[tc] <= 1) {
                        is_samples[tc] = false;
                        ++tc;
                    }
                }
                int i = 0;
                while (i < tc - 1) {
                    bool b = true;
                    int j = i + 1;
                    while (j < tc) {
                        if ((ts0[j] - ts0[i]) * (ts1[j] - ts1[i]) * dot01 < 0) {
                            double t = ts1[i];
                            ts1[i] = ts1[j];
                            ts1[j] = t;
                            if ((start_vt01 + vt1 * ts1[j] - vt0 * ts0[j]).SqrLength() > epsilon2) {
                                --tc;
                                if (tc > j) {
                                    ts0[j] = ts0[tc];
                                    ts1[j] = ts1[tc];
                                    is_samples[j] = is_samples[tc];
                                } 
                            }
                            if ((start_vt01 + vt1 * ts1[i] - vt0 * ts0[i]).SqrLength() > epsilon2) {
                                --tc;
                                if (tc > i) {
                                    ts0[i] = ts0[tc];
                                    ts1[i] = ts1[tc];
                                    is_samples[i] = is_samples[tc];
                                }
                            }
                            b = false;
                            break;
                        }
                        ++j;
                    }
                    if (b) {
                        ++i;
                    }
                }
                if (tc == 0) {
                    return 0;
                }
                if (tc == 1) {
                    intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
                    return 1;
                }
                double t00 = ts0[0];
                double t01 = ts0[0];
                double t10 = ts1[0];
                double t11 = ts1[0];
                for (int i = 0; i < tc; ++i) {
                    if (ts0[i] < t00) {
                        t00 = ts0[i];
                    }
                    else if (ts0[i] > t01) {
                        t01 = ts0[i];
                    }
                    if (ts1[i] < t10) {
                        t10 = ts1[i];
                    }
                    else if (ts1[i] > t11) {
                        t11 = ts1[i];
                    }
                }
                bool b = false;
                double t0, t1;
                for (int i = 0; i < tc; ++i) {
                    if (!is_samples[i]) {
                        t0 = ts0[i];
                        t1 = ts1[i];
                        b = true;
                    }
                }
                if (dot01 < 0) {
                    double t = t10;
                    t10 = t11;
                    t11 = t;
                }
                int n = 0;
                if (b) {
                    if (t00 != t0 || t10 != t1) {
                        intersections[n] = CurveCurveIntersection(t00, t10, true, t0, t1, false);
                        ++n;
                    }
                    if (t01 != t0 || t11 != t1) {
                        intersections[n] = CurveCurveIntersection(t0, t1, false, t01, t11, true);
                        ++n;
                    }
                    if (n == 0) {
                        intersections[n] = CurveCurveIntersection(t0, t1, false);
                        ++n;
                    }
                }
                else {
                    if (t00 != t01 || t10 != t11) {
                        intersections[n] = CurveCurveIntersection(t00, t10, true, t01, t11, true);
                        ++n;
                    }
                    if (n == 0) {
                        intersections[n] = CurveCurveIntersection(t00, t10, true);
                        ++n;
                    }
                }
                return n;
            }
        }
    }
}

int WGIntersectHelper2d::LineArcIntersect(const WGVector2d& line_start_point, const WGVector2d& line_end_point,
    const WGVector2d& arc_center, double arc_radius, double arc_start_angle, double arc_delta_angle, double distance_epsilon,
    bool& is_line_singularity, bool& is_arc_singularity, CurveCurveIntersection intersections[10]) {
    WGVector2d line_vt = line_end_point - line_start_point;
    is_line_singularity = is_zero(line_vt.X, g_double_epsilon) && is_zero(line_vt.Y, g_double_epsilon);
    is_arc_singularity = arc_radius <= 0 || is_zero(arc_delta_angle, g_double_epsilon);
    if (is_line_singularity || is_arc_singularity) {
        return 0;
    }
    double line_length = line_vt.Length();
    WGVector2d line_unit_vt = line_vt / line_length;
    WGVector2d vt0 = arc_center - line_start_point;
    double d = vt0.Cross(line_unit_vt);
    double center_line_distance = abs(d);
    if (center_line_distance > arc_radius + distance_epsilon) {
        return 0;
    }
    double ts0[10];
    double ts1[10];
    bool is_samples[10];
    WGVector2d vt(-line_unit_vt.Y, line_unit_vt.X);
    double vt_angle = vt.UnitAngle();
    double extreme_angle1 = vt_angle;
    double extreme_angle2 = extreme_angle1 + g_pi;
    if (extreme_angle2 >= g_pi * 2) {
        extreme_angle2 -= g_pi * 2;
    }
    if (d < 0) {
        double t = extreme_angle1;
        extreme_angle1 = extreme_angle2;
        extreme_angle2 = t;
    }
    bool is_extreme_ok1 = center_line_distance >= arc_radius - distance_epsilon;
    int tc = 0;
    if (center_line_distance >= arc_radius) {
        WGVector2d pt = arc_center + vt * d;
        double t0 = (pt - line_start_point).Dot(line_unit_vt) / line_length;
        if (t0 >= 0 && t0 <= 1) {
            double t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, vt_angle);
            if (t1 >= 0 && t1 <= 1) {
                ts0[tc] = t0;
                ts1[tc] = t1;
                is_samples[tc] = false;
                ++tc;
            }
        }
    }
    else {
        double d1 = vt0.Dot(line_unit_vt);
        double d2 = sqrt(arc_radius * arc_radius - d * d);
        {
            double c = d1 - d2;
            double t0 = c / line_length;
            if (t0 >= 0 && t0 <= 1) {
                WGVector2d vt = (line_unit_vt * c - vt0) / arc_radius;
                double t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, vt.UnitAngle());
                if (t1 >= 0 && t1 <= 1) {
                    ts0[tc] = t0;
                    ts1[tc] = t1;
                    is_samples[tc] = false;
                    ++tc;
                }
            }
        }
        {
            double c = d1 + d2;
            double t0 = c / line_length;
            if (t0 >= 0 && t0 <= 1) {
                WGVector2d vt = (line_unit_vt * c - vt0) / arc_radius;
                double t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, vt.UnitAngle());
                if (t1 >= 0 && t1 <= 1) {
                    ts0[tc] = t0;
                    ts1[tc] = t1;
                    is_samples[tc] = false;
                    ++tc;
                }
            }
        }
    }
    int tc1 = tc;
    double epsilon2 = distance_epsilon * distance_epsilon;
    WGVector2d arc_start_point = WGArc2d::CalculatePoint(arc_center, arc_radius, arc_start_angle);
    WGVector2d arc_end_point = WGArc2d::CalculatePoint(arc_center, arc_radius, arc_start_angle + arc_delta_angle);
    vt = arc_start_point - line_start_point;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 0;
        ts1[tc] = 0;
        is_samples[tc] = true;
        ++tc;
    }
    if (abs(vt.Cross(line_unit_vt)) <= distance_epsilon) {
        double t0 = vt.Dot(line_unit_vt) / line_length;
        if (t0 >= 0 && t0 <= 1) {
            ts0[tc] = t0;
            ts1[tc] = 0;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = arc_end_point - line_start_point;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 0;
        ts1[tc] = 1;
        is_samples[tc] = true;
        ++tc;
    }
    if (abs(vt.Cross(line_unit_vt)) <= distance_epsilon) {
        double t0 = vt.Dot(line_unit_vt) / line_length;
        if (t0 >= 0 && t0 <= 1) {
            ts0[tc] = t0;
            ts1[tc] = 1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = -vt0;
    d = vt.Normalize(g_double_epsilon);
    if (d > g_double_epsilon && abs(d - arc_radius) <= distance_epsilon) {
        double t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, vt.UnitAngle());
        if (t1 >= 0 && t1 <= 1) {
            ts0[tc] = 0;
            ts1[tc] = t1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = line_end_point - arc_center;
    d = vt.Normalize(g_double_epsilon);
    if (d > g_double_epsilon && abs(d - arc_radius) <= distance_epsilon) {
        double t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, vt.UnitAngle());
        if (t1 >= 0 && t1 <= 1) {
            ts0[tc] = 1;
            ts1[tc] = t1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = arc_start_point - line_end_point;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 1;
        ts1[tc] = 0;
        is_samples[tc] = true;
        ++tc;
    }
    vt = arc_end_point - line_end_point;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 1;
        ts1[tc] = 1;
        is_samples[tc] = true;
        ++tc;
    }
    if (tc == 0) {
        return 0;
    }
    if (tc == tc1) {
        if (tc == 1) {
            intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
            return 1;
        }
        else {
            double extreme_t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, extreme_angle1);
            if ((extreme_t1 < ts1[0] && extreme_t1 > ts1[1]) || (extreme_t1 < ts1[1] && extreme_t1 > ts1[0])) {
                if (is_extreme_ok1) {
                    intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0], ts0[1], ts1[1], is_samples[1]);
                    return 1;
                }
            }
            intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
            intersections[1] = CurveCurveIntersection(ts0[1], ts1[1], is_samples[1]);
            return 2;
        }
    }
    RemoveSameTs(ts0, ts1, is_samples, tc);
    if (tc == 1) {
        intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
        return 1;
    }
    WGVector2d points0[10];
    WGVector2d points1[10];
    for (int i = 0; i < tc; ++i) {
        points0[i] = line_start_point + line_vt * ts0[i];
        points1[i] = WGArc2d::CalculatePoint(arc_center, arc_radius, arc_start_angle + ts1[i] * arc_delta_angle);
    }
    RepairCross(ts0, ts1, is_samples, points0, points1, tc, distance_epsilon);
    double extreme_t1 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, extreme_angle1);
    double extreme_t2 = WGArc2d::CalculateT(arc_start_angle, arc_delta_angle, extreme_angle2);
    int indices[10];
    for (int i = 0; i < tc; ++i) {
        indices[i] = i;
    }
    int count = 0;
    SortTs(ts0, ts1, indices, tc);
    for (int j = 1; j < tc; ++j) {
        int i = j - 1;
        double t00 = ts0[indices[i]];
        double t01 = ts0[indices[j]];
        double t10 = ts1[indices[i]];
        double t11 = ts1[indices[j]];
        if (t00 != t01 || t10 != t11) {
            bool b = false;
            if (abs(arc_delta_angle * (t11 - t10)) < g_pi) {
                if ((extreme_t2 > t10 || extreme_t2 < t11) && (extreme_t2 > t11 || extreme_t2 < t10)) {
                    if ((extreme_t1 < t10 && extreme_t1 > t11) || (extreme_t1 < t11 && extreme_t1 > t10)) {
                        if (is_extreme_ok1) {
                            b = true;
                        }
                    }
                    else {
                        b = true;
                    }
                }
            }
            if (b) {
                intersections[count] = CurveCurveIntersection(t00, t10, is_samples[indices[i]], t01, t11, is_samples[indices[j]]);
                ++count;
            }
        }
    }
    for (int i = 0; i < tc; ++i) {
        bool b = true;
        for (int j = 0; j < count; ++j) {
            const CurveCurveIntersection& intersection = intersections[j];
            if (intersection.PointCount == 1) {
                if (intersection.Ts0[0] != ts0[i] || intersection.Ts1[0] != ts1[i]) {
                    continue;
                }
            }
            else {
                if (ts0[i] < intersection.Ts0[0] || ts0[i] > intersection.Ts0[1]) {
                    continue;
                }
                if (intersection.Ts1[0] <= intersection.Ts1[1]) {
                    if (ts1[i] < intersection.Ts1[0] || ts1[i] > intersection.Ts1[1]) {
                        continue;
                    }
                }
                else {
                    if (ts1[i] < intersection.Ts1[1] || ts1[i] > intersection.Ts1[0]) {
                        continue;
                    }
                }
            }
            b = false;
            break;
        }
        if (b) {
            intersections[count] = CurveCurveIntersection(ts0[i], ts1[i], is_samples[i]);
            ++count;
        }
    }
    MarkNotSample(ts0, ts1, is_samples, tc, intersections, count);
    count = MergeOverlapIntersections(intersections, count);
    return count;
}

int WGIntersectHelper2d::ArcArcIntersect(const WGVector2d& center0, double radius0, double start_angle0, double delta_angle0,
    const WGVector2d& center1, double radius1, double start_angle1, double delta_angle1, double distance_epsilon,
    bool& is_singularity0, bool& is_singularity1, CurveCurveIntersection intersections[10]) {
    is_singularity0 = radius0 <= g_double_epsilon || is_zero(delta_angle0, g_double_epsilon);
    is_singularity1 = radius1 <= g_double_epsilon || is_zero(delta_angle1, g_double_epsilon);
    if (is_singularity0 || is_singularity1) {
        return 0;
    }
    WGVector2d center_center_vt = center1 - center0;
    double center_center_distance = center_center_vt.Length();
    if (abs(radius1 - radius0) > center_center_distance + distance_epsilon) {
        return 0;
    }
    if (radius1 + radius0 < center_center_distance - distance_epsilon) {
        return 0;
    }
    double ts0[10];
    double ts1[10];
    bool is_samples[10];
    int tc = 0;
    if (center_center_distance != 0) {
        if (radius1 - radius0 >= center_center_distance) {
            WGVector2d vt = -center_center_vt / center_center_distance;
            double angle = vt.UnitAngle();
            double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, angle);
            if (t0 >= 0 && t0 <= 1) {
                double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, angle);
                if (t1 >= 0 && t1 <= 1) {
                    ts0[tc] = t0;
                    ts1[tc] = t1;
                    is_samples[tc] = false;
                    ++tc;
                }
            }
        }
        else if (radius0 - radius1 >= center_center_distance) {
            WGVector2d vt = center_center_vt / center_center_distance;
            double angle = vt.UnitAngle();
            double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, angle);
            if (t0 >= 0 && t0 <= 1) {
                double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, angle);
                if (t1 >= 0 && t1 <= 1) {
                    ts0[tc] = t0;
                    ts1[tc] = t1;
                    is_samples[tc] = false;
                    ++tc;
                }
            }
        }
        else if (radius1 + radius0 <= center_center_distance) {
            WGVector2d vt = center_center_vt / center_center_distance;
            double angle0 = vt.UnitAngle();
            double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, angle0);
            if (t0 >= 0 && t0 <= 1) {
                double angle1 = angle0 + g_pi;
                if (angle1 >= g_pi * 2) {
                    angle1 -= g_pi * 2;
                }
                double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, angle1);
                if (t1 >= 0 && t1 <= 1) {
                    ts0[tc] = t0;
                    ts1[tc] = t1;
                    is_samples[tc] = false;
                    ++tc;
                }
            }
        }
        else {
            double a = (center_center_distance * center_center_distance + radius0 * radius0 - radius1 * radius1) / (2 * center_center_distance);
            WGVector2d vt = center_center_vt / center_center_distance;
            WGVector2d pt = center0 + vt * a;
            double d = sqrt(radius0 * radius0 - a * a);
            vt = WGVector2d(-vt.Y, vt.X);
            WGVector2d pt1 = pt - vt * d;
            WGVector2d pt2 = pt + vt * d;
            {
                double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, ((pt1 - center0) / radius0).UnitAngle());
                if (t0 >= 0 && t0 <= 1) {
                    double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, ((pt1 - center1) / radius1).UnitAngle());
                    if (t1 >= 0 && t1 <= 1) {
                        ts0[tc] = t0;
                        ts1[tc] = t1;
                        is_samples[tc] = false;
                        ++tc;
                    }
                }
            }
            {
                double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, ((pt2 - center0) / radius0).UnitAngle());
                if (t0 >= 0 && t0 <= 1) {
                    double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, ((pt2 - center1) / radius1).UnitAngle());
                    if (t1 >= 0 && t1 <= 1) {
                        ts0[tc] = t0;
                        ts1[tc] = t1;
                        is_samples[tc] = false;
                        ++tc;
                    }
                }
            }
        }    
    }
    int tc1 = tc;
    WGVector2d start_point0 = WGArc2d::CalculatePoint(center0, radius0, start_angle0);
    WGVector2d end_point0 = WGArc2d::CalculatePoint(center0, radius0, start_angle0 + delta_angle0);
    WGVector2d start_point1 = WGArc2d::CalculatePoint(center1, radius1, start_angle1);
    WGVector2d end_point1 = WGArc2d::CalculatePoint(center1, radius1, start_angle1 + delta_angle1);
    double epsilon2 = distance_epsilon * distance_epsilon;
    double end_ts0[4];
    double end_ts1[4];
    int end_tc = 0;
    WGVector2d vt = start_point1 - start_point0;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 0;
        ts1[tc] = 0;
        is_samples[tc] = true;
        ++tc;
        end_ts0[end_tc] = 0;
        end_ts1[end_tc] = 0;
        ++end_tc;
    }
    vt = end_point1 - start_point0;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 0;
        ts1[tc] = 1;
        is_samples[tc] = true;
        ++tc;
        end_ts0[end_tc] = 0;
        end_ts1[end_tc] = 1;
        ++end_tc;
    }
    vt = start_point1 - end_point0;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 1;
        ts1[tc] = 0;
        is_samples[tc] = true;
        ++tc;
        end_ts0[end_tc] = 1;
        end_ts1[end_tc] = 0;
        ++end_tc;
    }
    vt = end_point1 - end_point0;
    if (vt.SqrLength() <= epsilon2) {
        ts0[tc] = 1;
        ts1[tc] = 1;
        is_samples[tc] = true;
        ++tc;
        end_ts0[end_tc] = 1;
        end_ts1[end_tc] = 1;
        ++end_tc;
    }
    vt = start_point0 - center1;
    double d = vt.Length();
    if (d > 0 && abs(d - radius1) <= distance_epsilon) {
        double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, (vt / d).UnitAngle());
        if (t1 >= 0 && t1 <= 1) {
            ts0[tc] = 0;
            ts1[tc] = t1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = end_point0 - center1;
    d = vt.Length();
    if (d > g_double_epsilon && abs(d - radius1) <= distance_epsilon) {
        double t1 = WGArc2d::CalculateT(start_angle1, delta_angle1, (vt / d).UnitAngle());
        if (t1 >= 0 && t1 <= 1) {
            ts0[tc] = 1;
            ts1[tc] = t1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = start_point1 - center0;
    d = vt.Length();
    if (d > 0 && abs(d - radius0) <= distance_epsilon) {
        double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, (vt / d).UnitAngle());
        if (t0 >= 0 && t0 <= 1) {
            ts0[tc] = t0;
            ts1[tc] = 0;
            is_samples[tc] = true;
            ++tc;
        }
    }
    vt = end_point1 - center0;
    d = vt.Length();
    if (d > 0 && abs(d - radius0) <= distance_epsilon) {
        double t0 = WGArc2d::CalculateT(start_angle0, delta_angle0, (vt / d).UnitAngle());
        if (t0 >= 0 && t0 <= 1) {
            ts0[tc] = t0;
            ts1[tc] = 1;
            is_samples[tc] = true;
            ++tc;
        }
    }
    if (tc == 0) {
        return 0;
    }
    if (tc == tc1) {
        if (tc == 1) {
            intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
            return 1;
        }
        else {
            double extreme_ts0[4];
            double extreme_ts1[4];
            int extreme_count = BuildArcArcDisjointExtremes(center_center_vt / center_center_distance, center_center_distance,
                radius0, start_angle0, delta_angle0, radius1, start_angle1, delta_angle1, distance_epsilon, extreme_ts0, extreme_ts1);
            bool b = true;
            for (int i = 0; i < extreme_count; ++i) {
                if (((extreme_ts0[i] < ts0[0] && extreme_ts0[i] > ts0[1]) || (extreme_ts0[i] < ts0[1] && extreme_ts0[i] > ts0[0])) &&
                    ((extreme_ts1[i] < ts1[0] && extreme_ts1[i] > ts1[1]) || (extreme_ts1[i] < ts1[1] && extreme_ts1[i] > ts1[0]))) {
                    b = false;
                    break;
                }
            }
            if (b) {
                intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0], ts0[1], ts1[1], is_samples[1]);
                return 1;
            }
            else {
                intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
                intersections[1] = CurveCurveIntersection(ts0[1], ts1[1], is_samples[1]);
                return 2;
            }
        }
    }
    RemoveSameTs(ts0, ts1, is_samples, tc);
    if (tc == 1) {
        intersections[0] = CurveCurveIntersection(ts0[0], ts1[0], is_samples[0]);
        return 1;
    }
    WGVector2d points0[10];
    WGVector2d points1[10];
    for (int i = 0; i < tc; ++i) {
        points0[i] = WGArc2d::CalculatePoint(center0, radius0, start_angle0 + ts0[i] * delta_angle0);
        points1[i] = WGArc2d::CalculatePoint(center1, radius1, start_angle1 + ts1[i] * delta_angle1);
    }
    RepairCross(ts0, ts1, is_samples, points0, points1, tc, distance_epsilon);
    double extreme_ts0[4];
    double extreme_ts1[4];
    int extreme_count = BuildArcArcDisjointExtremes(center_center_vt / center_center_distance, center_center_distance,
        radius0, start_angle0, delta_angle0, radius1, start_angle1, delta_angle1, distance_epsilon, extreme_ts0, extreme_ts1);
    int indices[10];
    for (int i = 0; i < tc; ++i) {
        indices[i] = i;
    }
    int count = 0;
    SortTs(ts0, ts1, indices, tc);
    for (int j = 1; j < tc; ++j) {
        int i = j - 1;
        double t00 = ts0[indices[i]];
        double t01 = ts0[indices[j]];
        double t10 = ts1[indices[i]];
        double t11 = ts1[indices[j]];
        if (t00 != t01 || t10 != t11) {
            bool b = true;
            if (t00 == t01 && abs(delta_angle1 * (t11 - t10)) > g_pi) {
                b = false;
            }
            else if (t10 == t11 && abs(delta_angle1 * (t01 - t00)) > g_pi) {
                b = false;
            }
            else {
                for (int i = 0; i < extreme_count; ++i) {
                    if (((extreme_ts0[i] < t00 && extreme_ts0[i] > t01) || (extreme_ts0[i] < t01 && extreme_ts0[i] > t00)) &&
                        ((extreme_ts1[i] < t10 && extreme_ts1[i] > t11) || (extreme_ts1[i] < t11 && extreme_ts1[i] > t10))) {
                        b = false;
                        break;
                    }
                }
            }
            if (b) {
                intersections[count] = CurveCurveIntersection(t00, t10, is_samples[indices[i]], t01, t11, is_samples[indices[j]]);
                ++count;
            }
        }
    }
    for (int i = 0; i < tc; ++i) {
        bool b = true;
        for (int j = 0; j < count; ++j) {
            const CurveCurveIntersection& intersection = intersections[j];
            if (intersection.PointCount == 1) {
                if (intersection.Ts0[0] != ts0[i] || intersection.Ts1[0] != ts1[i]) {
                    continue;
                }
            }
            else {
                if (ts0[i] < intersection.Ts0[0] || ts0[i] > intersection.Ts0[1]) {
                    continue;
                }
                if (intersection.Ts1[0] <= intersection.Ts1[1]) {
                    if (ts1[i] < intersection.Ts1[0] || ts1[i] > intersection.Ts1[1]) {
                        continue;
                    }
                }
                else {
                    if (ts1[i] < intersection.Ts1[1] || ts1[i] > intersection.Ts1[0]) {
                        continue;
                    }
                }
            }
            b = false;
            break;
        }
        if (b) {
            intersections[count] = CurveCurveIntersection(ts0[i], ts1[i], is_samples[i]);
            ++count;
        }
    }
    MarkNotSample(ts0, ts1, is_samples, tc, intersections, count);
    count = MergeOverlapIntersections(intersections, count);
    RepairByEndTs(end_ts0, end_ts1, end_tc, intersections, count);
    return count;
}

int WGIntersectHelper2d::BuildArcArcDisjointExtremes(const WGVector2d& center_direction, double center_distance,
    double radius0, double start_angle0, double delta_angle0, double radius1, double start_angle1, double delta_angle1,
    double distance_epsilon, double extreme_ts0[4], double extreme_ts1[4]) {
    bool is_extreme_oks[4];
    int extreme_count = 0;
    double extreme_ds0[4];
    int extreme_count0 = 0;
    double angle1 = center_direction.UnitAngle();
    double angle2 = angle1 + g_pi;
    if (angle2 >= g_pi * 2) {
        angle2 -= g_pi * 2;
    }
    extreme_ts0[extreme_count0] = WGArc2d::CalculateT(start_angle0, delta_angle0, angle1);
    if (extreme_ts0[extreme_count0] >= 0 && extreme_ts0[extreme_count0] <= 1) {
        extreme_ds0[extreme_count0] = radius0;
        ++extreme_count0;
    }
    extreme_ts0[extreme_count0] = WGArc2d::CalculateT(start_angle0, delta_angle0, angle2);
    if (extreme_ts0[extreme_count0] >= 0 && extreme_ts0[extreme_count0] <= 1) {
        extreme_ds0[extreme_count0] = -radius0;
        ++extreme_count0;
    }
    if (extreme_count0 > 0) {
        double extreme_ds1[4];
        int extreme_count1 = 0;
        extreme_ts1[extreme_count1] = WGArc2d::CalculateT(start_angle1, delta_angle1, angle1);
        if (extreme_ts1[extreme_count1] >= 0 && extreme_ts1[extreme_count1] <= 1) {
            extreme_ds1[extreme_count1] = radius1;
            ++extreme_count1;
        }
        extreme_ts1[extreme_count1] = WGArc2d::CalculateT(start_angle1, delta_angle1, angle2);
        if (extreme_ts1[extreme_count1] >= 0 && extreme_ts1[extreme_count1] <= 1) {
            extreme_ds1[extreme_count1] = -radius1;
            ++extreme_count1;
        }
        if (extreme_count1 > 0) {
            if (extreme_count0 == 1) {
                if (extreme_count1 == 1) {
                    is_extreme_oks[0] = abs(extreme_ds1[0] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    extreme_count = 1;
                }
                else {
                    extreme_ts0[1] = extreme_ts0[0];
                    is_extreme_oks[0] = abs(extreme_ds1[0] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    is_extreme_oks[1] = abs(extreme_ds1[1] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    extreme_count = 2;
                }
            }
            else {
                if (extreme_count1 == 1) {
                    extreme_ts1[1] = extreme_ts1[0];
                    is_extreme_oks[0] = abs(extreme_ds1[0] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    is_extreme_oks[1] = abs(extreme_ds1[0] + center_distance - extreme_ds0[1]) <= distance_epsilon;
                    extreme_count = 2;
                }
                else {
                    extreme_ts0[2] = extreme_ts0[1];
                    extreme_ts0[3] = extreme_ts0[1];
                    extreme_ts0[1] = extreme_ts0[0];
                    extreme_ts1[2] = extreme_ts1[0];
                    extreme_ts1[3] = extreme_ts1[1];
                    is_extreme_oks[0] = abs(extreme_ds1[0] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    is_extreme_oks[1] = abs(extreme_ds1[1] + center_distance - extreme_ds0[0]) <= distance_epsilon;
                    is_extreme_oks[2] = abs(extreme_ds1[0] + center_distance - extreme_ds0[1]) <= distance_epsilon;
                    is_extreme_oks[3] = abs(extreme_ds1[1] + center_distance - extreme_ds0[1]) <= distance_epsilon;
                    extreme_count = 4;
                }
            }
            int i = 0;
            int n = 0;
            while (i < extreme_count) {
                if (!is_extreme_oks[i]) {
                    if (i != n) {
                        extreme_ts0[n] = extreme_ts0[i];
                        extreme_ts1[n] = extreme_ts1[i];
                    }
                    ++n;
                }
                ++i;
            }
            extreme_count = n;
        }
    }
    return extreme_count;
}

bool WGIntersectHelper2d::Sample(PointCurveEquationSystem* equations, const WGVector2d& point, const WGVector2d& sample_direction,
    const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance) {
    {
        sample_distance = 0;
        WGVector2d p0 = point;
        sample_t = sample_domain.Middle();
        WGVector2d p1;
        equations->CalculateD0(sample_t, p1);
        WGVector2d v = p1 - p0;
        double d = v.Length();
        if (d <= epsilon) {
            return true;
        }
        WGVector2d v1;
        equations->CalculateD1(sample_t, v1);
        for (int itr_count = 0; itr_count <= 30; ++itr_count) {
            double a = v1.Cross(sample_direction);
            if (a == 0) {
                return false;
            }
            double dt0 = v1.Cross(v) / a;
            double dt1 = sample_direction.Cross(v) / a;
            double t1 = sample_t + dt1;
            if (t1 < sample_domain.Min) {
                t1 = sample_domain.Min;
                double dt = t1 - sample_t;
                dt0 *= dt / dt1;
                dt1 = dt;
            }
            else if (t1 > sample_domain.Max) {
                t1 = sample_domain.Max;
                double dt = t1 - sample_t;
                dt0 *= dt / dt1;
                dt1 = dt;
            }
            if (abs(dt0) <= 1E-16 && abs(dt1) <= 1E-16) {
                break;
            }
            double nt0 = sample_distance + dt0;
            double nt1 = sample_t + dt1;
            p0 = point + sample_direction * nt0;
            equations->CalculateD0(nt1, p1);
            v = p1 - p0;
            double nd = v.Length();
            bool b = false;
            while (nd > d) {
                dt0 *= 0.5;
                dt1 *= 0.5;
                if (abs(dt0) <= 1E-16 && abs(dt1) <= 1E-16) {
                    b = true;
                    break;
                }
                nt0 = sample_distance + dt0;
                nt1 = sample_t + dt1;
                p0 = point + sample_direction * nt0;
                equations->CalculateD0(nt1, p1);
                v = p1 - p0;
                nd = v.Length();
            }
            if (b) {
                break;
            }
            sample_distance = nt0;
            sample_t = nt1;
            d = nd;
            if (d <= epsilon) {
                return true;
            }
            equations->CalculateD1(sample_t, v1);
        }
    }
    {
        double t0 = sample_domain.Min;
        double t1 = sample_domain.Max;
        WGVector2d p0, p1;
        equations->CalculateD0(t0, p0);
        equations->CalculateD0(t1, p1);
        WGVector2d v0 = p0 - point;
        double d0 = v0.Cross(sample_direction);
        WGVector2d v1 = p1 - point;
        double d1 = v1.Cross(sample_direction);
        double a0 = abs(d0);
        double a1 = abs(d1);
        if (a0 < a1) {
            if (a0 <= epsilon) {
                sample_t = t0;
                sample_distance = v0.Dot(sample_direction);
                return true;
            }
        }
        else {
            if (a1 <= epsilon) {
                sample_t = t1;
                sample_distance = v1.Dot(sample_direction);
                return true;
            }
        }
        if (d0 * d1 > 0) {
            return false;
        }
        while (true) {
            double t = (t0 + t1) * 0.5;
            WGVector2d p;
            equations->CalculateD0(t, p);
            WGVector2d v = p - point;
            double d = v.Cross(sample_direction);
            if (abs(d) <= epsilon) {
                sample_t = t;
                sample_distance = v.Dot(sample_direction);
                return true;
            }
            if (d0 * d < 0) {
                t1 = t;
            }
            else {
                t0 = t;
            }
            if (t1 - t0 <= g_double_epsilon) {
                sample_t = t;
                sample_distance = v.Dot(sample_direction);
                return true;
            }
        }
    }
}

void WGIntersectHelper2d::BuildFuzzyIntersection(PointCurveEquationSystem* equations, const WGVector2d& point, 
    const WSIntervalVector* variable, double distance_epsilon, std::vector<PointCurveIntersection>& intersections) {
    WSInterval dx, dy;
    equations->CalculateD1(variable, dx, dy);
    WGVector2d vt(dx.Middle(), dy.Middle());
    vt.Normalize(0);
    WGVector2d sample_direction(-vt.Y, vt.X);
    double sample_epsilon = distance_epsilon * 0.001;
    double sample_t;
    double sample_distance;
    if (Sample(equations, point, sample_direction, variable->Get(0), sample_epsilon, sample_t, sample_distance)) {
        double d2 = sample_distance * sample_distance;
        WGVector2d point2;
        const WSInterval& t = variable->Get(0);
        equations->CalculateD0(t.Min, point2);
        double d = (point2 - point).SqrLength();
        if (d < d2) {
            d2 = d;
            sample_t = t.Min;
        }
        equations->CalculateD0(t.Max, point2);
        d = (point2 - point).SqrLength();
        if (d < d2) {
            d2 = d;
            sample_t = t.Max;
        }
        if (d2 <= distance_epsilon * distance_epsilon) {
            intersections.push_back(PointCurveIntersection(sample_t));
        }
    }
}

WSInterval WGIntersectHelper2d::CalculateBezierD1(int degree, const double* coefs, const WSInterval& t) {
    assert(degree < 16);
    if (degree < 1) {
        return 0;
    }
    double d = coefs[1] - coefs[0];
    if (degree == 1) {
        return d;
    }
    double coefs1[16];
    coefs1[0] = d;
    for (int i = 2; i <= degree; ++i) {
        coefs1[i - 1] = coefs[i] - coefs[i - 1];
    }
    int degree1 = degree - 1;
    WSBernsteinCalculator::SubSection(degree1, coefs1, t);
    WSInterval v = coefs1[0];
    for (int i = 1; i <= degree1; ++i) {
        v.Merge(coefs1[i]);
    }
    return v * (double)degree;
}

WSInterval WGIntersectHelper2d::CalculateBezierD2(int degree, const double* coefs, const WSInterval& t) {
    assert(degree < 16);
    if (degree < 2) {
        return 0;
    }
    double d = coefs[2] - 2 * coefs[1] + coefs[0];
    if (degree == 2) {
        return d * 2;
    }
    double coefs2[16];
    coefs2[0] = d;
    for (int i = 3; i <= degree; ++i) {
        coefs2[i - 2] = coefs[i] - 2 * coefs[i - 1] + coefs[i - 2];
    }
    int degree2 = degree - 2;
    WSBernsteinCalculator::SubSection(degree2, coefs2, t);
    WSInterval v = coefs2[0];
    for (int i = 1; i <= degree2; ++i) {
        v.Merge(coefs2[i]);
    }
    return v * ((double)degree * (degree - 1));
}

WSInterval WGIntersectHelper2d::CalculateBezierD3(int degree, const double* coefs, const WSInterval& t) {
    assert(degree < 16);
    if (degree < 3) {
        return 0;
    }
    double d = coefs[3] - 3 * coefs[2] + 3 * coefs[1] - coefs[0];
    if (degree == 3) {
        return d * 6;
    }
    double coefs3[16];
    coefs3[0] = d;
    for (int i = 4; i <= degree; ++i) {
        coefs3[i - 3] = coefs[i] - 3 * coefs[i - 1] + 3 * coefs[i - 2] - coefs[i - 3];
    }
    int degree3 = degree - 3;
    WSBernsteinCalculator::SubSection(degree3, coefs3, t);
    WSInterval v = coefs3[0];
    for (int i = 1; i <= degree3; ++i) {
        v.Merge(coefs3[i]);
    }
    return v * ((double)degree * (degree - 1) * (degree - 2));
}

WSInterval WGIntersectHelper2d::CalculateBezierD4(int degree, const double* coefs, const WSInterval& t) {
    assert(degree < 16);
    if (degree < 4) {
        return 0;
    }
    double d = coefs[4] - 4 * coefs[3] + 6 * coefs[2] - 4 * coefs[1] + coefs[0];
    if (degree == 4) {
        return d * 24;
    }
    double coefs4[16];
    coefs4[0] = d;
    for (int i = 5; i <= degree; ++i) {
        coefs4[i - 4] = coefs[i] - 4 * coefs[i - 1] + 6 * coefs[i - 2] - 4 * coefs[i - 3] + coefs[i - 4];
    }
    int degree4 = degree - 4;
    WSBernsteinCalculator::SubSection(degree4, coefs4, t);
    WSInterval v = coefs4[0];
    for (int i = 1; i <= degree4; ++i) {
        v.Merge(coefs4[i]);
    }
    return v * ((double)degree * (degree - 1) * (degree - 2) * (degree - 3));
}

void WGIntersectHelper2d::PreparePointCurveIntersect(PointCurveEquationSystem* equations, PointCurveSolver* solver,
    const WSIntervalVector* variable_domain, double distance_epsilon) {
    solver->Execute(equations, variable_domain, 1024, 1024);
}

void WGIntersectHelper2d::BuildPointCurveIntersections(PointCurveEquationSystem* equations, PointCurveSolver* solver,
    const WGVector2d& point, double distance_epsilon, std::vector<PointCurveIntersection>& intersections) {
    for (int i = 0; i < solver->GetClearRootCount(); ++i) {
        WGVector2d point2;
        const WSInterval& t = solver->GetClearRoot(i)->Get(0);
        equations->CalculateD0(t.Min, point2);
        double d1 = (point2 - point).SqrLength();
        equations->CalculateD0(t.Max, point2);
        double d2 = (point2 - point).SqrLength();
        if (d1 < d2) {
            intersections.push_back(PointCurveIntersection(t.Min));
        }
        else {
            intersections.push_back(PointCurveIntersection(t.Max));
        }
    }
    for (int i = 0; i < solver->GetFuzzyDomainCount(); ++i) {
        if ((solver->GetFuzzyTag(i) & 1) == 1) {
            BuildFuzzyIntersection(equations, point, solver->GetFuzzyDomain(i), distance_epsilon, intersections);
        }
    }
}

void WGIntersectHelper2d::PointBezierCurveIntersect(const WGVector2d& point, int degree, const WGVector2d* control_points,
    double distance_epsilon, IntersectCache* intersect_cache, bool& is_singularity, std::vector<PointCurveIntersection>& intersections) {
    PointBezierCurveEquationSystem* equations = intersect_cache->AllocPointBezierCurveEquations();
    PointCurveSolver* solver = intersect_cache->AllocPointCurveSolver();
    equations->Update(point, degree, control_points, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    PreparePointCurveIntersect(equations, solver, &variable, distance_epsilon);
    BuildPointCurveIntersections(equations, solver, point, distance_epsilon, intersections);
    is_singularity = solver->GetIsSingularity();
    intersect_cache->FreePointCurveSolver(solver);
    intersect_cache->FreePointBezierCurveEquations(equations);
}

void WGIntersectHelper2d::PointRationalBezierCurveIntersect(const WGVector2d& point, int degree, const WGVector2d* control_points, const double* weights,
    double distance_epsilon, IntersectCache* intersect_cache, bool& is_singularity, std::vector<PointCurveIntersection>& intersections) {
    PointRationalBezierCurveEquationSystem* equations = intersect_cache->AllocPointRationalBezierCurveEquations();
    PointCurveSolver* solver = intersect_cache->AllocPointCurveSolver();
    equations->Update(point, degree, control_points, weights, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    PreparePointCurveIntersect(equations, solver, &variable, distance_epsilon);
    BuildPointCurveIntersections(equations, solver, point, distance_epsilon, intersections);
    is_singularity = solver->GetIsSingularity();
    intersect_cache->FreePointCurveSolver(solver);
    intersect_cache->FreePointRationalBezierCurveEquations(equations);
}

int WGIntersectHelper2d::PointLineIntersect(const WGVector2d& point, const WGVector2d& start_point, const WGVector2d& end_point,
    double distance_epsilon, bool& is_singularity, PointCurveIntersection intersections[1]) {
    WGVector2d vt = end_point - start_point;
    if (is_zero(vt.X, g_double_epsilon) && is_zero(vt.Y, g_double_epsilon)) {
        is_singularity = true;
        return 0;
    }
    is_singularity = false;
    double d = vt.Normalize(0);
    WGVector2d vt2 = point - start_point;
    double t = vt2.Dot(vt) / d;
    if (t >= 0 && t <= 1) {
        if (is_zero(vt2.Cross(vt), distance_epsilon)) {
            intersections[0] = PointCurveIntersection(t);
            return 1;
        }
    }
    return 0;
}

int WGIntersectHelper2d::PointArcIntersect(const WGVector2d& point, const WGVector2d& center, double radius, double start_angle, double delta_angle,
    double distance_epsilon, bool& is_singularity, PointCurveIntersection intersections[2]) {
    is_singularity = radius <= g_double_epsilon || is_zero(delta_angle, g_double_epsilon);
    if (is_singularity) {
        return 0;
    }
    WGVector2d vt = point - center;
    if (is_zero(vt.X, g_double_epsilon) && is_zero(vt.Y, g_double_epsilon)) {
        is_singularity = true;
        return 0;
    }
    double d = vt.Normalize(0);
    if (is_zero(d - radius, distance_epsilon)) {
        double t = WGArc2d::CalculateT(start_angle, delta_angle, vt.UnitAngle());
        if (t == 0) {
            intersections[0] = PointCurveIntersection(0);
            if (abs(delta_angle) == g_pi * 2) {
                intersections[1] = PointCurveIntersection(1);
                return 2;
            }
            return 1;
        }
        if (t > 0 && t <= 1) {
            intersections[0] = PointCurveIntersection(t);
            return 1;
        }
    }
    return 0;
}

void WGIntersectHelper2d::PrepareCurveCurveIntersect(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
    const WSIntervalVector* variable_domain, double distance_epsilon) {
    solver->PrepareExecute(equations, variable_domain, 1024, 1024);
    solver->PrimaryExecute();
}

void WGIntersectHelper2d::BuildCurveCurveIntersections(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
    WSIntervalVector* runtime_variable, double distance_epsilon, IntersectCache* intersect_cache, 
    std::vector<CurveCurveIntersection>& intersections) {
    for (int i = 0; i < solver->GetClearRootCount(); ++i) {
        intersections.push_back(CurveCurveIntersection(solver->GetClearRoot(i)->Get(0).Middle(), solver->GetClearRoot(i)->Get(1).Middle(), false));
    }
    CurveCurveFuzzyPair* first_fuzzy_pair = nullptr;
    for (int i = 0; i < solver->GetFuzzyDomainCount(); ++i) {
        if ((solver->GetFuzzyTag(i) & 1) != 0) {
            BuildMonotonicIntersection(equations, distance_epsilon, solver->GetFuzzyDomain(i)->Get(0),
                solver->GetFuzzyDomain(i)->Get(1), intersections);
        }
        else {
            CurveCurveFuzzyPair* fuzzy_pair = BuildFuzzyPair(equations, runtime_variable, distance_epsilon,
                solver->GetFuzzyDomain(i)->Get(0), solver->GetFuzzyDomain(i)->Get(1), intersect_cache);
            if (fuzzy_pair) {
                if (abs(fuzzy_pair->SamplePoints[0].Distance) > distance_epsilon || abs(fuzzy_pair->SamplePoints[1].Distance) > distance_epsilon) {
                    WSIterateResult iterate_result = solver->FinallyIterate(i, distance_epsilon);
                    if (iterate_result == WSIterateResult::NoRoot) {
                        intersect_cache->FreePair(fuzzy_pair);
                    }
                    else if (iterate_result == WSIterateResult::ClearRoot) {
                        intersections.push_back(CurveCurveIntersection(solver->GetClearRoot(i)->Get(0).Middle(),
                            solver->GetClearRoot(i)->Get(1).Middle(), false));
                        intersect_cache->FreePair(fuzzy_pair);
                    }
                    else {
                        WSInterval v0 = solver->GetFuzzyDomain(i)->Get(0);
                        WSInterval v1 = solver->GetFuzzyDomain(i)->Get(1);
                        if (fuzzy_pair->SamplePoints[0].T0 >= v0.Min && fuzzy_pair->SamplePoints[0].T0 <= v0.Max &&
                            fuzzy_pair->SamplePoints[0].T1 >= v0.Min && fuzzy_pair->SamplePoints[0].T1 <= v0.Max &&
                            fuzzy_pair->SamplePoints[1].T0 >= v0.Min && fuzzy_pair->SamplePoints[1].T0 <= v0.Max &&
                            fuzzy_pair->SamplePoints[1].T1 >= v0.Min && fuzzy_pair->SamplePoints[1].T1 <= v0.Max) {
                            fuzzy_pair->Next = first_fuzzy_pair;
                            first_fuzzy_pair = fuzzy_pair;
                        }
                        else {
                            intersect_cache->FreePair(fuzzy_pair);
                            fuzzy_pair = BuildFuzzyPair(equations, runtime_variable, distance_epsilon, v0, v1, intersect_cache);
                            if (fuzzy_pair) {
                                fuzzy_pair->Next = first_fuzzy_pair;
                                first_fuzzy_pair = fuzzy_pair;
                            }
                        }
                    }
                }
                else {
                    fuzzy_pair->Next = first_fuzzy_pair;
                    first_fuzzy_pair = fuzzy_pair;
                }
            }
        }
    }
    BuildFuzzyIntersection(equations, solver, runtime_variable, distance_epsilon, 
        first_fuzzy_pair, intersect_cache, intersections);
}

void WGIntersectHelper2d::MergeNotSampleIntersections(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable, 
    std::vector<CurveCurveIntersection>& intersections) {
    int n = 0;
    for (int i = 0; i < (int)intersections.size(); ++i) {
        const CurveCurveIntersection& intersection1 = intersections.at(i);
        bool b = true;
        if (intersection1.PointCount == 1 && !intersection1.IsSamples[0]) {
            for (int j = i + 1; j < (int)intersections.size(); ++j) {
                CurveCurveIntersection& intersection2 = intersections.at(j);
                if (intersection2.PointCount == 1 && !intersection2.IsSamples[0]) {
                    if (intersection1.Ts0[0] < intersection2.Ts0[0]) {
                        runtime_variable->Set(0, WSInterval(intersection1.Ts0[0], intersection2.Ts0[0]));
                    }
                    else {
                        runtime_variable->Set(0, WSInterval(intersection2.Ts0[0], intersection1.Ts0[0]));
                    }
                    if (intersection1.Ts1[0] < intersection2.Ts1[0]) {
                        runtime_variable->Set(1, WSInterval(intersection1.Ts1[0], intersection2.Ts1[0]));
                    }
                    else {
                        runtime_variable->Set(1, WSInterval(intersection2.Ts1[0], intersection1.Ts1[0]));
                    }
                    if (CheckMonotonic(equations, runtime_variable)) {
                        WGVector2d point0;
                        WGVector2d point1;
                        equations->CalculateD0(0, intersection1.Ts0[0], point0);
                        equations->CalculateD0(1, intersection1.Ts1[0], point1);
                        double d1 = (point1 - point0).SqrLength();
                        equations->CalculateD0(0, intersection2.Ts0[0], point0);
                        equations->CalculateD0(1, intersection2.Ts1[0], point1);
                        double d2 = (point1 - point0).SqrLength();
                        if (d1 < d2) {
                            intersection2.Ts0[0] = intersection1.Ts0[0];
                            intersection2.Ts1[0] = intersection1.Ts1[0];
                        }
                        b = false;
                        break;
                    }
                }
            }
        }
        if (b) {
            if (i != n) {
                intersections.at(n) = intersection1;
            }
            ++n;
        }
    }
    intersections.resize(n);
}

void WGIntersectHelper2d::AdjustIntersectionsWithEndpoints(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable,
    const std::vector<PointCurveIntersection>& intersections00, const std::vector<PointCurveIntersection>& intersections01,
    const std::vector<PointCurveIntersection>& intersections10, const std::vector<PointCurveIntersection>& intersections11,
    double distance_epsilon, std::vector<CurveCurveIntersection>& intersections) {
    int n = (int)(intersections00.size() + intersections01.size() + intersections10.size() + intersections11.size());
    if (n == 0) {
        return;
    }
    int c = n + (int)intersections.size() * 2;
    WSCache cache;
    cache.Resize(c * (sizeof(double) + sizeof(double) + sizeof(bool) + sizeof(WGVector2d) + sizeof(WGVector2d) + sizeof(int)));
    int cache_offset = 0;
    double* ts0 = (double*)cache.Data(cache_offset);
    cache_offset += sizeof(double) * c;
    double* ts1 = (double*)cache.Data(cache_offset);
    cache_offset += sizeof(double) * c;
    bool* is_samples = (bool*)cache.Data(cache_offset);
    cache_offset += sizeof(bool) * c;
    WGVector2d* points0 = (WGVector2d*)cache.Data(cache_offset);
    cache_offset += sizeof(WGVector2d) * c;
    WGVector2d* points1 = (WGVector2d*)cache.Data(cache_offset);
    cache_offset += sizeof(WGVector2d) * c;
    int* indices = (int*)cache.Data(cache_offset);
    int tc = 0;
    for (int i = 0; i < (int)intersections00.size(); ++i) {
        ts0[tc] = 0;
        ts1[tc] = intersections00.at(i).T;
        is_samples[tc] = true;
        ++tc;
    }
    for (int i = 0; i < (int)intersections01.size(); ++i) {
        ts0[tc] = 1;
        ts1[tc] = intersections01.at(i).T;
        is_samples[tc] = true;
        ++tc;
    }
    for (int i = 0; i < (int)intersections10.size(); ++i) {
        ts0[tc] = intersections10.at(i).T;
        ts1[tc] = 0;
        is_samples[tc] = true;
        ++tc;
    }
    for (int i = 0; i < (int)intersections11.size(); ++i) {
        ts0[tc] = intersections11.at(i).T;
        ts1[tc] = 1;
        is_samples[tc] = true;
        ++tc;
    }
    for (int i = 0; i < (int)intersections.size(); ++i) {
        const CurveCurveIntersection& intersection = intersections.at(i);
        ts0[tc] = intersection.Ts0[0];
        ts1[tc] = intersection.Ts1[0];
        is_samples[tc] = intersection.IsSamples[0];
        ++tc;
        if (intersection.PointCount == 2) {
            ts0[tc] = intersection.Ts0[1];
            ts1[tc] = intersection.Ts1[1];
            is_samples[tc] = intersection.IsSamples[1];
            ++tc;
        }
    }
    RemoveSameTs(ts0, ts1, is_samples, tc);
    for (int i = 0; i < tc; ++i) {
        equations->CalculateD0(0, ts0[i], points0[i]);
        equations->CalculateD0(1, ts1[i], points1[i]);
    }
    RepairCross(ts0, ts1, is_samples, points0, points1, tc, distance_epsilon);
    int i = 0;
    while (i < tc) {
        bool b = true;
        if (!is_samples[i]) {
            for (int j = 0; j < tc; ++j) {
                if (i != j) {
                    if (ts0[j] == 0 || ts0[j] == 1 || ts1[j] == 0 || ts1[j] == 1) {
                        if (ts0[i] <= ts0[j]) {
                            runtime_variable->Set(0, WSInterval(ts0[i], ts0[j]));
                        }
                        else {
                            runtime_variable->Set(0, WSInterval(ts0[j], ts0[i]));
                        }
                        if (ts1[i] <= ts1[j]) {
                            runtime_variable->Set(1, WSInterval(ts1[i], ts1[j]));
                        }
                        else {
                            runtime_variable->Set(1, WSInterval(ts1[j], ts1[i]));
                        }
                        if (CheckMonotonic(equations, runtime_variable)) {
                            double d1 = (points1[i] - points0[i]).SqrLength();
                            double d2 = (points1[j] - points0[j]).SqrLength();
                            if (d2 <= d1) {
                                is_samples[j] = false;
                                b = false;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (b) {
            ++i;
        }
        else {
            --tc;
            ts0[i] = ts0[tc];
            ts1[i] = ts1[tc];
            is_samples[i] = is_samples[tc];
            points0[i] = points0[tc];
            points1[i] = points1[tc];
        }
    }
    std::vector<CurveCurveIntersection> old_intersections = std::move(intersections);
    intersections.reserve(tc);
    for (int i = 0; i < tc; ++i) {
        indices[i] = i;
    }
    SortTs(ts0, ts1, indices, tc);
    for (int j = 1; j < tc; ++j) {
        int i = j - 1;
        double t00 = ts0[indices[i]];
        double t10 = ts1[indices[i]];
        double t01 = ts0[indices[j]];
        double t11 = ts1[indices[j]];
        if (t00 != t01 || t10 != t11) {
            bool is_sample0 = is_samples[indices[i]];
            bool is_sample1 = is_samples[indices[j]];
            bool b = true;
            for (int k = 0; k < (int)old_intersections.size(); ++k) {
                const CurveCurveIntersection& old_intersection = old_intersections.at(k);
                if (old_intersection.PointCount == 1) {
                    continue;
                }
                if (t00 >= old_intersection.Ts0[0] && t01 <= old_intersection.Ts0[1]) {
                    if (((t01 - t00) * (t11 - t10) >= 0) ==
                        ((old_intersection.Ts0[1] - old_intersection.Ts0[0]) * (old_intersection.Ts1[1] - old_intersection.Ts1[0]) >= 0)) {
                        if (old_intersection.Ts1[0] < old_intersection.Ts1[1]) {
                            if (t10 >= old_intersection.Ts1[0] && t10 <= old_intersection.Ts1[1] &&
                                t11 >= old_intersection.Ts1[0] && t11 <= old_intersection.Ts1[1]) {
                                intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                                b = false;
                                break;
                            }
                        }
                        else {
                            if (t10 >= old_intersection.Ts1[1] && t10 <= old_intersection.Ts1[0] &&
                                t11 >= old_intersection.Ts1[1] && t11 <= old_intersection.Ts1[0]) {
                                intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                                b = false;
                                break;
                            }
                        }
                    }
                }
            }
            if (b) {
                if (t10 <= t11) {
                    runtime_variable->Set(0, WSInterval(t00, t01));
                    runtime_variable->Set(1, WSInterval(t10, t11));
                }
                else {
                    runtime_variable->Set(0, WSInterval(t00, t01));
                    runtime_variable->Set(1, WSInterval(t11, t10));
                }
                if (CheckMonotonicOrShortLinear(equations, runtime_variable, g_curve_primary_linear_limit_angle, distance_epsilon)) {
                    intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                }
            }
        }
    }
    SortTs(ts1, ts0, indices, tc);
    for (int j = 1; j < tc; ++j) {
        int i = j - 1;
        double t00 = ts0[indices[i]];
        double t10 = ts1[indices[i]];
        double t01 = ts0[indices[j]];
        double t11 = ts1[indices[j]];
        if (t00 != t01 || t10 != t11) {
            bool is_sample0 = is_samples[indices[i]];
            bool is_sample1 = is_samples[indices[j]];
            if (t00 > t01) {
                double t = t00;
                t00 = t01;
                t01 = t;
                t = t10;
                t10 = t11;
                t11 = t;
                bool b = is_sample0;
                is_sample0 = is_sample1;
                is_sample1 = b;
            }
            bool b = true;
            for (int k = 0; k < (int)intersections.size(); ++k) {
                const CurveCurveIntersection& intersection = intersections.at(k);
                if (intersection.PointCount == 1) {
                    continue;
                }
                if (t00 == intersection.Ts0[0] && t01 == intersection.Ts0[1]) {
                    if ((t10 == intersection.Ts1[0] && t11 == intersection.Ts1[1]) ||
                        (t10 == intersection.Ts1[1] && t11 == intersection.Ts1[0])) {
                        b = false;
                        break;
                    }
                }
            }
            if (b) {
                for (int k = 0; k < (int)old_intersections.size(); ++k) {
                    const CurveCurveIntersection& old_intersection = old_intersections.at(k);
                    if (old_intersection.PointCount == 1) {
                        continue;
                    }
                    if (t00 >= old_intersection.Ts0[0] && t01 <= old_intersection.Ts0[1]) {
                        if (((t01 - t00) * (t11 - t10) >= 0) ==
                            ((old_intersection.Ts0[1] - old_intersection.Ts0[0]) * (old_intersection.Ts1[1] - old_intersection.Ts1[0]) >= 0)) {
                            if (old_intersection.Ts1[0] < old_intersection.Ts1[1]) {
                                if (t10 >= old_intersection.Ts1[0] && t10 <= old_intersection.Ts1[1] &&
                                    t11 >= old_intersection.Ts1[0] && t11 <= old_intersection.Ts1[1]) {
                                    intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                                    b = false;
                                    break;
                                }
                            }
                            else {
                                if (t10 >= old_intersection.Ts1[1] && t10 <= old_intersection.Ts1[0] &&
                                    t11 >= old_intersection.Ts1[1] && t11 <= old_intersection.Ts1[0]) {
                                    intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                                    b = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if (b) {
                if (t10 <= t11) {
                    runtime_variable->Set(0, WSInterval(t00, t01));
                    runtime_variable->Set(1, WSInterval(t10, t11));
                }
                else {
                    runtime_variable->Set(0, WSInterval(t00, t01));
                    runtime_variable->Set(1, WSInterval(t11, t10));
                }
                if (CheckMonotonicOrShortLinear(equations, runtime_variable, g_curve_primary_linear_limit_angle, distance_epsilon)) {
                    intersections.push_back(CurveCurveIntersection(t00, t10, is_sample0, t01, t11, is_sample1));
                }
            }
        }
    }
    for (int i = 0; i < tc; ++i) {
        bool b = true;
        for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
            const CurveCurveIntersection& intersection = *intersection_itr;
            if (intersection.PointCount == 1) {
                if (intersection.Ts0[0] != ts0[i] || intersection.Ts1[0] != ts1[i]) {
                    continue;
                }
            }
            else {
                if (ts0[i] < intersection.Ts0[0] || ts0[i] > intersection.Ts0[1]) {
                    continue;
                }
                if (intersection.Ts1[0] <= intersection.Ts1[1]) {
                    if (ts1[i] < intersection.Ts1[0] || ts1[i] > intersection.Ts1[1]) {
                        continue;
                    }
                }
                else {
                    if (ts1[i] < intersection.Ts1[1] || ts1[i] > intersection.Ts1[0]) {
                        continue;
                    }
                }
            }
            b = false;
            break;
        }
        if (b) {
            intersections.push_back(CurveCurveIntersection(ts0[i], ts1[i], is_samples[i]));
        }
    }
    MarkNotSample(ts0, ts1, is_samples, tc, intersections.data(), (int)intersections.size());
    n = MergeOverlapIntersections(intersections.data(), (int)intersections.size());
    intersections.resize(n);
}

void WGIntersectHelper2d::LineBezierCurveIntersect(const WGVector2d& start_point, const WGVector2d& end_point,
    int degree, const WGVector2d* control_points, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    LineBezierCurveEquationSystem* equations = intersect_cache->AllocLineBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(start_point, end_point, degree, control_points, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointBezierCurveIntersect(start_point, degree, control_points, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointBezierCurveIntersect(end_point, degree, control_points, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointCurveIntersection temp_intersections[2];
        int n = PointLineIntersect(control_points[0], start_point, end_point, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections10.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections10.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        n = PointLineIntersect(control_points[degree], start_point, end_point, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections11.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections11.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((start_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((start_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((end_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((end_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeLineBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::LineRationalBezierCurveIntersect(const WGVector2d& start_point, const WGVector2d& end_point,
    int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    LineRationalBezierCurveEquationSystem* equations = intersect_cache->AllocLineRationalBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(start_point, end_point, degree, control_points, weights, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointRationalBezierCurveIntersect(start_point, degree, control_points, weights, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointRationalBezierCurveIntersect(end_point, degree, control_points, weights, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointCurveIntersection temp_intersections[2];
        int n = PointLineIntersect(control_points[0], start_point, end_point, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections10.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections10.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        n = PointLineIntersect(control_points[degree], start_point, end_point, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections11.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections11.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((start_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((start_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((end_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((end_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeLineRationalBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::ArcBezierCurveIntersect(const WGVector2d& center, double radius, double start_angle, double delta_angle,
    int degree, const WGVector2d* control_points, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    ArcBezierCurveEquationSystem* equations = intersect_cache->AllocArcBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(center, radius, start_angle, delta_angle, degree, control_points, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        WGVector2d start_point = WGArc2d::CalculatePoint(center, radius, start_angle);
        WGVector2d end_point = WGArc2d::CalculatePoint(center, radius, start_angle + delta_angle);
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointBezierCurveIntersect(start_point, degree, control_points, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointBezierCurveIntersect(end_point, degree, control_points, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointCurveIntersection temp_intersections[2];
        int n = PointArcIntersect(control_points[0], center, radius, start_angle, delta_angle, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections10.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections10.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        n = PointArcIntersect(control_points[degree], center, radius, start_angle, delta_angle, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections11.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections11.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((start_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((start_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((end_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((end_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeArcBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::ArcRationalBezierCurveIntersect(const WGVector2d& center, double radius, double start_angle, double delta_angle,
    int degree, const WGVector2d* control_points, const double* weights, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    ArcRationalBezierCurveEquationSystem* equations = intersect_cache->AllocArcRationalBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(center, radius, start_angle, delta_angle, degree, control_points, weights, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        WGVector2d start_point = WGArc2d::CalculatePoint(center, radius, start_angle);
        WGVector2d end_point = WGArc2d::CalculatePoint(center, radius, start_angle + delta_angle);
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointRationalBezierCurveIntersect(start_point, degree, control_points, weights, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointRationalBezierCurveIntersect(end_point, degree, control_points, weights, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointCurveIntersection temp_intersections[2];
        int n = PointArcIntersect(control_points[0], center, radius, start_angle, delta_angle, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections10.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections10.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        n = PointArcIntersect(control_points[degree], center, radius, start_angle, delta_angle, distance_epsilon, is_singularity, temp_intersections);
        if (n > 0) {
            intersections11.reserve(n);
            for (int i = 0; i < n; ++i) {
                intersections11.push_back(temp_intersections[i]);
            }
        }
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((start_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((start_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((end_point - control_points[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((end_point - control_points[degree]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeArcRationalBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::BezierCurveBezierCurveIntersect(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    BezierCurveBezierCurveEquationSystem* equations = intersect_cache->AllocBezierCurveBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(degree0, control_points0, degree1, control_points1, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointBezierCurveIntersect(control_points0[0], degree1, control_points1, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointBezierCurveIntersect(control_points0[degree0], degree1, control_points1, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointBezierCurveIntersect(control_points1[0], degree0, control_points0, distance_epsilon,
            intersect_cache, is_singularity, intersections10);
        is_singularity0 |= is_singularity;
        PointBezierCurveIntersect(control_points1[degree1], degree0, control_points0, distance_epsilon,
            intersect_cache, is_singularity, intersections11);
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((control_points0[0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeBezierCurveBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::BezierCurveRationalBezierCurveIntersect(int degree0, const WGVector2d* control_points0,
    int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    BezierCurveRationalBezierCurveEquationSystem* equations = intersect_cache->AllocBezierCurveRationalBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(degree0, control_points0, degree1, control_points1, weights1, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointRationalBezierCurveIntersect(control_points0[0], degree1, control_points1, weights1, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointRationalBezierCurveIntersect(control_points0[degree0], degree1, control_points1, weights1, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointBezierCurveIntersect(control_points1[0], degree0, control_points0, distance_epsilon,
            intersect_cache, is_singularity, intersections10);
        is_singularity0 |= is_singularity;
        PointBezierCurveIntersect(control_points1[degree1], degree0, control_points0, distance_epsilon,
            intersect_cache, is_singularity, intersections11);
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((control_points0[0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeBezierCurveRationalBezierCurveEquations(equations);
    }
}

void WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveIntersect(int degree0, const WGVector2d* control_points0, const double* weights0,
    int degree1, const WGVector2d* control_points1, const double* weights1, double distance_epsilon, IntersectCache* intersect_cache,
    bool& is_singularity0, bool& is_singularity1, std::vector<CurveCurveIntersection>& intersections) {
    RationalBezierCurveRationalBezierCurveEquationSystem* equations = intersect_cache->AllocRationalBezierCurveRationalBezierCurveEquations();
    CurveCurveSolver* solver = intersect_cache->AllocCurveCurveSolver();
    equations->Update(degree0, control_points0, weights0, degree1, control_points1, weights1, distance_epsilon);
    WSIntervalVector variable(equations->GetVariableCount());
    variable.Set(0, WSInterval(0, 1));
    variable.Set(1, WSInterval(0, 1));
    equations->InitializeIntermediateVariables(&variable);
    PrepareCurveCurveIntersect(equations, solver, &variable, distance_epsilon);
    is_singularity0 = solver->GetIsSingularity0();
    is_singularity1 = solver->GetIsSingularity1();
    if (!is_singularity0 && !is_singularity1) {
        std::vector<PointCurveIntersection> intersections00;
        std::vector<PointCurveIntersection> intersections01;
        std::vector<PointCurveIntersection> intersections10;
        std::vector<PointCurveIntersection> intersections11;
        bool is_singularity;
        PointRationalBezierCurveIntersect(control_points0[0], degree1, control_points1, weights1, distance_epsilon,
            intersect_cache, is_singularity, intersections00);
        is_singularity1 |= is_singularity;
        PointRationalBezierCurveIntersect(control_points0[degree0], degree1, control_points1, weights1, distance_epsilon,
            intersect_cache, is_singularity, intersections01);
        is_singularity1 |= is_singularity;
        PointRationalBezierCurveIntersect(control_points1[0], degree0, control_points0, weights0, distance_epsilon,
            intersect_cache, is_singularity, intersections10);
        is_singularity0 |= is_singularity;
        PointRationalBezierCurveIntersect(control_points1[degree1], degree0, control_points0, weights0, distance_epsilon,
            intersect_cache, is_singularity, intersections11);
        is_singularity0 |= is_singularity;
        if (!is_singularity0 && !is_singularity1) {
            BuildCurveCurveIntersections(equations, solver, &variable, distance_epsilon, intersect_cache, intersections);
            MergeNotSampleIntersections(equations, &variable, intersections);
            double epsilon2 = distance_epsilon * distance_epsilon;
            double end_ts0[4];
            double end_ts1[4];
            int end_tc = 0;
            if ((control_points0[0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections00.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(0));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 0;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[0]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(0));
                intersections10.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 0;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            if ((control_points0[degree0] - control_points1[degree1]).SqrLength() <= epsilon2) {
                intersections01.push_back(PointCurveIntersection(1));
                intersections11.push_back(PointCurveIntersection(1));
                end_ts0[end_tc] = 1;
                end_ts1[end_tc] = 1;
                ++end_tc;
            }
            AdjustIntersectionsWithEndpoints(equations, &variable, intersections00, intersections01,
                intersections10, intersections11, distance_epsilon, intersections);
            RepairByEndTs(end_ts0, end_ts1, end_tc, intersections);
        }
        intersect_cache->FreeCurveCurveSolver(solver);
        intersect_cache->FreeRationalBezierCurveRationalBezierCurveEquations(equations);
    }
}

WGIntersectHelper2d::IntersectCache::IntersectCache() :
    m_first_point_bezier_curve_equations(nullptr),
    m_first_point_rational_bezier_curve_equations(nullptr),
    m_first_line_bezier_curve_equations(nullptr),
    m_first_line_rational_bezier_curve_equations(nullptr),
    m_first_arc_bezier_curve_equations(nullptr),
    m_first_arc_rational_bezier_curve_equations(nullptr),
    m_first_bezier_curve_bezier_curve_equations(nullptr),
    m_first_bezier_curve_rational_bezier_curve_equations(nullptr),
    m_first_rational_bezier_curve_rational_bezier_curve_equations(nullptr),
    m_point_curve_solver(nullptr),
    m_curve_curve_solver(nullptr),
    m_first_pair(nullptr),
    m_first_pair_list(nullptr) {
}

WGIntersectHelper2d::IntersectCache::~IntersectCache() {
    while (m_first_point_bezier_curve_equations) {
        PointBezierCurveEquationSystem* temp = m_first_point_bezier_curve_equations;
        m_first_point_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_point_rational_bezier_curve_equations) {
        PointRationalBezierCurveEquationSystem* temp = m_first_point_rational_bezier_curve_equations;
        m_first_point_rational_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_line_bezier_curve_equations) {
        LineBezierCurveEquationSystem* temp = m_first_line_bezier_curve_equations;
        m_first_line_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_line_rational_bezier_curve_equations) {
        LineRationalBezierCurveEquationSystem* temp = m_first_line_rational_bezier_curve_equations;
        m_first_line_rational_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_arc_bezier_curve_equations) {
        ArcBezierCurveEquationSystem* temp = m_first_arc_bezier_curve_equations;
        m_first_arc_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_arc_rational_bezier_curve_equations) {
        ArcRationalBezierCurveEquationSystem* temp = m_first_arc_rational_bezier_curve_equations;
        m_first_arc_rational_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_bezier_curve_bezier_curve_equations) {
        BezierCurveBezierCurveEquationSystem* temp = m_first_bezier_curve_bezier_curve_equations;
        m_first_bezier_curve_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_bezier_curve_rational_bezier_curve_equations) {
        BezierCurveRationalBezierCurveEquationSystem* temp = m_first_bezier_curve_rational_bezier_curve_equations;
        m_first_bezier_curve_rational_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_first_rational_bezier_curve_rational_bezier_curve_equations) {
        RationalBezierCurveRationalBezierCurveEquationSystem* temp = m_first_rational_bezier_curve_rational_bezier_curve_equations;
        m_first_rational_bezier_curve_rational_bezier_curve_equations = temp->Next;
        delete temp;
    }
    while (m_point_curve_solver) {
        PointCurveSolver* temp = m_point_curve_solver;
        m_point_curve_solver = temp->Next;
        delete temp;
    }
    while (m_curve_curve_solver) {
        CurveCurveSolver* temp = m_curve_curve_solver;
        m_curve_curve_solver = temp->Next;
        delete temp;
    }
    while (m_first_pair_list) {
        CurveCurveFuzzyPairList* temp = m_first_pair_list;
        m_first_pair_list = temp->Next;
        delete temp;
    }
    while (m_first_pair) {
        CurveCurveFuzzyPair* temp = m_first_pair;
        m_first_pair = temp->Next;
        delete temp;
    }
}

WGIntersectHelper2d::PointBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocPointBezierCurveEquations() {
    PointBezierCurveEquationSystem* equations;
    if (m_first_point_bezier_curve_equations) {
        equations = m_first_point_bezier_curve_equations;
        m_first_point_bezier_curve_equations = m_first_point_bezier_curve_equations->Next;
    }
    else {
        equations = new PointBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreePointBezierCurveEquations(PointBezierCurveEquationSystem* equations) {
    equations->Next = m_first_point_bezier_curve_equations;
    m_first_point_bezier_curve_equations = equations;
}

WGIntersectHelper2d::PointRationalBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocPointRationalBezierCurveEquations() {
    PointRationalBezierCurveEquationSystem* equations;
    if (m_first_point_rational_bezier_curve_equations) {
        equations = m_first_point_rational_bezier_curve_equations;
        m_first_point_rational_bezier_curve_equations = m_first_point_rational_bezier_curve_equations->Next;
    }
    else {
        equations = new PointRationalBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreePointRationalBezierCurveEquations(PointRationalBezierCurveEquationSystem* equations) {
    equations->Next = m_first_point_rational_bezier_curve_equations;
    m_first_point_rational_bezier_curve_equations = equations;
}

WGIntersectHelper2d::LineBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocLineBezierCurveEquations() {
    LineBezierCurveEquationSystem* equations;
    if (m_first_line_bezier_curve_equations) {
        equations = m_first_line_bezier_curve_equations;
        m_first_line_bezier_curve_equations = m_first_line_bezier_curve_equations->Next;
    }
    else {
        equations = new LineBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeLineBezierCurveEquations(LineBezierCurveEquationSystem* equations) {
    equations->Next = m_first_line_bezier_curve_equations;
    m_first_line_bezier_curve_equations = equations;
}

WGIntersectHelper2d::LineRationalBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocLineRationalBezierCurveEquations() {
    LineRationalBezierCurveEquationSystem* equations;
    if (m_first_line_rational_bezier_curve_equations) {
        equations = m_first_line_rational_bezier_curve_equations;
        m_first_line_rational_bezier_curve_equations = m_first_line_rational_bezier_curve_equations->Next;
    }
    else {
        equations = new LineRationalBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeLineRationalBezierCurveEquations(LineRationalBezierCurveEquationSystem* equations) {
    equations->Next = m_first_line_rational_bezier_curve_equations;
    m_first_line_rational_bezier_curve_equations = equations;
}

WGIntersectHelper2d::ArcBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocArcBezierCurveEquations() {
    ArcBezierCurveEquationSystem* equations;
    if (m_first_arc_bezier_curve_equations) {
        equations = m_first_arc_bezier_curve_equations;
        m_first_arc_bezier_curve_equations = m_first_arc_bezier_curve_equations->Next;
    }
    else {
        equations = new ArcBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeArcBezierCurveEquations(ArcBezierCurveEquationSystem* equations) {
    equations->Next = m_first_arc_bezier_curve_equations;
    m_first_arc_bezier_curve_equations = equations;
}

WGIntersectHelper2d::ArcRationalBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocArcRationalBezierCurveEquations() {
    ArcRationalBezierCurveEquationSystem* equations;
    if (m_first_arc_rational_bezier_curve_equations) {
        equations = m_first_arc_rational_bezier_curve_equations;
        m_first_arc_rational_bezier_curve_equations = m_first_arc_rational_bezier_curve_equations->Next;
    }
    else {
        equations = new ArcRationalBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeArcRationalBezierCurveEquations(ArcRationalBezierCurveEquationSystem* equations) {
    equations->Next = m_first_arc_rational_bezier_curve_equations;
    m_first_arc_rational_bezier_curve_equations = equations;
}

WGIntersectHelper2d::BezierCurveBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocBezierCurveBezierCurveEquations() {
    BezierCurveBezierCurveEquationSystem* equations;
    if (m_first_bezier_curve_bezier_curve_equations) {
        equations = m_first_bezier_curve_bezier_curve_equations;
        m_first_bezier_curve_bezier_curve_equations = m_first_bezier_curve_bezier_curve_equations->Next;
    }
    else {
        equations = new BezierCurveBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeBezierCurveBezierCurveEquations(BezierCurveBezierCurveEquationSystem* equations) {
    equations->Next = m_first_bezier_curve_bezier_curve_equations;
    m_first_bezier_curve_bezier_curve_equations = equations;
}

WGIntersectHelper2d::BezierCurveRationalBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocBezierCurveRationalBezierCurveEquations() {
    BezierCurveRationalBezierCurveEquationSystem* equations;
    if (m_first_bezier_curve_rational_bezier_curve_equations) {
        equations = m_first_bezier_curve_rational_bezier_curve_equations;
        m_first_bezier_curve_rational_bezier_curve_equations = m_first_bezier_curve_rational_bezier_curve_equations->Next;
    }
    else {
        equations = new BezierCurveRationalBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeBezierCurveRationalBezierCurveEquations(BezierCurveRationalBezierCurveEquationSystem* equations) {
    equations->Next = m_first_bezier_curve_rational_bezier_curve_equations;
    m_first_bezier_curve_rational_bezier_curve_equations = equations;
}

WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveEquationSystem* WGIntersectHelper2d::IntersectCache::AllocRationalBezierCurveRationalBezierCurveEquations() {
    RationalBezierCurveRationalBezierCurveEquationSystem* equations;
    if (m_first_rational_bezier_curve_rational_bezier_curve_equations) {
        equations = m_first_rational_bezier_curve_rational_bezier_curve_equations;
        m_first_rational_bezier_curve_rational_bezier_curve_equations = m_first_rational_bezier_curve_rational_bezier_curve_equations->Next;
    }
    else {
        equations = new RationalBezierCurveRationalBezierCurveEquationSystem();
    }
    return equations;
}

void WGIntersectHelper2d::IntersectCache::FreeRationalBezierCurveRationalBezierCurveEquations(RationalBezierCurveRationalBezierCurveEquationSystem* equations) {
    equations->Next = m_first_rational_bezier_curve_rational_bezier_curve_equations;
    m_first_rational_bezier_curve_rational_bezier_curve_equations = equations;
}

WGIntersectHelper2d::PointCurveSolver* WGIntersectHelper2d::IntersectCache::AllocPointCurveSolver() {
    PointCurveSolver* solver;
    if (m_point_curve_solver) {
        solver = m_point_curve_solver;
        m_point_curve_solver = m_point_curve_solver->Next;
    }
    else {
        solver = new PointCurveSolver();
    }
    return solver;
}
void WGIntersectHelper2d::IntersectCache::FreePointCurveSolver(PointCurveSolver* solver) {
    solver->Next = m_point_curve_solver;
    m_point_curve_solver = solver;
}

WGIntersectHelper2d::CurveCurveSolver* WGIntersectHelper2d::IntersectCache::AllocCurveCurveSolver() {
    CurveCurveSolver* solver;
    if (m_curve_curve_solver) {
        solver = m_curve_curve_solver;
        m_curve_curve_solver = m_curve_curve_solver->Next;
    }
    else {
        solver = new CurveCurveSolver();
    }
    return solver;
}

void WGIntersectHelper2d::IntersectCache::FreeCurveCurveSolver(CurveCurveSolver* solver) {
    solver->Next = m_curve_curve_solver;
    m_curve_curve_solver = solver;
}

WGIntersectHelper2d::CurveCurveFuzzyPairList* WGIntersectHelper2d::IntersectCache::AllocPairList() {
    CurveCurveFuzzyPairList* fuzzy_pair_list;
    if (m_first_pair_list) {
        fuzzy_pair_list = m_first_pair_list;
        m_first_pair_list = fuzzy_pair_list->Next;
    }
    else {
        fuzzy_pair_list = new CurveCurveFuzzyPairList();
        fuzzy_pair_list->TailPair = nullptr;
    }
    return fuzzy_pair_list;
}

void WGIntersectHelper2d::IntersectCache::FreePairList(CurveCurveFuzzyPairList* fuzzy_pair_list) {
    if (fuzzy_pair_list->TailPair) {
        CurveCurveFuzzyPair* fuzzy_pair = fuzzy_pair_list->TailPair;
        do {
            CurveCurveFuzzyPair* temp = fuzzy_pair;
            fuzzy_pair = temp->Next;
            FreePair(temp);
        } while (fuzzy_pair != fuzzy_pair_list->TailPair);
    }
    fuzzy_pair_list->Next = m_first_pair_list;
    m_first_pair_list = fuzzy_pair_list;
}

WGIntersectHelper2d::CurveCurveFuzzyPair* WGIntersectHelper2d::IntersectCache::AllocPair() {
    CurveCurveFuzzyPair* fuzzy_pair;
    if (m_first_pair) {
        fuzzy_pair = m_first_pair;
        m_first_pair = fuzzy_pair->Next;
    }
    else {
        fuzzy_pair = new CurveCurveFuzzyPair();
    }
    return fuzzy_pair;
}

void WGIntersectHelper2d::IntersectCache::FreePair(CurveCurveFuzzyPair* fuzzy_pair) {
    fuzzy_pair->Next = m_first_pair;
    m_first_pair = fuzzy_pair;
}

bool WGIntersectHelper2d::CheckMonotonic(CurveCurveEquationSystem* equations, const WSIntervalVector* variable) {
    WSInterval dx0, dy0, dx1, dy1;
    equations->CalculateD1(0, variable, dx0, dy0);
    equations->CalculateD1(1, variable, dx1, dy1);
    WSInterval d = dx0 * dy1 - dx1 * dy0;
    return d.Min > 0 || d.Max < 0;
}

bool WGIntersectHelper2d::CheckMonotonicOrShortLinear(CurveCurveEquationSystem* equations,
    const WSIntervalVector* variable, double limit_angle, double distance_epsilon) {
    WSInterval dx0, dy0, dx1, dy1;
    equations->CalculateD1(0, variable, dx0, dy0);
    equations->CalculateD1(1, variable, dx1, dy1);
    WSInterval d = dx0 * dy1 - dx1 * dy0;
    if (d.Min > 0 || d.Max < 0) {
        return true;
    }
    double a0 = calculate_tangent_cone_angle(dx0, dy0);
    double a1 = calculate_tangent_cone_angle(dx1, dy1);
    if (a0 <= limit_angle && a1 <= limit_angle) {
        const WSInterval& t0 = variable->Get(0);
        if (t0.Length() != 0) {
            WGVector2d point0, point1;
            equations->CalculateD0(0, t0.Min, point0);
            equations->CalculateD0(0, t0.Max, point1);
            if ((point0 - point1).SqrLength() >= distance_epsilon * distance_epsilon) {
                return false;
            }
        }
        const WSInterval& t1 = variable->Get(1);
        if (t1.Length() != 0) {
            WGVector2d point0, point1;
            equations->CalculateD0(1, t1.Min, point0);
            equations->CalculateD0(1, t1.Max, point1);
            if ((point0 - point1).SqrLength() >= distance_epsilon * distance_epsilon) {
                return false;
            }
        }
        return true;
    }
    return false;
}

void WGIntersectHelper2d::GetFuzzyDomainInfo(CurveCurveEquationSystem* equations, const WSIntervalVector* variable,
    WGVector2d& sample_direction, bool& same_direction) {
    WSInterval dx0, dy0, dx1, dy1;
    equations->CalculateD1(0, variable, dx0, dy0);
    equations->CalculateD1(1, variable, dx1, dy1);
    WGVector2d vt0(dx0.Middle(), dy0.Middle());
    WGVector2d vt1(dx1.Middle(), dy1.Middle());
    vt0.Normalize(0);
    vt1.Normalize(0);
    same_direction = vt0.Dot(vt1) > 0;
    if (same_direction) {
        sample_direction = vt0 + vt1;
    }
    else {
        sample_direction = vt0 - vt1;
    }
    sample_direction.Normalize(0);
    sample_direction = WGVector2d(-sample_direction.Y, sample_direction.X);
}

bool WGIntersectHelper2d::Sample(CurveCurveEquationSystem* equations, const WGVector2d& base_point, int sample_curve_index,
    const WGVector2d& sample_direction, const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance) {
    {
        sample_distance = 0;
        WGVector2d p0 = base_point;
        sample_t = sample_domain.Middle();
        WGVector2d p1;
        equations->CalculateD0(sample_curve_index, sample_t, p1);
        WGVector2d v = p1 - p0;
        double d = v.Length();
        if (d <= epsilon) {
            return true;
        }
        WGVector2d v1;
        equations->CalculateD1(sample_curve_index, sample_t, v1);
        for (int itr_count = 0; itr_count <= 30; ++itr_count) {
            double a = v1.Cross(sample_direction);
            if (a == 0) {
                return false;
            }
            double dt0 = v1.Cross(v) / a;
            double dt1 = sample_direction.Cross(v) / a;
            double t1 = sample_t + dt1;
            if (t1 < sample_domain.Min) {
                t1 = sample_domain.Min;
                double dt = t1 - sample_t;
                dt0 *= dt / dt1;
                dt1 = dt;
            }
            else if (t1 > sample_domain.Max) {
                t1 = sample_domain.Max;
                double dt = t1 - sample_t;
                dt0 *= dt / dt1;
                dt1 = dt;
            }
            if (abs(dt0) <= 1E-16 && abs(dt1) <= 1E-16) {
                break;
            }
            double nt0 = sample_distance + dt0;
            double nt1 = sample_t + dt1;
            p0 = base_point + sample_direction * nt0;
            equations->CalculateD0(sample_curve_index, nt1, p1);
            v = p1 - p0;
            double nd = v.Length();
            bool b = false;
            while (nd > d) {
                dt0 *= 0.5;
                dt1 *= 0.5;
                if (abs(dt0) <= 1E-16 && abs(dt1) <= 1E-16) {
                    b = true;
                    break;
                }
                nt0 = sample_distance + dt0;
                nt1 = sample_t + dt1;
                p0 = base_point + sample_direction * nt0;
                equations->CalculateD0(sample_curve_index, nt1, p1);
                v = p1 - p0;
                nd = v.Length();
            }
            if (b) {
                break;
            }
            sample_distance = nt0;
            sample_t = nt1;
            d = nd;
            if (d <= epsilon) {
                return true;
            }
            equations->CalculateD1(sample_curve_index, sample_t, v1);
        }
    }
    {
        double t0 = sample_domain.Min;
        double t1 = sample_domain.Max;
        WGVector2d p0, p1;
        equations->CalculateD0(sample_curve_index, t0, p0);
        equations->CalculateD0(sample_curve_index, t1, p1);
        WGVector2d v0 = p0 - base_point;
        double d0 = v0.Cross(sample_direction);
        WGVector2d v1 = p1 - base_point;
        double d1 = v1.Cross(sample_direction);
        double a0 = abs(d0);
        double a1 = abs(d1);
        if (a0 < a1) {
            if (a0 <= epsilon) {
                sample_t = t0;
                sample_distance = v0.Dot(sample_direction);
                return true;
            }
        }
        else {
            if (a1 <= epsilon) {
                sample_t = t1;
                sample_distance = v1.Dot(sample_direction);
                return true;
            }
        }
        if (d0 * d1 > 0) {
            return false;
        }
        while (true) {
            double t = (t0 + t1) * 0.5;
            WGVector2d p;
            equations->CalculateD0(sample_curve_index, t, p);
            WGVector2d v = p - base_point;
            double d = v.Cross(sample_direction);
            if (abs(d) <= epsilon) {
                sample_t = t;
                sample_distance = v.Dot(sample_direction);
                return true;
            }
            if (d0 * d < 0) {
                t1 = t;
            }
            else {
                t0 = t;
            }
            if (t1 - t0 <= g_double_epsilon) {
                sample_t = t;
                sample_distance = v.Dot(sample_direction);
                return true;
            }
        }
    }
}

bool WGIntersectHelper2d::Sample(CurveCurveEquationSystem* equations, int base_curve_index, double base_t,
    const WGVector2d& sample_direction, const WSInterval& sample_domain, double epsilon, double& sample_t, double& sample_distance) {
    WGVector2d base_point;
    equations->CalculateD0(base_curve_index, base_t, base_point);
    int sample_curve_index = base_curve_index ^ 1;
    return Sample(equations, base_point, sample_curve_index, sample_direction, sample_domain, epsilon, sample_t, sample_distance);
}

bool WGIntersectHelper2d::MonotonicIterate(CurveCurveEquationSystem* equations, double distance_epsilon,
    const WSInterval& domain0, const WSInterval& domain1, double& t0, double& t1) {
    t0 = domain0.Middle();
    t1 = domain1.Middle();
    WGVector2d point0, point1;
    equations->CalculateD0(0, t0, point0);
    equations->CalculateD0(1, t1, point1);
    WGVector2d vt = point1 - point0;
    double d = vt.Length();
    for (int itr_count = 0; itr_count < 100000; ++itr_count) {
        WGVector2d vt0, vt1;
        equations->CalculateD1(0, t0, vt0);
        equations->CalculateD1(1, t1, vt1);
        double a = vt0.Cross(vt1);
        if (a == 0) {
            return d <= distance_epsilon;
        }
        double dt0 = vt.Cross(vt1) / a;
        double dt1 = vt.Cross(vt0) / a;
        double nt0 = t0 + dt0;
        if (nt0 > domain0.Max) {
            nt0 = domain0.Max;
            dt0 = nt0 - t0;
        }
        else if (nt0 < domain0.Min) {
            nt0 = domain0.Min;
            dt0 = nt0 - t0;
        }
        double nt1 = t1 + dt1;
        if (nt1 > domain1.Max) {
            nt1 = domain1.Max;
            dt1 = nt1 - t1;
        }
        else if (nt1 < domain1.Min) {
            nt1 = domain1.Min;
            dt1 = nt1 - t1;
        }
        double et0, et1;
        if (d <= distance_epsilon) {
            double c = vt0.Length();
            if (c == 0) {
                et0 = 1;
            }
            else {
                et0 = distance_epsilon / c;
                if (et0 > 1) {
                    et0 = 1;
                }
            }
            et0 *= 0.01;
            c = vt1.Length();
            if (c == 0) {
                et1 = 1;
            }
            else {
                et1 = distance_epsilon / c;
                if (et1 > 1) {
                    et1 = 1;
                }
            }
            et1 *= 0.01;
        }
        else {
            et0 = 1E-16;
            et1 = 1E-16;
        }
        if (abs(dt0) <= et0 && abs(dt1) <= et1) {
            return d <= distance_epsilon;
        }
        WGVector2d np0, np1;
        equations->CalculateD0(0, nt0, np0);
        equations->CalculateD0(1, nt1, np1);
        WGVector2d nv = np1 - np0;
        double nd = nv.Length();
        while (nd > d) {
            dt0 *= 0.5;
            dt1 *= 0.5;
            if (abs(dt0) <= et0 && abs(dt1) <= et1) {
                return d <= distance_epsilon;
            }
            nt0 = t0 + dt0;
            nt1 = t1 + dt1;
            equations->CalculateD0(0, nt0, np0);
            equations->CalculateD0(1, nt1, np1);
            nv = np1 - np0;
            nd = nv.Length();
        }
        t0 = nt0;
        t1 = nt1;
        point0 = np0;
        point1 = np1;
        vt = nv;
        d = nd;
    }
    return d <= distance_epsilon;
}

void WGIntersectHelper2d::BuildMonotonicIntersection(CurveCurveEquationSystem* equations, double distance_epsilon,
    const WSInterval& domain0, const WSInterval& domain1, std::vector<CurveCurveIntersection>& intersections) {
    double t0, t1;
    if (MonotonicIterate(equations, distance_epsilon, domain0, domain1, t0, t1)) {
        intersections.push_back(CurveCurveIntersection(t0, t1, false));
    }
}

WGIntersectHelper2d::CurveCurveFuzzyPair* WGIntersectHelper2d::BuildFuzzyPair(CurveCurveEquationSystem* equations, WSIntervalVector* runtime_variable,
    double distance_epsilon, const WSInterval& domain0, const WSInterval& domain1, IntersectCache* intsect_cache) {
    double sample_epsilon = distance_epsilon * 0.001;
    CurveCurveFuzzyPair* fuzzy_pair = intsect_cache->AllocPair();
    runtime_variable->Set(0, domain0);
    runtime_variable->Set(1, domain1);
    GetFuzzyDomainInfo(equations, runtime_variable, fuzzy_pair->SampleDirection, fuzzy_pair->SameDirection);
    if (fuzzy_pair->SameDirection) {
        int sample_count = 0;
        double sample_t;
        double sample_distance;
        if (Sample(equations, 0, domain0.Min, fuzzy_pair->SampleDirection, domain1, sample_epsilon, sample_t, sample_distance)) {
            fuzzy_pair->SamplePoints[0].T0 = domain0.Min;
            fuzzy_pair->SamplePoints[0].T1 = sample_t;
            fuzzy_pair->SamplePoints[0].Distance = sample_distance;
            sample_count = 1;
        }
        else if (Sample(equations, 1, domain1.Min, fuzzy_pair->SampleDirection, domain0, sample_epsilon, sample_t, sample_distance)) {
            fuzzy_pair->SamplePoints[0].T0 = sample_t;
            fuzzy_pair->SamplePoints[0].T1 = domain1.Min;
            fuzzy_pair->SamplePoints[0].Distance = -sample_distance;
            sample_count = 1;
        }
        if (sample_count == 1) {
            if (Sample(equations, 0, domain0.Max, fuzzy_pair->SampleDirection, domain1, sample_epsilon, sample_t, sample_distance)) {
                fuzzy_pair->SamplePoints[1].T0 = domain0.Max;
                fuzzy_pair->SamplePoints[1].T1 = sample_t;
                fuzzy_pair->SamplePoints[1].Distance = sample_distance;
                sample_count = 2;
            }
            else if (Sample(equations, 1, domain1.Max, fuzzy_pair->SampleDirection, domain0, sample_epsilon, sample_t, sample_distance)) {
                fuzzy_pair->SamplePoints[1].T0 = sample_t;
                fuzzy_pair->SamplePoints[1].T1 = domain1.Max;
                fuzzy_pair->SamplePoints[1].Distance = -sample_distance;
                sample_count = 2;
            }
        }
        if (sample_count != 2) {
            intsect_cache->FreePair(fuzzy_pair);
            fuzzy_pair = nullptr;
        }
    }
    else {
        int sample_count = 0;
        double sample_t;
        double sample_distance;
        if (Sample(equations, 0, domain0.Min, fuzzy_pair->SampleDirection, domain1, sample_epsilon, sample_t, sample_distance)) {
            fuzzy_pair->SamplePoints[0].T0 = domain0.Min;
            fuzzy_pair->SamplePoints[0].T1 = sample_t;
            fuzzy_pair->SamplePoints[0].Distance = sample_distance;
            sample_count = 1;
        }
        else if (Sample(equations, 1, domain1.Max, fuzzy_pair->SampleDirection, domain0, sample_epsilon, sample_t, sample_distance)) {
            fuzzy_pair->SamplePoints[0].T0 = sample_t;
            fuzzy_pair->SamplePoints[0].T1 = domain1.Max;
            fuzzy_pair->SamplePoints[0].Distance = -sample_distance;
            sample_count = 1;
        }
        if (sample_count == 1) {
            if (Sample(equations, 0, domain0.Max, fuzzy_pair->SampleDirection, domain1, sample_epsilon, sample_t, sample_distance)) {
                fuzzy_pair->SamplePoints[1].T0 = domain0.Max;
                fuzzy_pair->SamplePoints[1].T1 = sample_t;
                fuzzy_pair->SamplePoints[1].Distance = sample_distance;
                sample_count = 2;
            }
            else if (Sample(equations, 1, domain1.Min, fuzzy_pair->SampleDirection, domain0, sample_epsilon, sample_t, sample_distance)) {
                fuzzy_pair->SamplePoints[1].T0 = sample_t;
                fuzzy_pair->SamplePoints[1].T1 = domain1.Min;
                fuzzy_pair->SamplePoints[1].Distance = -sample_distance;
                sample_count = 2;
            }
        }
        if (sample_count != 2) {
            intsect_cache->FreePair(fuzzy_pair);
            fuzzy_pair = nullptr;
        }
    }
    if (!fuzzy_pair) {
        return nullptr;
    }
    fuzzy_pair->Domain0 = domain0;
    fuzzy_pair->Domain1 = domain1;
    fuzzy_pair->Length = -1;
    return fuzzy_pair;
}

WGIntersectHelper2d::CurveCurveFuzzyPair* WGIntersectHelper2d::SplitFuzzyPair(CurveCurveFuzzyPair* fuzzy_pair,
    double split_t0, double split_t1, double distance, IntersectCache* intsect_cache) {
    CurveCurveFuzzyPair* fuzzy_pair2 = intsect_cache->AllocPair();
    fuzzy_pair2->SameDirection = fuzzy_pair->SameDirection;
    fuzzy_pair2->SampleDirection = fuzzy_pair->SampleDirection;
    fuzzy_pair2->Domain0 = fuzzy_pair->Domain0;
    fuzzy_pair2->Domain1 = fuzzy_pair->Domain1;
    fuzzy_pair2->SamplePoints[1] = fuzzy_pair->SamplePoints[1];
    fuzzy_pair->SamplePoints[1].T0 = split_t0;
    fuzzy_pair->SamplePoints[1].T1 = split_t1;
    fuzzy_pair->SamplePoints[1].Distance = distance;
    fuzzy_pair2->SamplePoints[0].T0 = split_t0;
    fuzzy_pair2->SamplePoints[0].T1 = split_t1;
    fuzzy_pair2->SamplePoints[0].Distance = distance;
    if (fuzzy_pair->SameDirection) {
        fuzzy_pair->Domain0.Max = split_t0;
        fuzzy_pair->Domain1.Max = split_t1;
        fuzzy_pair2->Domain0.Min = split_t0;
        fuzzy_pair2->Domain1.Min = split_t1;
    }
    else {
        fuzzy_pair->Domain0.Max = split_t0;
        fuzzy_pair->Domain1.Min = split_t1;
        fuzzy_pair2->Domain0.Min = split_t0;
        fuzzy_pair2->Domain1.Max = split_t1;
    }
    fuzzy_pair->Length = -1;
    fuzzy_pair2->Length = -1;
    return fuzzy_pair2;
}

WGIntersectHelper2d::CurveCurveFuzzyPairList* WGIntersectHelper2d::NewFuzzyPairList(
    IntersectCache* intsect_cache, CurveCurveFuzzyPair* fuzzy_pair) {
    CurveCurveFuzzyPairList* fuzzy_pair_list = intsect_cache->AllocPairList();
    fuzzy_pair_list->Domain0 = fuzzy_pair->Domain0;
    fuzzy_pair_list->Domain1 = fuzzy_pair->Domain1;
    fuzzy_pair_list->SameDirection = fuzzy_pair->SameDirection;
    fuzzy_pair_list->TailPair = fuzzy_pair;
    fuzzy_pair->Next = fuzzy_pair;
    return fuzzy_pair_list;
}

bool WGIntersectHelper2d::MergeFuzzyPairList(CurveCurveFuzzyPairList* fuzzy_pair_list, CurveCurveFuzzyPair* fuzzy_pair) {
    if (fuzzy_pair->SameDirection != fuzzy_pair_list->SameDirection) {
        return false;
    }
    if (fuzzy_pair->Domain0.Min > fuzzy_pair_list->Domain0.Max ||
        fuzzy_pair->Domain0.Max < fuzzy_pair_list->Domain0.Min ||
        fuzzy_pair->Domain1.Min > fuzzy_pair_list->Domain1.Max ||
        fuzzy_pair->Domain1.Max < fuzzy_pair_list->Domain1.Min) {
        return false;
    }
    fuzzy_pair->Next = fuzzy_pair_list->TailPair->Next;
    fuzzy_pair_list->TailPair->Next = fuzzy_pair;
    fuzzy_pair_list->Domain0.Merge(fuzzy_pair->Domain0);
    fuzzy_pair_list->Domain1.Merge(fuzzy_pair->Domain1);
    return true;
}

bool WGIntersectHelper2d::MergeFuzzyPairList(CurveCurveFuzzyPairList* dst_fuzzy_pair_list, CurveCurveFuzzyPairList* src_fuzzy_pair_list) {
    if (src_fuzzy_pair_list->SameDirection != dst_fuzzy_pair_list->SameDirection) {
        return false;
    }
    if (src_fuzzy_pair_list->Domain0.Min > dst_fuzzy_pair_list->Domain0.Max ||
        src_fuzzy_pair_list->Domain0.Max < dst_fuzzy_pair_list->Domain0.Min ||
        src_fuzzy_pair_list->Domain1.Min > dst_fuzzy_pair_list->Domain1.Max ||
        src_fuzzy_pair_list->Domain1.Max < dst_fuzzy_pair_list->Domain1.Min) {
        return false;
    }
    CurveCurveFuzzyPair* first_src_pair = src_fuzzy_pair_list->TailPair->Next;
    src_fuzzy_pair_list->TailPair->Next = dst_fuzzy_pair_list->TailPair->Next;
    dst_fuzzy_pair_list->TailPair->Next = first_src_pair;
    src_fuzzy_pair_list->TailPair = nullptr;
    dst_fuzzy_pair_list->Domain0.Merge(src_fuzzy_pair_list->Domain0);
    dst_fuzzy_pair_list->Domain1.Merge(src_fuzzy_pair_list->Domain1);
    return true;
}

void WGIntersectHelper2d::BuildFuzzyIntersection(CurveCurveEquationSystem* equations, CurveCurveSolver* solver,
    WSIntervalVector* runtime_variable, double distance_epsilon, CurveCurveFuzzyPair* first_fuzzy_pair, 
    IntersectCache* intsect_cache, std::vector<CurveCurveIntersection>& intersections) {
    double sample_epsilon = distance_epsilon * 0.001;
    CurveCurveFuzzyPairList* first_fuzzy_pair_list = nullptr;
    while (true) {
        if (first_fuzzy_pair) {
            CurveCurveFuzzyPair* fuzzy_pair = first_fuzzy_pair;
            first_fuzzy_pair = fuzzy_pair->Next;
            int n = 0;
            if (abs(fuzzy_pair->SamplePoints[0].Distance) <= distance_epsilon) {
                n |= 1;
            }
            if (abs(fuzzy_pair->SamplePoints[1].Distance) <= distance_epsilon) {
                n |= 2;
            }
            if (n == 3) {
                CurveCurveFuzzyPairList* fuzzy_pair_list = first_fuzzy_pair_list;
                while (fuzzy_pair_list) {
                    if (MergeFuzzyPairList(fuzzy_pair_list, fuzzy_pair)) {
                        break;
                    }
                    fuzzy_pair_list = fuzzy_pair_list->Next;
                }
                if (fuzzy_pair_list) {
                    CurveCurveFuzzyPairList* prev_fuzzy_pair_list = nullptr;
                    CurveCurveFuzzyPairList* curr_fuzzy_pair_list = first_fuzzy_pair_list;
                    while (curr_fuzzy_pair_list) {
                        if (fuzzy_pair_list != curr_fuzzy_pair_list) {
                            if (MergeFuzzyPairList(fuzzy_pair_list, curr_fuzzy_pair_list)) {
                                CurveCurveFuzzyPairList* next_fuzzy_pair_list = curr_fuzzy_pair_list->Next;
                                if (prev_fuzzy_pair_list) {
                                    prev_fuzzy_pair_list->Next = next_fuzzy_pair_list;
                                }
                                else {
                                    first_fuzzy_pair_list = next_fuzzy_pair_list;
                                }
                                intsect_cache->FreePairList(curr_fuzzy_pair_list);
                                curr_fuzzy_pair_list = next_fuzzy_pair_list;
                                continue;
                            }
                        }
                        prev_fuzzy_pair_list = curr_fuzzy_pair_list;
                        curr_fuzzy_pair_list = prev_fuzzy_pair_list->Next;
                    }
                }
                else {
                    fuzzy_pair_list = NewFuzzyPairList(intsect_cache, fuzzy_pair);
                    fuzzy_pair_list->Next = first_fuzzy_pair_list;
                    first_fuzzy_pair_list = fuzzy_pair_list;
                }
                continue;
            }
            if (fuzzy_pair->Length < 0) {
                WGVector2d point0;
                equations->CalculateD0(0, fuzzy_pair->SamplePoints[0].T0, point0);
                WGVector2d point1;
                equations->CalculateD0(0, fuzzy_pair->SamplePoints[1].T0, point1);
                fuzzy_pair->Length = (point1 - point0).Length();
            }
            if (fuzzy_pair->Length <= distance_epsilon) {
                if ((n & 1) != 0) {
                    intersections.push_back(CurveCurveIntersection(fuzzy_pair->SamplePoints[0].T0, fuzzy_pair->SamplePoints[0].T1, true));
                }
                if ((n & 2) != 0) {
                    intersections.push_back(CurveCurveIntersection(fuzzy_pair->SamplePoints[1].T0, fuzzy_pair->SamplePoints[1].T1, true));
                }
                continue;
            }
            runtime_variable->Set(0, WSInterval(fuzzy_pair->SamplePoints[0].T0, fuzzy_pair->SamplePoints[1].T0));
            if (fuzzy_pair->SameDirection) {
                runtime_variable->Set(1, WSInterval(fuzzy_pair->SamplePoints[0].T1, fuzzy_pair->SamplePoints[1].T1));
            }
            else {
                runtime_variable->Set(1, WSInterval(fuzzy_pair->SamplePoints[1].T1, fuzzy_pair->SamplePoints[0].T1));
            }
            equations->InitializeIntermediateVariables(runtime_variable);
            equations->SetAlgebraEnable(runtime_variable, distance_epsilon * 0.1);
            WSIterateResult iterate_result1 = solver->Iterate(runtime_variable, 0.5, 0.1);
            if (iterate_result1 == WSIterateResult::NoRoot) {
                if ((n & 1) != 0) {
                    WGVector2d vt0, vt1;
                    equations->CalculateD1(0, fuzzy_pair->SamplePoints[0].T0, vt0);
                    equations->CalculateD1(1, fuzzy_pair->SamplePoints[0].T1, vt1);
                    if (!fuzzy_pair->SameDirection) {
                        vt1 = -vt1;
                    }
                    if (vt0.Normalize(g_double_epsilon) <= g_double_epsilon || vt1.Normalize(g_double_epsilon) <= g_double_epsilon) {
                        intsect_cache->FreePair(fuzzy_pair);
                        continue;
                    }
                    if ((vt1 - vt0).Dot(fuzzy_pair->SampleDirection) * fuzzy_pair->SamplePoints[0].Distance >= 0) {
                        intsect_cache->FreePair(fuzzy_pair);
                        continue;
                    }
                }
                else if ((n & 2) != 0) {
                    WGVector2d vt0, vt1;
                    equations->CalculateD1(0, fuzzy_pair->SamplePoints[1].T0, vt0);
                    equations->CalculateD1(1, fuzzy_pair->SamplePoints[1].T1, vt1);
                    if (!fuzzy_pair->SameDirection) {
                        vt1 = -vt1;
                    }
                    if (vt0.Normalize(g_double_epsilon) <= g_double_epsilon || vt1.Normalize(g_double_epsilon) <= g_double_epsilon) {
                        intsect_cache->FreePair(fuzzy_pair);
                        continue;
                    }
                    if ((vt1 - vt0).Dot(fuzzy_pair->SampleDirection) * fuzzy_pair->SamplePoints[1].Distance <= 0) {
                        intsect_cache->FreePair(fuzzy_pair);
                        continue;
                    }
                }
                else {
                    intsect_cache->FreePair(fuzzy_pair);
                    continue;
                }
            }
            if (iterate_result1 == WSIterateResult::ClearRoot) {
                double t0 = runtime_variable->Get(0).Middle();
                WSInterval domain1;
                if (fuzzy_pair->SameDirection) {
                    domain1 = WSInterval(fuzzy_pair->SamplePoints[0].T1, fuzzy_pair->SamplePoints[1].T1);
                }
                else {
                    domain1 = WSInterval(fuzzy_pair->SamplePoints[1].T1, fuzzy_pair->SamplePoints[0].T1);
                }
                double t1;
                double d;
                Sample(equations, 0, t0, fuzzy_pair->SampleDirection, domain1, sample_epsilon, t1, d);
                if ((n & 1) != 0) {
                    fuzzy_pair->SamplePoints[1].T0 = t0;
                    fuzzy_pair->SamplePoints[1].T1 = t1;
                    fuzzy_pair->SamplePoints[1].Distance = d;
                    fuzzy_pair->Length = -1;
                    fuzzy_pair->Domain0.Max = t0;
                    if (fuzzy_pair->SameDirection) {
                        fuzzy_pair->Domain1.Max = t1;
                    }
                    else {
                        fuzzy_pair->Domain1.Min = t1;
                    }
                }
                else if ((n & 2) != 0) {
                    fuzzy_pair->SamplePoints[0].T0 = t0;
                    fuzzy_pair->SamplePoints[0].T1 = t1;
                    fuzzy_pair->SamplePoints[0].Distance = d;
                    fuzzy_pair->Length = -1;
                    fuzzy_pair->Domain0.Min = t0;
                    if (fuzzy_pair->SameDirection) {
                        fuzzy_pair->Domain1.Min = t1;
                    }
                    else {
                        fuzzy_pair->Domain1.Max = t1;
                    }
                }
                continue;
            }
            double t0 = (fuzzy_pair->SamplePoints[0].T0 + fuzzy_pair->SamplePoints[1].T0) * 0.5;
            WSInterval domain1;
            if (fuzzy_pair->SameDirection) {
                domain1 = WSInterval(fuzzy_pair->SamplePoints[0].T1, fuzzy_pair->SamplePoints[1].T1);
            }
            else {
                domain1 = WSInterval(fuzzy_pair->SamplePoints[1].T1, fuzzy_pair->SamplePoints[0].T1);
            }
            double t1;
            double d;
            Sample(equations, 0, t0, fuzzy_pair->SampleDirection, domain1, sample_epsilon, t1, d);
            CurveCurveFuzzyPair* fuzzy_pair2 = SplitFuzzyPair(fuzzy_pair, t0, t1, d, intsect_cache);
            fuzzy_pair2->Next = first_fuzzy_pair;
            fuzzy_pair->Next = fuzzy_pair2;
            first_fuzzy_pair = fuzzy_pair;
            continue;
        }
        if (first_fuzzy_pair_list) {
            CurveCurveFuzzyPairList* fuzzy_pair_list = first_fuzzy_pair_list;
            first_fuzzy_pair_list = fuzzy_pair_list->Next;
            int sample_count = equations->GetSampleCount();
            assert(sample_count > 1);
            double delta_t = fuzzy_pair_list->Domain0.Length() / (sample_count - 1);
            CurveCurveFuzzyPair* first_fuzzy_pair = fuzzy_pair_list->TailPair->Next;
            fuzzy_pair_list->TailPair->Next = nullptr;
            double split_t0, split_t1;
            bool splitting = false;
            CurveCurveFuzzyPair* prev_fuzzy_pair = nullptr;
            CurveCurveFuzzyPair* curr_fuzzy_pair = first_fuzzy_pair;
            while (curr_fuzzy_pair) {
                if (curr_fuzzy_pair->SamplePoints[1].T0 - curr_fuzzy_pair->SamplePoints[0].T0 <= delta_t) {
                    curr_fuzzy_pair = curr_fuzzy_pair->Next;
                }
                else {
                    split_t0 = (curr_fuzzy_pair->SamplePoints[1].T0 + curr_fuzzy_pair->SamplePoints[0].T0) * 0.5;
                    WSInterval domain1;
                    if (curr_fuzzy_pair->SameDirection) {
                        domain1 = WSInterval(curr_fuzzy_pair->SamplePoints[0].T1, curr_fuzzy_pair->SamplePoints[1].T1);
                    }
                    else {
                        domain1 = WSInterval(curr_fuzzy_pair->SamplePoints[1].T1, curr_fuzzy_pair->SamplePoints[0].T1);
                    }
                    double d;
                    Sample(equations, 0, split_t0, curr_fuzzy_pair->SampleDirection, domain1, sample_epsilon, split_t1, d);
                    if (abs(d) > distance_epsilon) {
                        if (prev_fuzzy_pair) {
                            prev_fuzzy_pair->Next = curr_fuzzy_pair->Next;
                        }
                        else {
                            first_fuzzy_pair = curr_fuzzy_pair->Next;
                        }
                        CurveCurveFuzzyPair* fuzzy_pair2 = SplitFuzzyPair(curr_fuzzy_pair, split_t0, split_t1, d, intsect_cache);
                        fuzzy_pair2->Next = first_fuzzy_pair;
                        curr_fuzzy_pair->Next = fuzzy_pair2;
                        first_fuzzy_pair = curr_fuzzy_pair;
                        splitting = true;
                        break;
                    }
                    else {
                        CurveCurveFuzzyPair* fuzzy_pair2 = SplitFuzzyPair(curr_fuzzy_pair, split_t0, split_t1, d, intsect_cache);
                        fuzzy_pair2->Next = nullptr;
                        fuzzy_pair_list->TailPair->Next = fuzzy_pair2;
                        fuzzy_pair_list->TailPair = fuzzy_pair2;
                    }
                }
            }
            if (splitting) {
                fuzzy_pair_list->TailPair = nullptr;
                CurveCurveFuzzyPairList* fuzzy_pair_list2 = intsect_cache->AllocPairList();
                fuzzy_pair_list2->SameDirection = fuzzy_pair_list->SameDirection;
                fuzzy_pair_list2->TailPair = nullptr;
                while (first_fuzzy_pair) {
                    curr_fuzzy_pair = first_fuzzy_pair;
                    first_fuzzy_pair = curr_fuzzy_pair->Next;
                    CurveCurveFuzzyPairList* selected_fuzzy_pair_list;
                    if ((curr_fuzzy_pair->SamplePoints[0].T0 + curr_fuzzy_pair->SamplePoints[1].T0) * 0.5 < split_t0) {
                        selected_fuzzy_pair_list = fuzzy_pair_list;
                    }
                    else {
                        selected_fuzzy_pair_list = fuzzy_pair_list2;
                    }
                    if (selected_fuzzy_pair_list->TailPair) {
                        curr_fuzzy_pair->Next = selected_fuzzy_pair_list->TailPair->Next;
                        selected_fuzzy_pair_list->TailPair->Next = curr_fuzzy_pair;
                        selected_fuzzy_pair_list->Domain0.Merge(curr_fuzzy_pair->Domain0);
                        selected_fuzzy_pair_list->Domain1.Merge(curr_fuzzy_pair->Domain1);
                    }
                    else {
                        selected_fuzzy_pair_list->Domain0 = curr_fuzzy_pair->Domain0;
                        selected_fuzzy_pair_list->Domain1 = curr_fuzzy_pair->Domain1;
                        selected_fuzzy_pair_list->TailPair = curr_fuzzy_pair;
                        curr_fuzzy_pair->Next = curr_fuzzy_pair;
                    }
                }
                fuzzy_pair_list2->Next = first_fuzzy_pair_list;
                fuzzy_pair_list->Next = fuzzy_pair_list2;
                first_fuzzy_pair_list = fuzzy_pair_list;
            }
            else {
                CurveCurveFuzzyPair* curr_fuzzy_pair = first_fuzzy_pair;
                double t00 = curr_fuzzy_pair->SamplePoints[0].T0;
                double t01 = curr_fuzzy_pair->SamplePoints[0].T1;
                double t10 = curr_fuzzy_pair->SamplePoints[1].T0;
                double t11 = curr_fuzzy_pair->SamplePoints[1].T1;
                curr_fuzzy_pair = curr_fuzzy_pair->Next;
                while (curr_fuzzy_pair) {
                    if (curr_fuzzy_pair->SamplePoints[0].T0 < t00) {
                        t00 = curr_fuzzy_pair->SamplePoints[0].T0;
                        t01 = curr_fuzzy_pair->SamplePoints[0].T1;
                    }
                    if (curr_fuzzy_pair->SamplePoints[1].T0 > t10) {
                        t10 = curr_fuzzy_pair->SamplePoints[1].T0;
                        t11 = curr_fuzzy_pair->SamplePoints[1].T1;
                    }
                    curr_fuzzy_pair = curr_fuzzy_pair->Next;
                }
                intersections.push_back(CurveCurveIntersection(t00, t01, true, t10, t11, true));
                fuzzy_pair_list->TailPair->Next = first_fuzzy_pair;
                intsect_cache->FreePairList(fuzzy_pair_list);
            }
            continue;
        }
        break;
    }
}

void WGIntersectHelper2d::SortTs(double* ts0, double* ts1, int* indices, int tc) {
    //todo 以后改为更高效的排序方式
    for (int i = 0; i < tc; ++i) {
        for (int j = i + 1; j < tc; ++j) {
            if (ts0[indices[i]] > ts0[indices[j]] ||
                (ts0[indices[i]] == ts0[indices[j]] && ts1[indices[i]] > ts1[indices[j]])) {
                int t = indices[i];
                indices[i] = indices[j];
                indices[j] = t;
            }
        }
    }
}

void WGIntersectHelper2d::RemoveSameTs(double* ts0, double* ts1, bool* is_samples, int& tc) {
    int n = 0;
    for (int i = 0; i < tc; ++i) {
        bool b = true;
        for (int j = i + 1; j < tc; ++j) {
            if (ts0[i] == ts0[j] && ts1[i] == ts1[j]) {
                is_samples[j] &= is_samples[i];
                b = false;
                break;
            }
        }
        if (b) {
            if (i != n) {
                ts0[n] = ts0[i];
                ts1[n] = ts1[i];
                is_samples[n] = is_samples[i];
            }
            ++n;
        }
    }
    tc = n;
}

void WGIntersectHelper2d::RepairCross(double* ts0, double* ts1, bool* is_samples, WGVector2d* points0, WGVector2d* points1, int& tc, double distance_epsilon) {
    double epsilon2 = distance_epsilon * distance_epsilon;
    int i = 0;
    while (i < tc - 1) {
        bool b = true;
        int j = i + 1;
        while (j < tc) {
            if (ts0[j] != ts0[i] && ts1[j] != ts1[i]) {
                WGVector2d vt0 = ts0[i] < ts0[j] ? points0[j] - points0[i] : points0[i] - points0[j];
                WGVector2d vt1 = ts1[i] < ts1[j] ? points1[j] - points1[i] : points1[i] - points1[j];
                double d = vt0.Dot(vt1);
                if (d != 0) {
                    if ((d > 0) != ((ts0[j] - ts0[i]) * (ts1[j] - ts1[i]) > 0)) {
                        double t = ts1[i];
                        ts1[i] = ts1[j];
                        ts1[j] = t;
                        WGVector2d p = points1[i];
                        points1[i] = points1[j];
                        points1[j] = p;
                        is_samples[i] = true;
                        is_samples[j] = true;
                        if ((ts0[j] == 0 || ts0[j] == 1 || ts1[j] == 0 || ts1[j] == 1) && (points1[j] - points0[j]).SqrLength() > epsilon2) {
                            --tc;
                            if (tc > j) {
                                ts0[j] = ts0[tc];
                                ts1[j] = ts1[tc];
                                is_samples[j] = is_samples[tc];
                                points0[j] = points0[tc];
                                points1[j] = points1[tc];
                            }
                        }
                        if ((ts0[i] == 0 || ts0[i] == 1 || ts1[i] == 0 || ts1[i] == 1) && (points1[i] - points0[i]).SqrLength() > epsilon2) {
                            --tc;
                            if (tc > i) {
                                ts0[i] = ts0[tc];
                                ts1[i] = ts1[tc];
                                is_samples[i] = is_samples[tc];
                                points0[i] = points0[tc];
                                points1[i] = points1[tc];
                            }
                        }
                        b = false;
                        break;
                    }
                }
            }
            ++j;
        }
        if (b) {
            ++i;
        }
    }
}

void WGIntersectHelper2d::MarkNotSample(double* ts0, double* ts1, bool* is_samples, int tc, CurveCurveIntersection* intersections, int intersection_count) {
    for (int i = 0; i < tc; ++i) {
        if (!is_samples[i]) {
            for (int j = 0; j < intersection_count; ++j) {
                if (intersections[j].PointCount == 1) {
                    if (intersections[j].Ts0[0] == ts0[i] && intersections[j].Ts1[0] == ts1[i]) {
                        intersections[j].IsSamples[0] = false;
                    }
                }
                else {
                    if (intersections[j].Ts0[0] == ts0[i] && intersections[j].Ts1[0] == ts1[i]) {
                        intersections[j].IsSamples[0] = false;
                    }
                    if (intersections[j].Ts0[1] == ts0[i] && intersections[j].Ts1[1] == ts1[i]) {
                        intersections[j].IsSamples[1] = false;
                    }
                }
            }
        }
    }
}

int WGIntersectHelper2d::MergeOverlapIntersections(CurveCurveIntersection* intersections, int intersection_count) {
    int n = 0;
    for (int i = 0; i < intersection_count; ++i) {
        bool b = true;
        if (intersections[i].PointCount == 2) {
            for (int j = i + 1; j < intersection_count; ++j) {
                if (intersections[j].PointCount != 2) {
                    continue;
                }
                if (intersections[i].Ts0[0] == intersections[j].Ts0[0] && intersections[i].Ts1[0] == intersections[j].Ts1[0]) {
                    if (!intersections[i].IsSamples[0] || !intersections[j].IsSamples[0]) {
                        continue;
                    }
                    double d = abs(intersections[j].Ts0[1] - intersections[i].Ts0[1]);
                    if (d < abs(intersections[j].Ts0[1] - intersections[j].Ts0[0]) || d < abs(intersections[i].Ts0[1] - intersections[i].Ts0[0])) {
                        continue;
                    }
                    d = abs(intersections[j].Ts1[1] - intersections[i].Ts1[1]);
                    if (d < abs(intersections[j].Ts1[1] - intersections[j].Ts1[0]) || d < abs(intersections[i].Ts1[1] - intersections[i].Ts1[0])) {
                        continue;
                    }
                    if (intersections[i].Ts0[1] < intersections[j].Ts0[1]) {
                        intersections[j].Ts0[0] = intersections[i].Ts0[1];
                        intersections[j].Ts1[0] = intersections[i].Ts1[1];
                        intersections[j].IsSamples[0] = intersections[i].IsSamples[1];
                    }
                    else {
                        intersections[j].Ts0[0] = intersections[j].Ts0[1];
                        intersections[j].Ts1[0] = intersections[j].Ts1[1];
                        intersections[j].IsSamples[0] = intersections[j].IsSamples[1];
                        intersections[j].Ts0[1] = intersections[i].Ts0[1];
                        intersections[j].Ts1[1] = intersections[i].Ts1[1];
                        intersections[j].IsSamples[1] = intersections[i].IsSamples[1];
                    }
                    b = false;
                    break;
                }
                if (intersections[i].Ts0[1] == intersections[j].Ts0[1] && intersections[i].Ts1[1] == intersections[j].Ts1[1]) {
                    if (!intersections[i].IsSamples[1] || !intersections[j].IsSamples[1]) {
                        continue;
                    }
                    double d = abs(intersections[j].Ts0[0] - intersections[i].Ts0[0]);
                    if (d < abs(intersections[j].Ts0[1] - intersections[j].Ts0[0]) || d < abs(intersections[i].Ts0[1] - intersections[i].Ts0[0])) {
                        continue;
                    }
                    d = abs(intersections[j].Ts1[0] - intersections[i].Ts1[0]);
                    if (d < abs(intersections[j].Ts1[1] - intersections[j].Ts1[0]) || d < abs(intersections[i].Ts1[1] - intersections[i].Ts1[0])) {
                        continue;
                    }
                    if (intersections[i].Ts0[0] > intersections[j].Ts0[0]) {
                        intersections[j].Ts0[1] = intersections[i].Ts0[0];
                        intersections[j].Ts1[1] = intersections[i].Ts1[0];
                        intersections[j].IsSamples[1] = intersections[i].IsSamples[0];
                    }
                    else {
                        intersections[j].Ts0[1] = intersections[j].Ts0[0];
                        intersections[j].Ts1[1] = intersections[j].Ts1[0];
                        intersections[j].IsSamples[1] = intersections[j].IsSamples[0];
                        intersections[j].Ts0[0] = intersections[i].Ts0[0];
                        intersections[j].Ts1[0] = intersections[i].Ts1[0];
                        intersections[j].IsSamples[0] = intersections[i].IsSamples[0];
                    }
                    b = false;
                    break;
                }
                if (intersections[i].Ts0[0] == intersections[j].Ts0[1] && intersections[i].Ts1[0] == intersections[j].Ts1[1]) {
                    if (!intersections[i].IsSamples[0] || !intersections[j].IsSamples[1]) {
                        continue;
                    }
                    double d = abs(intersections[j].Ts0[0] - intersections[i].Ts0[1]);
                    if (d < abs(intersections[j].Ts0[1] - intersections[j].Ts0[0]) || d < abs(intersections[i].Ts0[1] - intersections[i].Ts0[0])) {
                        continue;
                    }
                    d = abs(intersections[j].Ts1[0] - intersections[i].Ts1[1]);
                    if (d < abs(intersections[j].Ts1[1] - intersections[j].Ts1[0]) || d < abs(intersections[i].Ts1[1] - intersections[i].Ts1[0])) {
                        continue;
                    }
                    intersections[j].Ts0[1] = intersections[i].Ts0[1];
                    intersections[j].Ts1[1] = intersections[i].Ts1[1];
                    intersections[j].IsSamples[1] = intersections[i].IsSamples[1];
                    b = false;
                    break;
                }
                if (intersections[i].Ts0[1] == intersections[j].Ts0[0] && intersections[i].Ts1[1] == intersections[j].Ts1[0]) {
                    if (!intersections[i].IsSamples[1] || !intersections[j].IsSamples[0]) {
                        continue;
                    }
                    double d = abs(intersections[j].Ts0[1] - intersections[i].Ts0[0]);
                    if (d < abs(intersections[j].Ts0[1] - intersections[j].Ts0[0]) || d < abs(intersections[i].Ts0[1] - intersections[i].Ts0[0])) {
                        continue;
                    }
                    d = abs(intersections[j].Ts1[1] - intersections[i].Ts1[0]);
                    if (d < abs(intersections[j].Ts1[1] - intersections[j].Ts1[0]) || d < abs(intersections[i].Ts1[1] - intersections[i].Ts1[0])) {
                        continue;
                    }
                    intersections[j].Ts0[0] = intersections[i].Ts0[0];
                    intersections[j].Ts1[0] = intersections[i].Ts1[0];
                    intersections[j].IsSamples[0] = intersections[i].IsSamples[0];
                    b = false;
                    break;
                }
            }
        }
        if (b) {
            if (i != n) {
                intersections[n] = intersections[i];
            }
            ++n;
        }
    }
    int i = 0;
    while (i < n) {
        bool b = false;
        if (intersections[i].PointCount == 2) {
            for (int j = 0; j < i; ++j) {
                if (intersections[j].PointCount == 2 &&
                    intersections[i].Ts0[0] == intersections[j].Ts0[0] && intersections[i].Ts1[0] == intersections[j].Ts1[0] &&
                    intersections[i].Ts0[1] == intersections[j].Ts0[1] && intersections[i].Ts1[1] == intersections[j].Ts1[1]) {
                    b = true;
                    break;
                }
            }
        }
        if (b) {
            --n;
            if (n > i) {
                intersections[i] = intersections[n];
            }
        }
        else {
            ++i;
        }
    }
    return n;
}

void WGIntersectHelper2d::RepairByEndTs(double* end_ts0, double* end_ts1, int end_tc, std::vector<CurveCurveIntersection>& intersections) {
    for (int i = 0; i < end_tc; ++i) {
        bool b = true;
        for (int j = 0; j < (int)intersections.size(); ++j) {
            const CurveCurveIntersection& intersection = intersections.at(j);
            if (intersection.PointCount == 1) {
                if (intersection.Ts0[0] == end_ts0[i] && intersection.Ts1[0] == end_ts1[i]) {
                    b = false;
                    break;
                }
            }
            else {
                if (intersection.Ts0[0] == end_ts0[i] && intersection.Ts1[0] == end_ts1[i]) {
                    b = false;
                    break;
                }
                if (intersection.Ts0[1] == end_ts0[i] && intersection.Ts1[1] == end_ts1[i]) {
                    b = false;
                    break;
                }
            }
        }
        if (b) {
            intersections.push_back(CurveCurveIntersection(end_ts0[i], end_ts1[i], true));
        }
    }
}

void WGIntersectHelper2d::RepairByEndTs(double* end_ts0, double* end_ts1, int end_tc, CurveCurveIntersection* intersections, int& intersection_count) {
    for (int i = 0; i < end_tc; ++i) {
        bool b = true;
        for (int j = 0; j < intersection_count; ++j) {
            const CurveCurveIntersection& intersection = intersections[j];
            if (intersection.PointCount == 1) {
                if (intersection.Ts0[0] == end_ts0[i] && intersection.Ts1[0] == end_ts1[i]) {
                    b = false;
                    break;
                }
            }
            else {
                if (intersection.Ts0[0] == end_ts0[i] && intersection.Ts1[0] == end_ts1[i]) {
                    b = false;
                    break;
                }
                if (intersection.Ts0[1] == end_ts0[i] && intersection.Ts1[1] == end_ts1[i]) {
                    b = false;
                    break;
                }
            }
        }
        if (b) {
            intersections[intersection_count] = CurveCurveIntersection(end_ts0[i], end_ts1[i], true);
            ++intersection_count;
        }
    }
}