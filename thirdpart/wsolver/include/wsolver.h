#ifndef _WSOLVER_
#define _WSOLVER_

#if defined(_WINDOWS)
    #if defined(WSOLVER_EXPORTS)
        #define WSOLVER_API __declspec(dllexport)
    #elif defined(WSOLVER_STATIC)
        #define WSOLVER_API
    #else
        #define WSOLVER_API __declspec(dllimport)
    #endif
#else
    #define WSOLVER_API
#endif

#include <stdint.h>

#pragma warning(push)
#pragma warning(disable: 4200)

inline double prev_double(double d) {
    if (d == 0) {
        return d;
    }
    int64_t* p = (int64_t*)&d;
    if ((*p & 0x7ff0000000000000) == 0x7ff0000000000000) {
        return d;
    }
    if (*p > 0) {
        --(*p);
    }
    else if (*p < 0) {
        ++(*p);
    }
    return d;
}

inline double next_double(double d) {
    if (d == 0) {
        return d;
    }
    int64_t* p = (int64_t*)&d;
    if ((*p & 0x7ff0000000000000) == 0x7ff0000000000000) {
        return d;
    }
    if (*p > 0) {
        ++(*p);
    }
    else {
        --(*p);
    }
    return d;
}

typedef double WSReal;
#define WS_REAL_EPSILON 1E-12
#define prev_real(d) prev_double(d)
#define next_real(d) next_double(d)

class WSCache;
class WSMatrix;
class WSInterval;
class WSVector;
class WSSliceVector;
class WSMatrix;
class WSSliceMatrix;
class WSIntervalVector;
class WSIntervalMatrix;
class WSEquationVariable;
class WSEquationBasis;
class WSTerm;
class WSTermCalculator;
class WSEquationSystem;
class WSEquationsCache;
class WSEquation;

#define WS_PI 3.14159265358979323846

class WSOLVER_API WSInterval {
public:
    WSReal Min;
    WSReal Max;
public:
    WSInterval();
    WSInterval(const WSReal& d);
    WSInterval(const WSInterval& d);
    WSInterval(const WSReal& min, const WSReal& max);
    void Merge(const WSReal& d);
    void Merge(const WSInterval& d);
    WSReal Middle() const;
    WSReal Length() const;
    bool IsIntersected(const WSReal& d, const WSReal& epsilon) const;
    bool IsIntersected(const WSInterval& d, const WSReal& epsilon) const;
};

WSOLVER_API WSInterval operator+(const WSInterval& d1, const WSInterval& d2);
WSOLVER_API WSInterval operator-(const WSInterval& d);
WSOLVER_API WSInterval operator-(const WSInterval& d1, const WSInterval& d2);
WSOLVER_API WSInterval operator*(const WSInterval& d1, const WSInterval& d2);
WSOLVER_API WSInterval operator*(const WSInterval& d1, const WSReal& d2);
WSOLVER_API WSInterval operator*(const WSReal& d1, const WSInterval& d2);
WSOLVER_API WSInterval operator/(const WSInterval& d1, const WSInterval& d2);
WSOLVER_API WSInterval operator/(const WSInterval& d1, const WSReal& d2);
WSOLVER_API WSInterval merge(const WSInterval& d1, const WSInterval& d2);
WSOLVER_API WSInterval cos(const WSInterval& x);
WSOLVER_API WSInterval sin(const WSInterval& x);
WSOLVER_API WSInterval tan(const WSInterval& x);
WSOLVER_API WSInterval acos(const WSInterval& x);
WSOLVER_API WSInterval asin(const WSInterval& x);
WSOLVER_API WSInterval atan(const WSInterval& x);
WSOLVER_API WSInterval abs(const WSInterval& x);
WSOLVER_API WSInterval limited_pow_1(const WSInterval& x, int y);
WSOLVER_API WSInterval limited_pow_2(const WSInterval& x, int y);
WSOLVER_API WSInterval pow(const WSInterval& x, int y);
WSOLVER_API WSInterval pow(const WSInterval& x, WSReal y);
WSOLVER_API WSInterval pow(const WSInterval& x, const WSInterval& y);
WSOLVER_API WSInterval log(const WSInterval& x);

class WSOLVER_API WSPolynomialTerm {
public:
    WSInterval Coef;
    int Powers[0];
};

class WSOLVER_API WSPolynomial {
public:
    WSPolynomial(int variable_count, int term_capacity);
    virtual ~WSPolynomial();
    int GetVariableCount();
    int GetTermCount();
    WSPolynomialTerm* GetTerm(int index);
    void ExchangeTerm(int i, int j);
    void SetCapacity(int capacity, bool reduce);
    int NewTerm();
    void Clear();
    void Reset(int variable_count);
    void CopyFrom(WSPolynomial* src);
    void CalculateInfo(const int* cls_list, int& leader_variable_index, int& degree, int& total_degree);
    void GetLeaderTerms(int variable_index, int& max_degree, int* term_indices, int& count);
    bool PRem(WSPolynomial* polynomial, int variable_index, int max_term_count, WSPolynomial* result, 
        WSPolynomial* temp_polynomial, int& max_degree1, int* term_indices1, int& count1, int max_degree2, int* term_indices2, int count2);
    bool PRem(WSPolynomial* polynomial, int variable_index, int max_term_count, WSPolynomial* result, 
        WSPolynomial* temp_polynomial, int* temp_indices1, int* temp_indices2);
    void Normalize();
    void Remove(int i);
    void Add(WSPolynomialTerm* term);
    void AddLast();
public:
    char* Print();
private:
    void* m_terms;
    int m_variable_count;
    int m_term_count;
    int m_term_capacity;
};

struct WSPolynomialInfo {
    WSPolynomial* Polynomial;
    //int LeaderIndex;
    //int LeaderDegree;
    //double LeaderComplexity;
    int Type;   //0: intermediate  1: basic  2: terminated
};

class WSPolynomialSystemElimination {
public:
    WSPolynomialSystemElimination(int polynomial_variable_count, int polynomial_capacity, int max_term_count,
        int* leader_variables, int leader_variable_count);
    ~WSPolynomialSystemElimination();
    void AddPolynomial(WSPolynomial* polynomial);
    int Execute(WSPolynomial**& algebra_polynomials, int& algebra_polynomial_count);
    //int Execute(int start_polynomial_index, int recommended_leader_index);
    int Execute();
    int Execute(bool* sanctioned_leaders, int& start_polynomial_index);
    int Execute(bool* sanctioned_leaders, int leader_index, int& start_polynomial_index);
private:
    void CalculateComplexity(WSPolynomial* polynomial, int variable_index, int& degree, int& item_count);
    void CalculateComplexity(int start_polynomial_index, int polynomial_count, int variable_index,
        int& polynomial_index, int& degree, double& complexity);
    //void CalculateInfo(WSPolynomialInfo* polynomial_info, bool* sanctioned_leaders);
private:
    int m_max_term_count;
    int m_polynomial_variable_count;
    WSPolynomialInfo* m_polynomial_infos;
    int m_polynomial_capacity;
    int m_polynomial_count;
    int* m_leader_variables;
    int m_leader_variable_count;
    WSPolynomial* m_temp_polynomial1;
    WSPolynomial* m_temp_polynomial2;
    WSPolynomial* m_temp_polynomial3;
    int* m_temp_indices1;
    int* m_temp_indices2;
};

enum class WSIterateResult {
    NoRoot,
    ClearRoot,
    TerminateEarly,
    Fuzzy
};

class WSOLVER_API WSEquationBasis {
public:
    WSEquationBasis();
    virtual ~WSEquationBasis() {}
public:
    virtual WSTermCalculator* GetTermCalculator() const = 0;
    virtual int GetVariableCount() const = 0;
    virtual int GetVariableIndex(int index) const = 0;
    virtual void SetVariableIndex(int index, int variable_index) = 0;
    virtual WSEquationBasis* Clone() const = 0;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const = 0;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const = 0;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const = 0;
    virtual bool Equals(WSEquationBasis* basis) const = 0;
};

class WSOLVER_API WSTerm {
public:
    WSReal GetCoef() const;
    void SetCoef(const WSReal& coef);
    int GetIndex() const;
    WSEquation* GetEquation() const;
    WSTermCalculator* GetCalculator() const;
    void SetCalculator(WSTermCalculator* calculator);
    int GetValueCacheIndex() const;
    int GetPartialDerivativeCacheIndex() const;
    int GetLinearCacheIndex() const;
private:
    WSEquation* m_equation;
    WSTermCalculator* m_calculator;
    WSReal m_coef;
    int m_index;
    int m_value_cache_index;
    int m_partial_derivative_cache_index;
    int m_linear_cache_index;
    friend class WSEquation;
    friend class WSEquationSystem;
    friend class WSTermCalculator;
};

class WSOLVER_API WSTermCalculator {
public:
    virtual ~WSTermCalculator() {}
    virtual int GetPriority() const = 0;
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable) = 0;
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index) = 0;
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b) = 0;
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable) = 0;
};

class WSOLVER_API WSEquation {
public:
    WSEquation(WSEquationSystem* equation_system, int term_capacity, 
        bool is_check_enable, bool is_iterate_enable, bool is_linear_enable);
    virtual ~WSEquation();
    WSEquationSystem* GetEquations() const;
    int GetIndex() const;
    int GetTermCapacity() const;
    void ClearTerms();
    WSEquation* AddTerm(int basis_index, const WSReal& term_coef);
    int GetTermCount() const;
    WSTerm* GetTerm(int index);
    bool IsCheckEnable() const;
    bool IsIterateEnable() const;
    int GetIterateEnableIndex() const;
    bool IsLinearEnable() const;
    int GetLinearEnableIndex() const;
    void SetLinearEnable(bool enable);
    WSInterval CalculateValue(WSEquationsCache* cache);
    WSInterval CalculatePartialDerivative(WSEquationsCache* cache, int variable_index);
    void CalculateLinear(WSEquationsCache* cache, WSVector* a, WSInterval* b);
    WSInterval CalculateValue(const WSIntervalVector* variable);
    WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index);
private:
    WSEquationSystem* m_equations;
    int m_index;
    WSTerm* m_terms;
    int m_term_capacity;
    int m_term_count;
    bool m_is_check_enable;
    bool m_is_iterate_enable;
    bool m_is_linear_enable;
    int m_iterate_enable_index;
    int m_linear_enable_index;
    friend class WSEquationSystem;
    friend class WSTermCalculator;
};

class WSOLVER_API WSCache {
public:
    WSCache();
    WSCache(int size);
    virtual ~WSCache();
    int GetSize();
    void Resize(int size);
    void* Data(int offset);
private:
    int m_size;
    void* m_data;
};

class WSOLVER_API WSEquationSystem {
public:
    WSEquationSystem();
    virtual ~WSEquationSystem();
    void BuildRuntime();
    void RebuildLinearRuntime();
    void ClearAlgebraRuntime();
    void RebuildAlgebraRuntime(const WSIntervalVector* variable_domain);
    bool GetDeterminedNoRoot() const;
    int GetIterateEnableEquationCount() const;
    int GetLinearEnableEquationCount() const;
    int GetTermValueCacheCount() const;
    int GetTermPartialDerivativeCacheCount() const;
    int GetTermLinearCacheCount() const;
    int GetVariableIterateEquationCount(int variable_index) const;
    int GetVariableIterateEquationIndex(int variable_index, int index) const;
    int GetVariableAlgebraPolynomialCount(int variable_index) const;
    int GetVariableAlgebraPolynomialIndex(int variable_index, int index) const;
    int GetVariableTermCount(int variable_index) const;
    WSTerm* GetVariableTerm(int variable_index, int index) const;
    int GetAlgebraPolynomialCount() const;
    WSPolynomial* GetAlgebraPolynomial(int index) const;
    WSInterval CalculateAlgebraPolynomialValue(WSPolynomial* polynomial, const WSIntervalVector* variable, WSReal& epsilon_coef);
    WSInterval CalculateAlgebraPolynomialPartialDerivative(WSPolynomial* polynomial, const WSIntervalVector* variable, int variable_index);
public:
    virtual int GetVariableCount() const = 0;
    virtual WSReal GetVariableEpsilon(const WSIntervalVector* variable, int index) const = 0;
    virtual int GetBasisCount() const = 0;
    virtual WSEquationBasis* GetBasis(int index) const = 0;
    virtual int GetEquationCount() const = 0;
    virtual WSEquation* GetEquation(int index) const = 0;
    virtual WSReal GetEquationCheckEpsilon(WSEquationsCache* cache, int index) const = 0;
    virtual WSReal GetEquationIterateEpsilon(WSEquationsCache* cache, int index) const = 0;
    virtual void BuildOriginalAlgebraEquations(const WSIntervalVector* variable_domain, 
        WSPolynomial**& polynomials, int& polynomial_count, int*& leader_variables, int& leader_variable_count, int& max_term_count) const = 0;
    virtual void BuildAlgebraWhenEliminationFail(const WSIntervalVector* variable_domain,
        WSPolynomial**& polynomials, int& polynomial_count) const {
        polynomials = nullptr;
        polynomial_count = 0;
    };
private:
    void RebuildIterateRuntime();
    void AddVariableIterateEquationIndex(int variable_index, int equation_index);
    void AddVariableAlgebraPolynomialIndex(int variable_index, int algebra_polynomial_index);
    void AddVariableTerm(int variable_index, WSTerm* term);
private:
    WSCache m_runtime_cache;
    int m_iterate_enable_equation_count;
    int m_linear_enable_equation_count;
    int m_term_value_cache_count;
    int m_term_partial_derivative_cache_count;
    int m_term_linear_cache_count;
    struct VariableRelations {
        int IterateEquationCount;
        int* IterateEquationIndices;
        int TermCount;
        WSTerm** Terms;
        int AlgebraPolynomialCount;
        int* AlgebraPolynomialIndices;
    };
    VariableRelations* m_variable_relations;
    int m_algebra_polynomial_count;
    WSPolynomial** m_algebra_polynomials;
    bool m_determined_no_root;
};

class WSOLVER_API WSVector {
public:
    WSVector();
    WSVector(int dimension);
    WSVector(int dimension, const WSReal* data);
    virtual ~WSVector();
    int GetDimension() const;
    const WSReal& Get(int index) const;
    void Set(int index, const WSReal& value);
public:
    void LoadZeros();
    void CopyFrom(const WSVector* src);
    void Add(const WSVector* vt);
    void Add(const WSVector* vt, WSVector* result) const;
    void Sub(const WSVector* vt);
    void Mul(const WSReal& d);
    void Mul(const WSReal& d, WSVector* result) const;
    void MulAdd(const WSReal& d, const WSVector* vt);
    WSReal Dot(const WSVector* vt);
protected:
    int m_dimension;
    WSReal* m_data;
};

class WSOLVER_API WSSliceVector : public WSVector {
public:
    WSSliceVector(WSCache* cache, int offset, int dimension);
    WSSliceVector(void* data, int dimension);
    virtual ~WSSliceVector();
};

class WSOLVER_API WSMatrix {
public:
    WSMatrix();
    WSMatrix(int row_count, int col_count);
    WSMatrix(int row_count, int col_count, const WSReal* data);
    virtual ~WSMatrix();
    int GetRowCount() const;
    int GetColCount() const;
    const WSReal& Get(int row, int col) const;
    void Set(int row, int col, const WSReal& value);
public:
    void LoadZeros();
    void LoadIdentity();
    void SetRowZeros(int row);
    void SetColZeros(int col);
    WSSliceVector GetRowVector(int row) const;
    void CopyFrom(const WSMatrix* src);
    void ExchangeRow(int row0, int row1, int start_col, int end_col);
    void ExchangeCol(int col0, int col1, int start_row, int end_row);
    void Transpose() const;
    void Transpose(WSMatrix* result) const;
    void Mul(const WSVector* vt, WSVector* result) const;
    void Mul(const WSIntervalVector* vt, WSIntervalVector* result) const;
    void Mul(const WSMatrix* matrix, WSMatrix* result) const;
    void LeftMul(const WSVector* vt, WSVector* result) const;
    void LeftMul(const WSIntervalVector* vt, WSIntervalVector* result) const;
    WSReal MulRow(int row, const WSVector* vt) const;
    WSReal MulCol(int col, const WSVector* vt) const;
    void QR(WSMatrix* q, WSMatrix* r) const;
protected:
    int m_row_count;
    int m_col_count;
    WSReal* m_data;
};

class WSOLVER_API WSSliceMatrix : public WSMatrix {
public:
    WSSliceMatrix();
    WSSliceMatrix(WSCache* cache, int offset, int row_count, int col_count);
    virtual ~WSSliceMatrix();
};

class WSOLVER_API WSIntervalVector {
public:
    WSIntervalVector();
    WSIntervalVector(int dimension);
    virtual ~WSIntervalVector();
    int GetDimension() const;
    WSInterval* GetPointer(int index);
    const WSInterval& Get(int index) const;
    void Set(int index, const WSInterval& value);
    void CopyFrom(const WSIntervalVector* src);
    bool Merge(const WSIntervalVector& src, const WSReal& epsilon);
public:
    void LoadZeros();
    void Neg();
protected:
    int m_dimension;
    WSInterval* m_data;
};

class WSOLVER_API WSSliceIntervalVector : public WSIntervalVector {
public:
    WSSliceIntervalVector();
    WSSliceIntervalVector(WSCache* cache, int cache_offset, int dimension);
    virtual ~WSSliceIntervalVector();
};

class WSOLVER_API WSIntervalMatrix {
public:
    WSIntervalMatrix();
    WSIntervalMatrix(int row_count, int col_count);
    virtual ~WSIntervalMatrix();
    int GetRowCount() const;
    int GetColCount() const;
    WSInterval* GetPointer(int row, int col);
    const WSInterval& Get(int row, int col) const;
    void Set(int row, int col, const WSInterval& value);
    void CopyFrom(const WSIntervalMatrix* src);
protected:
    int m_row_count;
    int m_col_count;
    WSInterval* m_data;
};

class WSOLVER_API WSSliceIntervalMatrix : public WSIntervalMatrix {
public:
    WSSliceIntervalMatrix();
    WSSliceIntervalMatrix(WSCache* cache, int cache_offset, int row_count, int col_count);
    virtual ~WSSliceIntervalMatrix();
};

class WSOLVER_API WSEquationsCache {
public:
    WSEquationsCache();
    WSEquationsCache(WSEquationSystem* equations, WSIntervalVector* variable, WSIntervalVector* terms_value, 
        WSIntervalMatrix* terms_partial_derivative, WSMatrix* terms_linear_a, WSIntervalVector* terms_linear_b, 
        int* terms_dirty_flag);
    WSEquationSystem* GetEquations();
    WSInterval GetVariable(int variable_index);
    WSIntervalVector* GetVariable();
    void SetVariable(int variable_index, const WSInterval& d);
    void SetTermsDirty(int changed_variable_index);
    WSInterval GetTermsValue(WSTerm* term);
    WSInterval GetTermsPartialDerivative(WSTerm* term, int variable_index);
    WSSliceVector GetTermsLinearA(WSTerm* term);
    WSInterval GetTermsLinearB(WSTerm* term);
    void CopyFrom(const WSEquationsCache* src);
private:
    WSEquationSystem* m_equations;
    WSIntervalVector* m_variable;
    WSIntervalVector* m_terms_value;
    WSIntervalMatrix* m_terms_partial_derivative;
    WSMatrix* m_terms_linear_a;
    WSIntervalVector* m_terms_linear_b;
    int* m_terms_dirty_flag;
};

class WSOLVER_API WSIterator {
public:
    WSIterator(WSEquationSystem* equations, WSCache* cache, int cache_offset);
    static int GetCacheSize(WSEquationSystem* equations);
    WSIterateResult Execute(WSEquationsCache* cache, WSReal min_rate, int max_iterate_count, 
        int& iterate_count, int& terminate_early_tag);
    WSReal GetSplitPriorityByLinear(WSEquationsCache* cache, int variable_index);
protected:
    virtual bool CheckTerminateEarly(WSEquationSystem* equations, WSEquationsCache* cache, int& tag);
private:
    WSEquationSystem* m_equations;
    WSSliceIntervalMatrix m_current_partial_derivative;
    WSSliceIntervalMatrix m_current_algebra_polynomial_partial_derivative;
    WSIntervalVector* m_old_variable;
    WSEquationsCache m_min_cache;
    WSEquationsCache m_max_cache;
    WSSliceMatrix m_current_linear_a;
    WSSliceMatrix m_current_linear_q;
    WSSliceMatrix m_current_linear_r;
    WSSliceIntervalMatrix m_current_linear_r2;
    WSSliceIntervalVector m_current_linear_b;
    WSSliceIntervalVector m_current_linear_b2;
    WSSliceIntervalVector m_temp_coef_variable;
};

template<class HeapItem> void HeapMoveDown(HeapItem* heap, int heap_item_count, int i) {
    while (true) {
        int j = i * 2 + 1;
        if (j >= heap_item_count) {
            break;
        }
        int k = j + 1;
        if (k < heap_item_count) {
            if (heap[k].LessThan(heap[j])) {
                j = k;
            }
        }
        if (heap[j].LessThan(heap[i])) {
            HeapItem temp_item = heap[i];
            heap[i] = heap[j];
            heap[j] = temp_item;
            i = j;
        }
        else {
            break;
        }
    }
}

template<class HeapItem> void HeapMoveUp(HeapItem* heap, int heap_item_count, int i) {
    while (i > 0) {
        int j = (i - 1) / 2;
        if (heap[i].LessThan(heap[j])) {
            HeapItem temp_item = heap[i];
            heap[i] = heap[j];
            heap[j] = temp_item;
            i = j;
        }
        else {
            break;
        }
    }
}

template<class HeapItem> void HeapPop(HeapItem* heap, int& heap_item_count) {
    --heap_item_count;
    if (heap_item_count == 0) {
        return;
    }
    HeapItem temp_item = heap[0];
    heap[0] = heap[heap_item_count];
    heap[heap_item_count] = temp_item;
    HeapMoveDown(heap, heap_item_count, 0);
}

WSOLVER_API bool real_to_fraction(WSReal x, int& numerator, int& denominator, int max_denominator = 1000);

#pragma warning(pop)

#endif