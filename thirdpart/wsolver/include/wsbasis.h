#ifndef _WSTERM_
#define _WSTERM_

#include "wsolver.h"

class WSOLVER_API WSConstCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSConstCalculator Instance;
};

class WSOLVER_API WSConstBasis : public WSEquationBasis {
public:
    WSConstBasis();
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
};

class WSOLVER_API WSPowerCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSPowerCalculator Instance;
};

class WSOLVER_API WSPowerBasis : public WSEquationBasis {
public:
    WSPowerBasis(int variable_index, int power);
    int GetPower() const;
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
    int m_power;
};

class WSOLVER_API WSMulCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSMulCalculator Instance;
};

class WSOLVER_API WSMulBasis : public WSEquationBasis {
public:
    WSMulBasis(int variable_index0, int variable_index1);
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_indices[2];
};

class WSOLVER_API WSSinCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSSinCalculator Instance;
};

class WSOLVER_API WSSinBasis : public WSEquationBasis {
public:
    WSSinBasis(int variable_index);
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
};

class WSOLVER_API WSCosCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSCosCalculator Instance;
};

class WSOLVER_API WSCosBasis : public WSEquationBasis {
public:
    WSCosBasis(int variable_index);
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
};

class WSOLVER_API WSLnCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSLnCalculator Instance;
};

class WSOLVER_API WSLnBasis : public WSEquationBasis {
public:
    WSLnBasis(int variable_index);
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
};

class WSOLVER_API WSAbsCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSAbsCalculator Instance;
};

class WSOLVER_API WSAbsBasis : public WSEquationBasis {
public:
    WSAbsBasis(int variable_index);
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
};

class WSOLVER_API WSBernsteinCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static void SubMinSection(int degree, WSReal* coefs, const WSReal& t);
    static void SubMaxSection(int degree, WSReal* coefs, const WSReal& t);
    static void SubSection(int degree, WSReal* coefs, const WSInterval& domain);
public:
    static WSBernsteinCalculator Instance;
};

class WSOLVER_API WSBernsteinBasis : public WSEquationBasis {
public:
    WSBernsteinBasis(int variable_index, int degree);
    WSBernsteinBasis(int variable_index, int degree, const WSReal* coefs);
    int GetDegree() const;
    void SetCoef(int index, const WSReal& coef);
    WSReal GetCoef(int index) const;
    const WSReal* GetCoefs() const;
    WSReal CalculateValue(const WSReal& variable) const;
    WSReal CalculateDerivative(const WSReal& variable) const;
    WSInterval CalculateValue(const WSInterval& variable) const;
    WSInterval CalculateDerivative(const WSInterval& variable) const;
    virtual WSTermCalculator* GetTermCalculator() const;
    virtual int GetVariableCount() const;
    virtual int GetVariableIndex(int index) const;
    virtual void SetVariableIndex(int index, int variable_index);
    virtual WSEquationBasis* Clone() const;
    virtual void AddPolynomialTerm(WSPolynomial* polynomial, const WSInterval& coef) const;
    virtual WSInterval CalculateValue(const WSIntervalVector* variable) const;
    virtual WSInterval CalculatePartialDerivative(const WSIntervalVector* variable, int variable_index) const;
    virtual bool Equals(WSEquationBasis* basis) const;
private:
    int m_variable_index;
    int m_degree;
    WSReal* m_coefs;
};

#endif
