#ifndef _WSTERM_
#define _WSTERM_

#include "wsolver.h"

class WSConstCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSConstCalculator Instance;
};

class WSConstBasis : public WSEquationBasis {
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

class WSPowerCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSPowerCalculator Instance;
};

class WSPowerBasis : public WSEquationBasis {
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

class WSMulCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSMulCalculator Instance;
};

class WSMulBasis : public WSEquationBasis {
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

class WSSinCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSSinCalculator Instance;
};

class WSSinBasis : public WSEquationBasis {
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

class WSCosCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSCosCalculator Instance;
};

class WSCosBasis : public WSEquationBasis {
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

class WSLnCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSLnCalculator Instance;
};

class WSLnBasis : public WSEquationBasis {
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

class WSAbsCalculator : public WSTermCalculator {
public:
    virtual int GetPriority() const { return 0; }
    virtual WSInterval CalculateValue(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable);
    virtual WSInterval CalculatePartialDerivative(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, int variable_index);
    virtual void CalculateLinear(WSEquationSystem* equations, WSTerm* term, const WSIntervalVector* variable, WSVector* a, WSInterval* b);
    virtual int GetSplitIndex(WSEquationSystem* equations, WSTerm* term, WSIntervalVector* variable);
public:
    static WSAbsCalculator Instance;
};

class WSAbsBasis : public WSEquationBasis {
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

class WSBernsteinCalculator : public WSTermCalculator {
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

class WSBernsteinBasis : public WSEquationBasis {
public:
    WSBernsteinBasis(int variable_index, int degree);
    WSBernsteinBasis(int variable_index, int degree, const WSReal* coefs);
    int GetDegree() const;
    void SetCoef(int index, const WSReal& coef);
    WSReal GetCoef(int index) const;
    const WSReal* GetCoefs() const;
    WSReal CalculateValue(const WSReal& variable) const;
    WSReal CalculateDerivative(const WSReal& variable) const;
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
