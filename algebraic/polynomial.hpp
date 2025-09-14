#pragma once

#include "complex.hpp"

#include <functional>

template<typename T = Complex<>>
class Polynomial {
private:
    std::string DebugStr = "0";

    void _UpdateDebugStr() {
        DebugStr = ToString(*this);
    }

public:
    struct Term {
        uint32_t Exp { 0 };
        T Cof { 0 };
    };

private:
    std::vector<Term> Terms;

    void normalize() {
        Terms.erase(
            std::remove_if(Terms.begin(), Terms.end(),
                        [](const Term& t) { return t.Cof.IsZero(); }),
            Terms.end()
        );
    }

    // May invalidate other coefficient references
    // You must normalize after using those, or at least ensure the value is not left at zero
    // Zero coefficients are not allowed
    T& GetCof(uint32_t Exp) {
        auto it = Terms.begin();
        while (it != Terms.end()) {
            if (it->Exp == Exp) {
                return it->Cof;
            }
            if (it->Exp > Exp) {
                break;
            }
            ++it;
        }
        return Terms.emplace(it, Term { Exp })->Cof;
    }

public:
    Polynomial() = default;
    Polynomial(T Cof, uint32_t Exp) : Terms{ { Exp, Cof } } {
        normalize();

        _UpdateDebugStr();
    }


    T GetCof(uint32_t Exp) const {
        for (Term const& t : Terms) {
            if (t.Exp == Exp) {
                return t.Cof;
            }
        }
        return T { 0 };
    }
    Term GetLeadingTerm() const {
        return Terms.empty() ? Term { } : Terms.back();
    }
    uint32_t Degree() const {
        return GetLeadingTerm().Exp;
    }
    bool IsZero() const {
        return Terms.empty();
    }

    // Set this value to be the remainder of ((*this) / Divisor)
    // Store the quotient in OutQuotient
    void ApplyRemainder(const Polynomial Divisor, Polynomial& OutQuotient) {
        if (Divisor.IsZero()) throw std::runtime_error("Polynomial division has zero quotient");

        const Term DivisorLeading = Divisor.GetLeadingTerm();

        OutQuotient = Polynomial { };
        
        while (!IsZero()) {
            const Term LeadingTerm = GetLeadingTerm();

            if (DivisorLeading.Exp > LeadingTerm.Exp) {
                break;
            }

            const Polynomial Factor {
                LeadingTerm.Cof / DivisorLeading.Cof,
                LeadingTerm.Exp - DivisorLeading.Exp
            };

            OutQuotient += Factor;

            const Polynomial Tmp = Divisor * Factor;

            (*this) -= Tmp;
        }
    }

    // Set this value to be the Nth derivative of itself
    void ApplyDerivative(size_t N) {
        if (N == 0) return;

        Terms.erase(
            std::remove_if(Terms.begin(), Terms.end(),
                        [N](const Term& t) { return t.Exp < N; }),
            Terms.end()
        );

        for (auto& t : Terms) {
            T CofNew = t.Cof;
            for (uint32_t k = 0; k < N; ++k) {
                CofNew *= T(t.Exp - k);
            }
            t.Cof = CofNew;
            t.Exp -= N;
        }

        _UpdateDebugStr();
    }

    

    T Evaluate(const T& Value) const {
        T Res;

        for (const Term& t : Terms) {
            Res += t.Cof * T::Power(Value, t.Exp);
        }
        
        return Res;
    }

    static std::vector<Polynomial> MakeSturmSequence(Polynomial Val) {
        std::vector<Polynomial> Res;

        Res.push_back(Val);

        Val.ApplyDerivative(1);

        Res.push_back(Val);

        while (!Res.back().IsZero()) {
            Res.push_back(
                -(*(Res.rbegin() + 1) % Res.back())
            );
        }

        Res.pop_back();

        return Res;
    }

    // Returns the number of sign changes, or -1 if the value is a root
    static int32_t CountSignChanges(const std::vector<Polynomial>& Sturm, const decltype(T::Real)& Value, bool& OutIsRoot) {
        OutIsRoot = false;
        int32_t Res = 0;

        int32_t PriorSign = 0;
        for (size_t i = 0; i < Sturm.size(); ++i) {
            T Val = Sturm[i].Evaluate(Value);

            if (!Val.IsReal()) {
                throw std::runtime_error("Imaginary detected!");
            }

            if (Val.IsZero()) {
                if (i == 0) {
                    OutIsRoot = true;
                }
            } else {
                int32_t NewSign = Val.Real < 0 ? -1 : 1;
                if (PriorSign != 0 && PriorSign != NewSign) {
                    ++Res;
                }
                PriorSign = NewSign;
            }
        }
        
        return Res;
    }

    // Should be inclusive on lower, exclusive on upper
    // Might not work perfectly if lower or upper is a root
    static int32_t MinNumRootsEnclosed(
        const std::vector<Polynomial>& Sturm,
        const decltype(T::Real)& Lower,
        const decltype(T::Real)& Upper)
    {
        if (Lower == Upper) throw std::runtime_error("Region of size 0");

        bool LowerIsRoot, Unused;
        int32_t LowerSignChange = CountSignChanges(Sturm, Lower, LowerIsRoot);
        int32_t UpperSignChange = CountSignChanges(Sturm, Upper, Unused);

        return (LowerIsRoot ? 1 : 0) + abs(abs(LowerSignChange) - abs(UpperSignChange));
    }

    static std::vector<decltype(T::Real)> EvaluateRootsInRange(
        const std::vector<Polynomial>& Sturm,
        const decltype(T::Real)& Lower,
        const decltype(T::Real)& Upper,
        const decltype(T::Real)& MaxError
    ) {
        std::vector<decltype(T::Real)> Roots;

        std::function<void(decltype(T::Real), decltype(T::Real))> Bisect =
            [&](decltype(T::Real) A, decltype(T::Real) B) {
                int32_t NumRoots = MinNumRootsEnclosed(Sturm, A, B);
                if (NumRoots == 0) return;
                if (NumRoots == 1) {
                    if ((B - A) <= MaxError)
                        Roots.push_back((A + B) / 2);
                    else
                        Bisect(A, (A + B) / 2), Bisect((A + B) / 2, B);
                    return;
                }
                decltype(T::Real) Mid = (A + B) / 2;
                Bisect(A, Mid);
                Bisect(Mid, B);
            };

        Bisect(Lower, Upper);
        std::sort(Roots.begin(), Roots.end());
        return Roots;
    }

    // Abs value of all roots should be <= this
    static decltype(T::Real) CauchyBounds(const Polynomial& Value) {
        decltype(T::Real) Largest = 0;

        decltype(T::Real) Tmp;
        for (const Term& t : Value.Terms) {
            if (!t.Cof.IsReal()) {
                throw std::runtime_error("Non real coefficient detected");
            }

            Tmp = t.Cof.Real;

            Tmp.ApplyAbs();

            if (Tmp > Largest) {
                Largest = std::move(Tmp);
            }
        }

        Tmp = 1;

        return Largest / Value.GetLeadingTerm().Cof.Real + Tmp;
    }


    // Serde
    static std::string ToString(const Polynomial& Val, int64_t MaxDigits = 3) {
        if (Val.IsZero()) {
            return "0";
        }

        std::string Res;
        bool InvertNextSign = false;
        for (size_t i = Val.Terms.size() - 1;; --i) {
            const Term& t = Val.Terms[i];
            const Term n = (i == 0) ? Term() : Val.Terms[i - 1];
            if (!t.Cof.IsReal()) {
                Res += "(";
            }

            if (t.Cof != 1 || t.Exp == 0) {
                Res += T::ToString(InvertNextSign ? -t.Cof : t.Cof);
            }

            if (!t.Cof.IsReal()) {
                Res += ")";
            }

            if (t.Exp > 0) {
                Res += "x";
            }

            if (t.Exp > 1) {
                Res += "^" + std::to_string(t.Exp);
            }
            
            if (i == 0) break;

            if (n.Cof.IsReal() && n.Cof.Real < 0) {
                InvertNextSign = true;
                Res += " - ";
            } else {
                InvertNextSign = false;
                Res += " + ";
            }
        }

        return Res;
    }


    Polynomial operator-() const {
        Polynomial Res = *this;
        for (Term& t : Res.Terms) {
            t.Cof.ApplyNegate();
        }
        Res._UpdateDebugStr();
        return Res;
    }
    Polynomial& operator-=(const Polynomial Other) {
        for (Term const& o : Other.Terms) {
            GetCof(o.Exp) -= o.Cof;
        }
        normalize();
        _UpdateDebugStr();
        return *this;
    }
    Polynomial operator-(const Polynomial Other) const {
        Polynomial Res = *this;
        Res -= Other;
        return Res;
    }
    Polynomial& operator+=(const Polynomial Other) {
        for (Term const& o : Other.Terms) {
            GetCof(o.Exp) += o.Cof;
        }
        normalize();
        _UpdateDebugStr();
        return *this;
    }
    Polynomial operator+(Polynomial Other) const {
        Other += *this;
        return Other;
    }
    Polynomial& operator*=(const Polynomial& Other) {
        Polynomial Result;
        for (Term const& t : Terms) {
            for (Term const& o : Other.Terms) {
                Result += Polynomial(t.Cof * o.Cof, t.Exp + o.Exp);
            }
        }
        *this = std::move(Result);
        _UpdateDebugStr();
        return *this;
    }
    Polynomial operator*(Polynomial Other) const {
        Other *= *this;
        return Other;
    }
    Polynomial& operator/=(const Polynomial Other) {
        Polynomial _ = *this;
        _.ApplyRemainder(Other, *this);
        _UpdateDebugStr();
        return *this;
    }
    Polynomial operator/(const Polynomial Other) const {
        Polynomial _ = *this;
        Polynomial Quotient;
        _.ApplyRemainder(Other, Quotient);
        return Quotient;
    }
    Polynomial& operator%=(const Polynomial Other) {
        Polynomial _;
        ApplyRemainder(Other, _);
        _UpdateDebugStr();
        return *this;
    }
    Polynomial operator%(const Polynomial& Other) const {
        Polynomial Remainder = *this;
        Polynomial _;
        Remainder.ApplyRemainder(Other, _);
        return Remainder;
    }
};