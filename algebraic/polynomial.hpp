#pragma once

#include "complex.hpp"

template<typename T = Complex<>>
class Polynomial {
public:
    struct Term {
        uint32_t Exp { 0 };
        T Cof { 0 };
    };

private:
    std::vector<Term> Terms;

    void sort() {
        std::sort(Terms.begin(), Terms.end(), [](const Term& a, const Term& b) {
            return a.Exp < b.Exp;
        });
    }

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

    Polynomial() = default;
    Polynomial(T Cof, uint32_t Exp) : Terms { { Exp, Cof } } {
        normalize();
    }

    // Set this value to be the remainder of ((*this) / Divisor)
    // Store the quotient in OutQuotient
    void ComputeRemainder(const Polynomial& Divisor, Polynomial& OutQuotient) {
        if (Divisor.IsZero()) throw std::runtime_error("Polynomial division has zero quotient");

        const Term DivisorLeading = Divisor.GetLeadingTerm();

        OutQuotient = Polynomial { };
        
        while (!IsZero() && DivisorLeading.Exp <= GetLeadingTerm().Exp) {
            const Polynomial Factor {
                GetLeadingTerm().Cof / DivisorLeading.Cof,
                GetLeadingTerm().Exp - DivisorLeading.Exp
            };
            OutQuotient += Factor;
            (*this) -= Divisor * Factor;
        }
    }

    // Set this value to be the Nth derivative of itself
    void ComputeDerivative(size_t N) {
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
    }


    // Serde
    static std::string ToString(const Polynomial& Val, int64_t MaxDigits = 20) {
        if (Val.IsZero()) {
            return "0";
        }

        std::string Res;
        for (size_t i = Val.Terms.size() - 1;; --i) {
            const Term& t = Val.Terms[i];
            if (!t.Cof.IsReal()) {
                Res += "(";
            }

            Res += T::ToString(t.Cof);

            if (!t.Cof.IsReal()) {
                Res += ")";
            }

            if (t.Exp != 0) {
                Res += "x^" + std::to_string(t.Exp);
            }
            
            if (i == 0) break;

            Res += " ";
        }

        return Res;
    }


    Polynomial operator-() const {
        Polynomial Res = *this;
        for (Term& t : Res.Terms) {
            t.Cof = T::Negate(t.Cof);
        }
        return Res;
    }
    Polynomial& operator -=(const Polynomial& Other) {
        for (Term const& o : Other.Terms) {
            GetCof(o.Exp) -= o.Cof;
        }
        normalize();
        return *this;
    }
    Polynomial operator-(const Polynomial& Other) {
        Polynomial Res = *this;
        Res -= Other;
        return Res;
    }
    Polynomial& operator +=(const Polynomial& Other) {
        for (Term const& o : Other.Terms) {
            GetCof(o.Exp) += o.Cof;
        }
        normalize();
        return *this;
    }
    Polynomial operator+(Polynomial Other) const {
        Other += *this;
        return Other;
    }
    Polynomial& operator*=(const Polynomial& Other) {
        for (Term const& t : Terms) {
            for (Term const& o : Other.Terms) {
                (*this) += Polynomial(t.Cof * o.Cof, t.Exp + o.Exp);
            }
        }
        return *this;
    }
    Polynomial operator*(Polynomial Other) const {
        Other *= *this;
        return Other;
    }
    Polynomial& operator/=(const Polynomial& Other) {
        Polynomial _ = *this;
        _.ComputeRemainder(Other, *this);
        return *this;
    }
    Polynomial operator/(const Polynomial& Other) const {
        Polynomial _ = *this;
        Polynomial Quotient;
        _.ComputeRemainder(Other, Quotient);
        return Quotient;
    }
    Polynomial& operator%=(const Polynomial& Other) {
        Polynomial _;
        ComputeRemainder(Other, _);
        return *this;
    }
    Polynomial operator%(const Polynomial& Other) const {
        Polynomial Remainder = *this;
        Polynomial _;
        Remainder.ComputeRemainder(Other, _);
        return Remainder;
    }
};