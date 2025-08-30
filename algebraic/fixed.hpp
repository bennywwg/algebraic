#pragma once

#include "bignum.hpp"

#include <map>

template<typename T = BigInt<>>
class Rational {
    // Val = A / B
    T A { 0 };
    T B { 1 };

    void normalize() {
        T G = T::GCD(A, B);
        A = A / G;
        B = B / G;
        if (B.Sign() < 0) {
            A = -A;
            B = -B;
        }
    }

public:
    Rational() = default;
    Rational(int64_t Val) : A { Val }, B { 1 } { }
    Rational(const T& Val) : A(Val) { }
    Rational(const T& Val, const T& Denom) : A(Val), B(Denom) {
        if (Denom.IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");
        normalize();
    }

    bool IsZero() const { return A.IsZero(); }

    // Production functions
    static Rational Reciprocal(Rational Val) {
        if (Val.IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");

        std::swap(Val.A, Val.B);

        // We could call normalize, but it would redundantly calculate the GCD, so just fix the signs
        if (Val.B.Sign() < 0) {
            Val.A = -Val.A;
            Val.B = -Val.B;
        }

        return Val;
    }

    static Rational Pow(Rational LHS, size_t RHS) {
        LHS.A = T::Pow(LHS.A, RHS);
        LHS.B = T::Pow(LHS.B, RHS);

        LHS.normalize();
        
        return LHS;
    }

    // Serde
    static Rational FromString(const std::string& Value) {
        bool Neg = Value[0] == '-';
        std::string Str = Neg ? Value.substr(1) : Value;
        auto Pos = Str.find('.');
        Rational R;
        if (Pos == std::string::npos) {
            R.A = T::FromString(Str);
            R.B = T(1);
        } else {
            std::string IntPart = Str.substr(0, Pos);
            std::string FracPart = Str.substr(Pos + 1);
            T Scale = T(1);
            for (size_t I = 0; I < FracPart.size(); ++I) Scale = Scale * T(10);
            R.A = T::FromString(IntPart + FracPart);
            R.B = Scale;
        }
        if (Neg) R.A = -R.A;
        R.normalize();
        return R;
    }

    static std::string ToString(Rational Val, int64_t MaxDigits = 20) {
        std::string Res = Val.A.Sign() < 0 ? "-" : "";
        if (Val.A.Sign() < 0) Val.A = -Val.A;
        T Quot = Val.A / Val.B;
        T Rem = Val.A % Val.B;
        Res += T::ToString(Quot);

        if (Rem == T(0)) return Res;

        Res += '.';
        std::string Fraction;
        std::map<T, size_t> RemPos;

        while (Rem != T(0) && (MaxDigits-- > 0)) {
            if (RemPos.count(Rem)) {
                size_t RepeatStart = RemPos[Rem];
                Fraction.insert(RepeatStart, "(");
                Fraction += ")";
                break;
            }
            RemPos[Rem] = Fraction.size();
            Rem = Rem * T(10);
            T Digit = Rem / Val.B;
            Fraction += char('0' + Digit[0]);
            Rem = Rem % Val.B;
        }

        Res += Fraction;
        return Res;
    }

    void ComputeReciprocal() {
        if (IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");

        std::swap(A, B);

        // We could call normalize, but it would redundantly calculate the GCD, so just fix the signs
        if (B.Sign() < 0) {
            A = -A;
            B = -B;
        }
    }

    // Operators
    Rational operator+() const { return *this; }
    Rational operator-() const { Rational r = *this; r.A = -r.A; return r; }
    Rational operator+(const Rational& o) const { Rational r; r.A = A * o.B + o.A * B; r.B = B * o.B; r.normalize(); return r; }
    Rational operator-(const Rational& o) const { Rational r; r.A = A * o.B - o.A * B; r.B = B * o.B; r.normalize(); return r; }
    Rational operator*(const Rational& o) const { Rational r; r.A = A * o.A; r.B = B * o.B; r.normalize(); return r; }
    Rational operator/(const Rational& o) const { return (*this) * Reciprocal(o); }
    bool operator==(const Rational& o) const { return A == o.A && B == o.B; }
    bool operator!=(const Rational& o) const { return !(*this == o); }
    bool operator<(const Rational& o) const { return A * o.B < o.A * B; }
    bool operator>(const Rational& o) const { return o < *this; }
    bool operator<=(const Rational& o) const { return !(*this > o); }
    bool operator>=(const Rational& o) const { return !(*this < o); }
};