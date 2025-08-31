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
    Rational(int32_t Val) : A { Val } { }
    Rational(const T& Val) : A { Val } {}
    Rational(const T& Val, const T& Denom) : A { Val }, B { Denom } {
        if (Denom.IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");
        normalize();
    }

    bool IsZero() const {
        return A.IsZero();
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

    static std::string ToString(Rational Val, int64_t MaxDigits = 3) {
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

    const std::wstring ToString() const {
        std::string Res = ToString(*this);
        return std::wstring(Res.begin(), Res.end());
    }

    void ApplyAbs() {
        A.ApplyAbs();
    }
    void ApplyNegate() {
        A.ApplyNegate();
    }
    void ApplyReciprocal() {
        if (IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");

        std::swap(A, B);

        // We could call normalize, but it would redundantly calculate the GCD, so just fix the signs
        if (B.Sign() < 0) {
            A = -A;
            B = -B;
        }
    }

    // Operators
    Rational operator+() const {
        return *this;
    }
    Rational operator-() const {
        Rational r = *this;
        r.A.ApplyNegate();
        return r;
    }
    Rational& operator+=(const Rational Other) {
        A = A * Other.B + Other.A * B;
        B = B * Other.B;
        return *this;
    }
    Rational operator+(Rational Other) const {
        Other += *this;
        return Other;
    }
    Rational& operator-=(const Rational Other) {
        A = A * Other.B - Other.A * B;
        B = B * Other.B;
        normalize();
        return *this;
    }
    Rational operator-(const Rational& Other) const {
        Rational Res = *this;
        Res -= Other;
        return Res;
    }
    Rational& operator*=(const Rational Other) {
        A *= Other.A;
        B *= Other.B;
        normalize();
        return *this;
    }
    Rational operator*(Rational Other) const {
        Other *= *this;
        return Other;
    }
    Rational& operator/=(Rational Other) {
        Other.ApplyReciprocal();
        *this *= Other;
        return *this;
    }
    Rational operator/(const Rational& Other) const {
        Rational Res = *this;
        Res /= Other;
        return Res;
    }
    bool operator==(const Rational& Other) const { return A == Other.A && B == Other.B; }
    bool operator!=(const Rational& Other) const { return !(*this == Other); }
    bool operator<(const Rational& Other) const { return A * Other.B < Other.A * B; }
    bool operator>(const Rational& Other) const { return Other < *this; }
    bool operator<=(const Rational& Other) const { return !(*this > Other); }
    bool operator>=(const Rational& Other) const { return !(*this < Other); }
};