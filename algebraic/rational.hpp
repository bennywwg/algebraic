#pragma once

#include "bignum.hpp"

#include <map>

template<typename F>
struct FloatTraits;

template<>
struct FloatTraits<float> {
    using bits_t = uint32_t;
    static constexpr bits_t exponent_mask = 0x7F800000u;
    static constexpr bits_t fraction_mask = 0x007FFFFFu;
    static constexpr uint32_t exponent_shift = 23;
    static constexpr uint32_t exponent_bias  = 127;
};

template<>
struct FloatTraits<double> {
    using bits_t = uint64_t;
    static constexpr bits_t exponent_mask = 0x7FF0000000000000ull;
    static constexpr bits_t fraction_mask = 0x000FFFFFFFFFFFFFull;
    static constexpr uint32_t exponent_shift = 52;
    static constexpr uint32_t exponent_bias  = 1023;
};

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
    Rational(int64_t Val) : A { Val } { }
    template<std::floating_point F>
    Rational(F Val) {
        static_assert(std::is_same_v<F, double> || std::is_same_v<F, float>);
        static_assert(std::numeric_limits<F>::is_iec559);
        
        using traits = FloatTraits<F>;
        using bits_t = traits::bits_t;

        const bits_t bits = std::bit_cast<bits_t>(Val);
        
        const bool sign = Val < 0;
        const bits_t fraction = bits & traits::fraction_mask;
        const bits_t exponent = (bits & traits::exponent_mask) >> traits::exponent_shift;

        const bool isDenormalized = (exponent == 0) && (fraction != 0);
        
        *this = Pow(2, static_cast<int64_t>(isDenormalized ? bits_t(1) : exponent) - static_cast<int64_t>(traits::exponent_bias));
        *this *= (Rational(isDenormalized ? 0 : 1u) + Rational(T(fraction), T::Power2(traits::exponent_shift)));
        A.ApplySign(sign);
    }
    Rational(const T& Val) : A { Val } {}
    Rational(const T& Val, const T& Denom) : A { Val }, B { Denom } {
        if (Denom.IsZero()) throw std::runtime_error("Attempting reciprocal of zero rational");
        normalize();
    }

    T Numerator() const {
        return A;
    }
    T Denominator() const {
        return B;
    }

    T Floor() const {
        return A / B - (A < 0 && A % B);
    }
    T Ceil() const {
        return A / B + (A > 0 && A % B);
    }
    T Round() const {
        return A >= 0 ? (A + (B >> 1)) / B : (A - (B >> 1)) / B;
    }

    bool IsZero() const {
        return A.IsZero();
    }

    static Rational Pow(Rational LHS, int64_t RHS) {
        LHS.A = T::Pow(LHS.A, abs(RHS));
        LHS.B = T::Pow(LHS.B, abs(RHS));

        if (RHS < 0) {
            LHS.ApplyReciprocal();
        }

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

    static std::string ToString(Rational Val, int64_t MaxDigits = 10) {
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