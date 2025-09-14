#pragma once

#include "rational.hpp"

template<typename T = Rational<>>
class Complex {
public:
    T Real;
    T Imag;

    Complex() = default;
    Complex(int32_t Val) : Real { Val } {}
    Complex(T InReal) : Real(InReal) { }
    Complex(T InReal, T InImag) : Real(InReal), Imag(InImag) { }

    bool IsZero() const {
        return Real.IsZero() && Imag.IsZero();
    }
    bool IsReal() const {
        return Imag.IsZero();
    }

    static Complex MakeReal(T Val) {
        return { Val, { } };
    }
    static Complex MakeImag(T Val) {
        return { { }, Val };
    }
    

    void ApplyNegate() {
        Real.ApplyNegate();
        Imag.ApplyNegate();
    }

    // Serde
    static std::string ToString(Complex Val, int64_t MaxDigits = 10) {
        if (Val.Real.IsZero() && Val.Imag.IsZero()) {
            return "0";
        }
        if (Val.Real.IsZero()) {
            return T::ToString(Val.Imag) + "i";
        } else {
            return T::ToString(Val.Real);
        }
    }

    const std::wstring ToString() const {
        return L"ABC";
        std::string Res = ToString(*this);
        return std::wstring(Res.begin(), Res.end());
    }

    static Complex Power(Complex LHS, size_t RHS) {
        Complex Result;
        Result.Real = T { 1 };
        Result.Imag = T { 0 };

        while (RHS > 0) {
            if (RHS % 2 == 1) {
                T RealTemp = Result.Real * LHS.Real - Result.Imag * LHS.Imag;
                Result.Imag = Result.Real * LHS.Imag + Result.Imag * LHS.Real;
                Result.Real = RealTemp;
            }
            T RealTemp = LHS.Real * LHS.Real - LHS.Imag * LHS.Imag;
            LHS.Imag = T(2) * LHS.Real * LHS.Imag;
            LHS.Real = RealTemp;
            RHS /= 2;
        }

        return Result;
    }

    Complex operator+() const {
        return *this;
    }
    Complex operator-() const {
        Complex C = *this;
        C.Real = -C.Real;
        C.Imag = -C.Imag;
        return C;
    }
    Complex& operator+=(const Complex& Other) {
        Real += Other.Real;
        Imag += Other.Imag;
        return *this;
    }
    Complex operator+(Complex Other) const {
        Other += *this;
        return Other;
    }
    Complex& operator-=(const Complex& Other) {
        Real -= Other.Real;
        Imag -= Other.Imag;
        return *this;
    }
    Complex operator-(Complex Other) const {
        Other.Real = Real - Other.Real;
        Other.Imag = Imag - Other.Imag;
        return Other;
    }
    Complex& operator*=(const Complex& Other) {
        T TmpReal = Real;
        Real = Real * Other.Real - Imag * Other.Imag;
        Imag = TmpReal * Other.Imag + Imag * Other.Real;
        return *this;
    }
    Complex operator*(Complex Other) const {
        Other *= *this;
        return Other;
    }
    Complex& operator /=(const Complex& Other) {
        // TODO could be more efficient
        T Den = Other.Real * Other.Real + Other.Imag * Other.Imag;
        Den.ApplyReciprocal();
        Complex C;
        C.Real = (Real * Other.Real + Imag * Other.Imag) * Den;
        C.Imag = (Imag * Other.Real - Real * Other.Imag) * Den;
        *this = C;
        return *this;
    }
    Complex operator/(const Complex& Other) const {
        Complex Res = *this;
        Res /= Other;
        return Res;
    }
    bool operator==(const Complex& Other) const { return Real == Other.Real && Imag == Other.Imag; }
    bool operator!=(const Complex& Other) const { return !(*this == Other); }
};