#include "bignum.hpp"

#include <map>

template<typename F = uint64_t, typename H = uint32_t>
requires WideEnough<F, H>
class Rational {
    using Num = BigInt<F, H>;

    // Val = A / B
    Num A;
    Num B { 1 };

    void normalize() {
        BigInt G = Num::GCD(A, B);
        A = A / G;
        B = B / G;
        if (B.Sign() < 0) {
            A = -A;
            B = -B;
        }
    }
public:
    static Rational FromString(const std::string& Value) {
        bool Neg = Value[0] == '-';
        std::string Str = Neg ? Value.substr(1) : Value;
        auto Pos = Str.find('.');
        Rational R;
        if (Pos == std::string::npos) {
            R.A = Num::FromString(Str);
            R.B = Num(1);
        } else {
            std::string IntPart = Str.substr(0, Pos);
            std::string FracPart = Str.substr(Pos + 1);
            Num Scale = Num(1);
            for (size_t I = 0; I < FracPart.size(); ++I) Scale = Scale * Num(10);
            R.A = Num::FromString(IntPart + FracPart);
            R.B = Scale;
        }
        if (Neg) R.A = -R.A;
        R.normalize();
        return R;
    }

    /*static std::string ToString(Rational Val, uint64_t Digits) {
        std::string Res = Val.A.Sign() < 0 ? "-" : "";
        Num Quot = Val.A / Val.B;
        Num Rem = Val.A % Val.B;
        Res += Num::ToString(Num::Abs(Quot));
        if (Digits > 0) {
            Res += '.';
            for (uint64_t i = 0; i < Digits; ++i) {
                Rem = Rem * Num(10);
                Num Digit = Rem / Val.B;
                Rem = Rem % Val.B;
                Res += char('0' + int64_t(Digit[0]));
            }
        }
        return Res;
    }*/

    static std::string ToString(Rational Val, int64_t MaxDigits = 20) {
        std::string Res = Val.A.Sign() < 0 ? "-" : "";
        if (Val.A.Sign() < 0) Val.A = -Val.A;
        Num Quot = Val.A / Val.B;
        Num Rem = Val.A % Val.B;
        Res += Num::ToString(Quot);

        if (Rem == Num(0)) return Res;

        Res += '.';
        std::string Fraction;
        std::map<Num, size_t> RemPos;

        while (Rem != Num(0) && (MaxDigits-- > 0)) {
            if (RemPos.count(Rem)) {
                size_t RepeatStart = RemPos[Rem];
                Fraction.insert(RepeatStart, "(");
                Fraction += ")";
                break;
            }
            RemPos[Rem] = Fraction.size();
            Rem = Rem * Num(10);
            Num Digit = Rem / Val.B;
            Fraction += char('0' + Digit[0]);
            Rem = Rem % Val.B;
        }

        Res += Fraction;
        return Res;
    }


    // Operators
    Rational operator+() const { return *this; }
    Rational operator-() const { Rational r = *this; r.A = -r.A; return r; }
    Rational operator+(const Rational& o) const { Rational r; r.A = A * o.B + o.A * B; r.B = B * o.B; r.normalize(); return r; }
    Rational operator-(const Rational& o) const { Rational r; r.A = A * o.B - o.A * B; r.B = B * o.B; r.normalize(); return r; }
    Rational operator*(const Rational& o) const { Rational r; r.A = A * o.A; r.B = B * o.B; r.normalize(); return r; }
    Rational operator/(const Rational& o) const { Rational r; r.A = A * o.B; r.B = B * o.A; r.normalize(); return r; }
    bool operator==(const Rational& o) const { return A == o.A && B == o.B; }
    bool operator!=(const Rational& o) const { return !(*this == o); }
};