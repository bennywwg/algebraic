#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <assert.h>

template<typename F, typename H>
concept WideEnough = std::unsigned_integral<F> && std::unsigned_integral<H> && (sizeof(F) == 2 * sizeof(H));

template<typename F = uint64_t, typename H = uint32_t>
requires WideEnough<F, H>
class BigInt {
    bool m_Sign = false;
    std::vector<H> m_Data;

    static constexpr H lsb(F x) { return static_cast<H>(x); }
    static constexpr H msb(F x) { return static_cast<H>(x >> (sizeof(H) * CHAR_BIT)); }

    void normalize() {
        while (!m_Data.empty()) {
            if (m_Data.back() != 0)
                break;
            m_Data.pop_back();
        }
        
        if (m_Data.empty()) {
            m_Sign = false;
        }
    }

public:
    BigInt() = default;
    explicit BigInt(H Val) {
        if (Val != 0) {
            m_Data = { Val };
        }
    }
    
    // Const functions
    size_t Size() const { return m_Data.size(); }
    int32_t Sign() const { return IsZero() ? 0 : (m_Sign ? -1 : 1); }
    H operator[](size_t Index) const { return Index < m_Data.size() ? m_Data[Index] : 0; }
    H& operator[](size_t Index) {
        if (Index >= m_Data.size()) {
            m_Data.resize(Index + 1, 0);
        }
        return m_Data[Index];
    }
    bool IsZero() const { return m_Data.empty(); }
    size_t CountBits() const {
        if (IsZero()) return 0;

        size_t LastIndex = Size() - 1;
        H TopWord = m_Data[LastIndex];

        size_t BitsInTop = 0;
        while (TopWord != 0) {
            TopWord >>= 1;
            ++BitsInTop;
        }

        return BitsInTop + LastIndex * sizeof(H) * 8;
    }
    size_t Log2Unsigned() const {
        if (IsZero()) throw std::runtime_error("Log2Unsigned(0) is undefined");
        return CountBits() - 1;
    }

    // Static, non mutating
    // return sign(abs(LHS) - abs(RHS))
    static int32_t CompareMagnitude(const BigInt& LHS, const BigInt& RHS) {
        if (LHS.Size() != RHS.Size()) {
            return LHS.Size() < RHS.Size() ? -1 : 1;
        }

        for (size_t i = LHS.Size(); i-- > 0;) {
            if (LHS[i] != RHS[i]) {
                return LHS[i] < RHS[i] ? -1 : 1;
            }
        }
        return 0;
    }
    

    // Production functions
    // Return 2^Exp
    static BigInt Power2(size_t Exp) {
        BigInt Res;
        Res.m_Data = { 1 };
        return ShiftLeft(Res, Exp);
    }

    static BigInt Negate(BigInt Val) {
        Val.m_Sign = !Val.m_Sign && !Val.IsZero();
        return Val;
    }

    static BigInt Abs(BigInt Val) {
        Val.m_Sign = false;
        return Val;
    }

    static BigInt ShiftWordsLeft(BigInt Val, size_t Amount) {
        if (Amount == 0) return Val;
        if (Val.IsZero()) return Val;

        BigInt Res;
        Res.m_Sign = Val.m_Sign;
        Res.m_Data.resize(Val.Size() + Amount, 0);

        for (size_t i = 0; i < Amount; ++i)
            Res.m_Data[i] = 0;

        for (size_t i = 0; i < Val.Size(); ++i)
            Res.m_Data[i + Amount] = Val[i];

        return Res;
    }

    static BigInt ShiftWordsRight(BigInt Val, size_t Amount) {
        if (Amount == 0) return Val;
        if (Val.IsZero()) return Val;

        if (Amount >= Val.Size()) return BigInt();

        BigInt Res;
        Res.m_Sign = Val.m_Sign;
        Res.m_Data.resize(Val.Size() - Amount, 0);

        for (size_t i = 0; i < Res.Size(); ++i)
            Res.m_Data[i] = Val[i + Amount];

        return Res;
    }

    static BigInt ShiftLeft(BigInt Val, size_t Bits) {
        if (Bits == 0 || Val.IsZero()) return Val;

        size_t WordBits = sizeof(H) * 8;
        size_t WordShift = Bits / WordBits;
        size_t BitShift = Bits % WordBits;

        BigInt Res = ShiftWordsLeft(Val, WordShift);

        if (BitShift == 0) return Res;

        H Carry = 0;
        for (size_t i = 0; i < Res.Size(); ++i) {
            F Temp = (static_cast<F>(Res[i]) << BitShift) | Carry;
            Res.m_Data[i] = lsb(Temp);
            Carry = msb(Temp);
        }
        if (Carry != 0) Res.m_Data.push_back(Carry);

        return Res;
    }

    static BigInt ShiftRight(BigInt Val, size_t Bits) {
        if (Bits == 0 || Val.IsZero()) return Val;

        size_t WordBits = sizeof(H) * 8;
        size_t WordShift = Bits / WordBits;
        size_t BitShift = Bits % WordBits;

        if (WordShift >= Val.Size()) return BigInt();

        BigInt Res = ShiftWordsRight(Val, WordShift);

        if (BitShift == 0) return Res;

        H Carry = 0;
        for (int64_t i = Res.Size() - 1; i >= 0; --i) {
            F Temp = (static_cast<F>(Res[i]) >> BitShift) | (static_cast<F>(Carry) << (WordBits - BitShift));
            Carry = Res[i] & ((1ULL << BitShift) - 1); // save bits that fell off
            Res.m_Data[i] = lsb(Temp);
        }

        Res.normalize();
        return Res;
    }

    static BigInt Add(const BigInt& LHS, const BigInt& RHS) {
        if (LHS.IsZero()) {
            return RHS;
        } else if (RHS.IsZero()) {
            return LHS;
        }

        const size_t Size = std::max(LHS.Size(), RHS.Size());

        BigInt Res;
        Res.m_Data.resize(Size, 0);
        
        if (LHS.m_Sign == RHS.m_Sign) {
            H Carry = 0;
            for (size_t i = 0; i < Size; ++i) {
                const F Sum
                    = static_cast<F>(LHS[i])
                    + static_cast<F>(RHS[i])
                    + Carry;
                
                Res.m_Data[i] = lsb(Sum);
                Carry = msb(Sum);
            }
            if (Carry != 0) {
                Res.m_Data.push_back(Carry);
            }
            Res.m_Sign = LHS.m_Sign;
            Res.normalize();
        } else {
            const bool Less = CompareMagnitude(LHS, RHS) < 0;
            const BigInt& Minuend = Less ? RHS : LHS;
            const BigInt& Subtrahend = Less ? LHS : RHS;

            H Borrow = 0;
            for (size_t i = 0; i < Size; ++i) {
                const F Diff
                    = static_cast<F>(Minuend[i])
                    - static_cast<F>(Subtrahend[i])
                    - Borrow;
                
                Res.m_Data[i] = lsb(Diff);
                Borrow = (Minuend[i] < Subtrahend[i] + Borrow) ? 1 : 0;
            }
            Res.m_Sign = Minuend.m_Sign;
            Res.normalize();
        }
        
        return Res;
    }

    static BigInt Multiply(const BigInt& LHS, const BigInt& RHS) {
        BigInt Res;
        Res.m_Data.resize(LHS.Size() + RHS.Size(), 0);

        for (size_t i = 0; i < LHS.Size(); ++i) {
            const F LHSTerm = static_cast<F>(LHS[i]);

            H Carry = 0;
            for (size_t j = 0; j < RHS.Size(); ++j) {
                const F Sum
                    = static_cast<F>(Res.m_Data[i + j])
                    + (LHSTerm * static_cast<F>(RHS[j]))
                    + Carry;
                
                Res.m_Data[i + j] = lsb(Sum);
                Carry = msb(Sum);
            }
            Res.m_Data[i + RHS.Size()] = Carry;
        }
        Res.m_Sign = LHS.m_Sign ^ RHS.m_Sign;
        Res.normalize();

        return Res;
    }

    static BigInt Divide(BigInt& OutRemainder, const BigInt& LHS, BigInt Divisor) {
        if (Divisor.IsZero()) throw std::runtime_error("Divide by zero");

        const bool FinalSign = LHS.m_Sign ^ Divisor.m_Sign;
        BigInt Res;

        // Make divisor negative so we can subtract it from the remainder
        Divisor.m_Sign = true;

        OutRemainder = LHS;
        OutRemainder.m_Sign = false;

        const size_t DivisorBits = Divisor.CountBits();

        while (CompareMagnitude(OutRemainder, Divisor) >= 0) {
            if (OutRemainder.CountBits() > DivisorBits + 1) {
                size_t BitDiff = OutRemainder.CountBits() - DivisorBits - 1;
                BigInt Operand = Divisor;
                Operand = ShiftLeft(Operand, BitDiff);

                Res = Add(Res, Power2(BitDiff));
                OutRemainder = Add(OutRemainder, Operand);
            } else {
                Res = Add(Res, BigInt(1));
                OutRemainder = Add(OutRemainder, Divisor);
            }

            assert(OutRemainder.m_Sign == false);

            if (OutRemainder.IsZero() || CompareMagnitude(OutRemainder, Divisor) < 0) {
                break;
            }
        }

        Res.m_Sign = FinalSign;
        OutRemainder.m_Sign = LHS.m_Sign;

        Res.normalize();
        OutRemainder.normalize();
        return Res;
    }

    // Slow!
    static BigInt Mod(const BigInt& LHS, const BigInt& RHS) {
        BigInt Res;
        Divide(Res, LHS, RHS);
        return Res;
    }
    
    static BigInt GCD(BigInt LHS, BigInt RHS) {
        while (RHS) {
            BigInt Temp = RHS;
            RHS = Mod(LHS, RHS);
            LHS = Temp;
        }
        return LHS;
    }


    // Serde
    static BigInt FromString(const std::string& Str) {
        BigInt Res;
        Res.m_Sign = false;

        size_t Start = 0;
        if (!Str.empty() && Str[0] == '-') {
            Res.m_Sign = true;
            Start = 1;
        }

        for (size_t i = Start; i < Str.size(); ++i) {
            char c = Str[i];
            if (c < '0' || c > '9') throw std::runtime_error("Invalid digit");
            Res = Multiply(Res, BigInt(10));
            Res = Add(Res, BigInt(c - '0'));
        }

        return Res;
    }

    static std::string ToString(const BigInt& Val) {
        if (Val.IsZero()) return "0";

        BigInt Tmp = Val;
        Tmp.m_Sign = false;

        std::string Str;
        BigInt Ten = BigInt(10);
        BigInt Remainder;

        while (!Tmp.IsZero()) {
            Tmp = Divide(Remainder, Tmp, Ten);
            Str.push_back('0' + static_cast<char>(Remainder[0]));
        }

        if (Val.m_Sign) Str.push_back('-');
        std::reverse(Str.begin(), Str.end());
        return Str;
    }

    static std::string ToHexString(const BigInt& Val) {
        if (Val.IsZero()) return "0";

        static const char HexDigits[] = "0123456789ABCDEF";
        std::string Str;
        if (Val.m_Sign) Str.push_back('-');
        Str += "0x";

        bool started = false;
        int64_t totalBits = Val.Size() * sizeof(H) * 8;

        for (int64_t bit = totalBits - 4; bit >= 0; bit -= 4) {
            size_t wordIndex = bit / (sizeof(H) * 8);
            size_t bitIndex = bit % (sizeof(H) * 8);

            H word = Val[wordIndex];
            char hexChar = HexDigits[(word >> bitIndex) & 0xF];

            if (!started && hexChar == '0') continue; // skip leading zeros
            started = true;
            Str.push_back(hexChar);
            //if (bit % (sizeof(H) * 4) == 0) Str.push_back('_');
        }

        return Str;
    }

    
    // Operators
    BigInt operator+() const { return *this; }
    BigInt operator-() const { return Negate(*this); }
    BigInt operator+(const BigInt& Other) const { return Add(*this, Other); }
    BigInt operator-(const BigInt& Other) const { return Add(*this, Negate(Other)); }
    BigInt operator*(const BigInt& Other) const { return Multiply(*this, Other); }
    BigInt operator/(const BigInt& Other) const { BigInt _; return Divide(_, *this, Other); }
    BigInt operator%(const BigInt& Other) const { BigInt Res; Divide(Res, *this, Other); return Res; }
    bool operator==(const BigInt& Other) const { return (m_Sign == Other.m_Sign) && (CompareMagnitude(*this, Other) == 0); }
    bool operator!=(const BigInt& Other) const { return !(*this == Other); }
    bool operator<(const BigInt& Other) const {
        if (Sign() > Other.Sign()) return false;
        if (Sign() < Other.Sign()) return true;
        return Sign() >= 0 ? CompareMagnitude(*this, Other) < 0 : CompareMagnitude(*this, Other) > 0;
    }
    bool operator>(const BigInt& Other) const { return Other < *this; }
    bool operator<=(const BigInt& Other) const { return !(Other < *this); }
    bool operator>=(const BigInt& Other) const { return !(*this < Other); }
    explicit operator bool() const { return !IsZero(); }
};
