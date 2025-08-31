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
    bool m_Sign { false };
    std::vector<H> m_Data;

    static constexpr H lsb(F x) { return static_cast<H>(x); }
    static constexpr H msb(F x) { return static_cast<H>(x >> (sizeof(H) * CHAR_BIT)); }
    static constexpr size_t m_wordBits = sizeof(H) * 8;

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
    template<std::integral T>
    requires (sizeof(T) <= sizeof(H))
    BigInt(T Val) {
        using U = std::make_unsigned_t<T>;
        if constexpr (std::signed_integral<T>) {
            m_Sign = Val < 0;
            m_Data = { static_cast<U>(m_Sign ? -Val : Val) };
        } else {
            m_Sign = false;
            m_Data = { static_cast<U>(Val) };
        }
        normalize();
    }
    
    // Basic functions
    size_t Size() const {
        return m_Data.size();
    }
    int32_t Sign() const {
        return IsZero() ? 0 : (m_Sign ? -1 : 1);
    }
    H operator[](size_t Index) const {
        return Index < m_Data.size() ? m_Data[Index] : 0;
    }
    H& operator[](size_t Index) {
        if (Index >= m_Data.size()) {
            m_Data.resize(Index + 1, 0);
        }
        return m_Data[Index];
    }
    bool IsZero() const {
        return m_Data.empty();
    }
    size_t TopBitIndex() const {
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
        return TopBitIndex() - 1;
    }
    // return sign(abs(LHS) - abs(RHS))
    static int32_t DiffMagnitude(const BigInt& LHS, const BigInt& RHS) {
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
    

    // Static functions
    // Return 2^Exp
    static BigInt Power2(size_t Exp) {
        BigInt Res;
        Res.m_Data = { 1 };
        Res.ApplyShiftLeft(Exp);
        return Res;
    }
    static BigInt GCD(BigInt LHS, BigInt RHS) {
        while (RHS) {
            BigInt Temp = RHS;
            LHS.ApplyRemainder(RHS);
            RHS = std::move(LHS);
            LHS = std::move(Temp);
        }
        return LHS;
    }
    static BigInt Pow(BigInt LHS, size_t RHS) {
        BigInt Res { 1 };

        while (RHS > 0) {
            if (RHS & 1) {
                Res *= LHS;
            }
            LHS *= LHS;
            RHS >>= 1;
        }

        return Res;
    }


    // Mutating functions
    void ApplyZero() {
        m_Data.clear();
        m_Sign = false;
    }
    void ApplyAbs() {
        m_Sign = false;
    }
    void ApplyNegate() {
        m_Sign = !m_Sign && !IsZero();
    }
    void ApplySign(bool Negative) {
        m_Sign = IsZero() ? false : Negative;
    }
    void ApplyShiftWordsLeft(size_t Amount) {
        if (Amount == 0 || IsZero()) return;

        m_Data.insert(m_Data.begin(), Amount, H{});
    }
    void ApplyShiftWordsRight(size_t Amount) {
        if (Amount == 0 || IsZero()) return;

        if (Amount >= m_Data.size()) {
            m_Data.clear();
            m_Sign = false;
        } else {
            m_Data.erase(m_Data.begin(), m_Data.begin() + Amount);
        }
    }
    void ApplyShiftLeft(size_t Amount) {
        if (Amount == 0 || IsZero()) return;

        const size_t BitShift = Amount % m_wordBits;

        ApplyShiftWordsLeft(Amount / m_wordBits);

        if (BitShift == 0) {
            return;
        }

        H Carry = 0;
        for (size_t i = 0; i < m_Data.size(); ++i) {
            F Temp = (static_cast<F>(m_Data[i]) << BitShift) | Carry;
            m_Data[i] = lsb(Temp);
            Carry = msb(Temp);
        }
        if (Carry != 0) m_Data.push_back(Carry);
    }
    void ApplyShiftRight(size_t Amount) {
        if (Amount == 0 || IsZero()) return;

        const size_t BitShift = Amount % m_wordBits;

        ApplyShiftWordsRight(Amount / m_wordBits);

        if (BitShift == 0 || IsZero()) {
            return;
        }

        H Carry = 0;
        for (int64_t i = m_Data.size() - 1; i >= 0; --i) {
            F Temp = (static_cast<F>(m_Data[i]) >> BitShift) | (static_cast<F>(Carry) << (m_wordBits - BitShift));
            Carry = m_Data[i] & ((1ULL << BitShift) - 1);
            m_Data[i] = lsb(Temp);
        }

        normalize();
    }
    // Compute the remainder of *this / Divisor, and assign to this
    // Output the quotient in the final parameter, if specified
    void ApplyRemainder(const BigInt& Divisor, BigInt* OutQuotient = nullptr) {
        if (this == &Divisor || this == OutQuotient) throw std::runtime_error("Can't perform ApplyRemainder with itself as an operand");
        if (Divisor.IsZero()) throw std::runtime_error("Divide by zero");

        const bool QuotientSign = m_Sign ^ Divisor.m_Sign;
        const bool RemainderSign = m_Sign;

        if (OutQuotient) OutQuotient->ApplyZero(); 
        m_Sign = false;

        const size_t DivisorBits = Divisor.TopBitIndex();

        BigInt Operand;
        while (DiffMagnitude(*this, Divisor) >= 0) {
            const size_t RemainderBits = TopBitIndex();

            Operand = Divisor;
            Operand.ApplySign(true);

            if (RemainderBits > DivisorBits + 1) {
                const size_t BitDiff = RemainderBits - DivisorBits - 1;
                Operand.ApplyShiftLeft(BitDiff);

                if (OutQuotient) (*OutQuotient) += Power2(BitDiff);
                (*this) += Operand;
            } else {
                if (OutQuotient) (*OutQuotient) += BigInt(1);
                (*this) += Operand;
            }

            if (IsZero()) {
                break;
            }
        }

        if (OutQuotient) {
            OutQuotient->m_Sign = QuotientSign;
            OutQuotient->normalize();
        }

        m_Sign = RemainderSign;
        normalize();
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
            Res *= BigInt(10);
            Res += BigInt(c - '0');
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
            Remainder = Tmp;
            Remainder.ApplyRemainder(Ten, &Tmp);
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
    BigInt operator+() const {
        return *this;
    }
    BigInt operator-() const {
        BigInt Res = *this;
        Res.ApplyNegate();
        return Res;
    }
    BigInt& operator+=(const BigInt Other) {
        if (IsZero()) {
            *this = Other;

            return *this;
        } else if (Other.IsZero()) {
            return *this;
        }
        
        const int32_t DiffMag = DiffMagnitude(*this, Other);

        const size_t MaxSize = std::max(m_Data.size(), Other.m_Data.size());

        if (m_Sign == Other.m_Sign) {
            if (DiffMag == 0) {
                ApplyShiftLeft(1);

                return *this;
            }

            H Carry = 0;
            for (size_t i = 0; i < MaxSize; ++i) {
                const F Sum
                    = static_cast<F>((*this)[i])
                    + static_cast<F>(Other[i])
                    + Carry;
                
                m_Data[i] = lsb(Sum);
                Carry = msb(Sum);
            }
            if (Carry != 0) {
                m_Data.push_back(Carry);
            }
        } else {
            if (DiffMag == 0) {
                m_Data.clear();
                m_Sign = false;
            } else {
                BigInt Tmp;
                const BigInt* Subtrahend = &Other;
                if (DiffMag < 0) { // *this < Other, use the tmp data and swap to not underflow
                    m_Sign = !m_Sign;
                    Tmp.m_Data = std::move(m_Data);
                    m_Data = Other.m_Data;
                    Subtrahend = &Tmp;
                }

                H Borrow = 0;
                for (size_t i = 0; i < m_Data.size(); ++i) {
                    F Diff
                        = static_cast<F>((*this)[i])
                        - static_cast<F>((*Subtrahend)[i])
                        - Borrow;
                    m_Data[i] = lsb(Diff);
                    Borrow = (Diff >> (sizeof(H) * 8)) & 1;
                }

                normalize();
            }
        }
        
        return *this;
    }
    BigInt operator+(BigInt Other) const {
        Other += *this;
        return Other;
    }
    BigInt& operator-=(const BigInt Other) {
        *this += -Other;
        return *this;
    }
    BigInt operator-(BigInt Other) const {
        Other -= *this;
        Other.ApplyNegate();
        return Other;
    }
    BigInt& operator*=(const BigInt Other) {
        std::vector<H> TmpData;
        TmpData.resize(m_Data.size() + Other.m_Data.size(), 0);

        for (size_t i = 0; i < m_Data.size(); ++i) {
            const F Term = static_cast<F>(m_Data[i]);

            H Carry = 0;
            for (size_t j = 0; j < Other.m_Data.size(); ++j) {
                const F Sum
                    = static_cast<F>(TmpData[i + j])
                    + (Term * static_cast<F>(Other.m_Data[j]))
                    + Carry;
                
                TmpData[i + j] = lsb(Sum);
                Carry = msb(Sum);
            }
            TmpData[i + Other.m_Data.size()] = Carry;
        }
        
        m_Data = std::move(TmpData);
        m_Sign = m_Sign ^ Other.m_Sign;
        normalize();

        return *this;
    }
    BigInt operator*(BigInt Other) const {
        Other *= *this;
        return Other;
    }
    BigInt& operator/=(const BigInt Other) {
        BigInt Quotient;
        ApplyRemainder(Other, &Quotient);
        *this = Quotient;
        return *this;
    }
    BigInt operator/(const BigInt& Other) const {
        BigInt Quotient;
        BigInt Remainder = *this;
        Remainder.ApplyRemainder(Other, &Quotient);
        return Quotient;
    }
    BigInt& operator%=(const BigInt Other) {
        ApplyRemainder(Other);
        return *this;
    }
    BigInt operator%(const BigInt& Other) const {
        BigInt Res = *this;
        Res.ApplyRemainder(Other);
        return Res;
    }
    bool operator==(const BigInt& Other) const {
        return (m_Sign == Other.m_Sign) && (DiffMagnitude(*this, Other) == 0);
    }
    bool operator!=(const BigInt& Other) const {
        return !(*this == Other);
    }
    bool operator<(const BigInt& Other) const {
        if (Sign() > Other.Sign()) return false;
        if (Sign() < Other.Sign()) return true;
        return Sign() >= 0 ? DiffMagnitude(*this, Other) < 0 : DiffMagnitude(*this, Other) > 0;
    }
    bool operator>(const BigInt& Other) const {
        return Other < *this;
    }
    bool operator<=(const BigInt& Other) const {
        return !(Other < *this);
    }
    bool operator>=(const BigInt& Other) const {
        return !(*this < Other);
    }
    explicit operator bool() const {
        return !IsZero();
    }
};
