#pragma once

#include "bignum.hpp"
#include "rational.hpp"
#include "complex.hpp"
#include "polynomial.hpp"

#include <vector>

// C(x) = M(x) + D(x)
// C(1) = 0, M(1) = 0, D(1) = 0
// C(2) = 1, M(2) = 0, D(2) = 1
// C(4) = 2, M(4) = 0, D(4) = 2
// C(x*2^s + y) = 

namespace Cz {
    using I = BigInt<>;
    using R = Rational<I>;

    const char* const Zero = "□";
    const char* const One = "◪";
    const char* const Two = "■";

    static const double InvLog2Of3 = 1.0 / std::log2(3.0);

    struct Result {
        size_t M = 0; // Number of multiply steps
        size_t D = 0; // Number of divide steps

        // The highest bit index reached during the process
        // If you add D, you get the maximum bit index needed to represent any intermediate value
        size_t B = 0;
    };

    size_t CountDist(I Val) {
        size_t MaxBit = Val.TopBitIndex();

        bool Found = false;

        size_t Res = 0;
        for (size_t i = 0; i <= MaxBit; ++i) {
            if (Val.GetBit(i)) {
                if (Found) {
                    return Res;
                } else {
                    Found = true;
                }
            } else if (Found) {
                ++Res;
            }
        }

        return 0;
    }

    std::string AsBase3(I Val, size_t Num) {
        if (Num == 0) {
            std::cout << Zero;
            return "";
        }

        std::vector<const char*> Digits;

        while (!Val.IsZero()) {
            I Remainder = Val % I(3);
            Val /= I(3);

            if (Remainder == I(0)) {
                Digits.push_back(Zero);
            } else if (Remainder == I(1)) {
                Digits.push_back(One);
            } else {
                Digits.push_back(Two);
            }
        }

        while (Digits.size() < Num) {
            Digits.push_back(Zero);
        }

        std::string Res;

        for (size_t i = 0; i < Num; ++i) {
            Res += Digits[Num - 1 - i];
            if (i % 8 == 7 && i != Num - 1) {
                Res += " ";
            }
        }

        return Res;
    }

    std::string AsBase2(I Val, size_t Num) {
        std::string Res;

        for(size_t i = Num; i > 0; --i) {
            if (Val.GetBit(i - 1)) {
                Res += Two;
            } else {
                Res += Zero;
            }

            if ((i - 1) % 8 == 0 && i != 1) {
                Res += " ";
            }
        }

        return Res;
    }

    bool GetNextStep(I Val) {
        if (Val.GetBit(0) == 0) {
            return true;
        } else {
            return false;
        }
    }

    std::string ToString(I Val, size_t decW = 10, size_t binW = 64, size_t triW = 64) {
        std::string out;
        if (decW > 0) {
            std::string dec = I::ToString(Val);
            if (dec.size() < decW) dec = std::string(decW - dec.size(), '_') + dec;
            out += dec;
        }
        if (binW > 0) out += " -> 2" + AsBase2(Val, binW);
        if (triW > 0) out += " 3" + AsBase3(Val, triW);
        out += GetNextStep(Val) ? " ( ÷2 )" : " (×3+1)";
        out += " <- " + std::to_string(CountDist(Val));
        return out;
    }

    void Apply(I& Val) {
        if (Val.GetBit(0) == 0) {
            Val.ApplyShiftRight(1);
        } else {
            Val = Val * I(3) + I(1);
        }
    }

    static Result C(I Val) {
        Result Res = {0, 0, Val.TopBitIndex()};

        while (Val != I(1)) {
            if (GetNextStep(Val)) {
                Val.ApplyShiftRight(1);
                Res.D += 1;
            } else {
                Val = Val * I(3) + I(1);
                Res.M += 1;
                Res.B = std::max(Res.B, Val.TopBitIndex() + Res.D);
            }
        }

        Res.B += 1;

        return Res;
    }
    
    void PrintStepsInfo(I Val) {
        Result Res = C(Val);
        std::cout
            << I::ToString(Val)
            << "->C="<< (Res.M + Res.D)
            << ",M=" << Res.M
            << ",D=" << Res.D
            << ",B=" << Res.B
            << "\n";
    }

    // Sort of the inverse of log2(Val)
    // This returns the number of times Val can be multiplied by 3
    // before exceeding bit length N
    size_t Pow3Exponentiations(I Val, size_t B) {
        assert(!Val.IsZero());

        size_t TopBit = Val.TopBitIndex();
        if (TopBit >= B) {
            return 0;
        }

        size_t LowerBound = static_cast<size_t>(std::max(std::floor(InvLog2Of3 * (B - TopBit - 1)), 0.0));

        Val = Val * I::Pow(I(3), LowerBound);

        //std::cout << "Initial lower bound: " << LowerBound << ", " << TopBit << "\n";

        assert(Val.TopBitIndex() < B);

        do {
            Val = Val * I(3);

            //std::cout << "Checked\n";

            if (Val.TopBitIndex() >= B) break;

            LowerBound += 1;
        } while (Val.TopBitIndex() < B);

        return LowerBound;
    }
};

struct MutliCollatz {
    struct Entry {
        Cz::I BaseValue = 1;

        size_t Exponentiations = 0;

        // This will be updated lazily when applying
        size_t Shift = 0;

        // How many rightward shifts the parent had at the time of last update
        size_t LastUpdatedNumShifts = 0;
    };

    std::vector<Entry> Entries;

    void Apply() {

    }

    static MutliCollatz From(Cz::I Val) {
        if (Val.IsZero()) {
            throw std::runtime_error("Cannot create MultiCollatz from zero");
        }

        MutliCollatz Res;

        const size_t TopBit = Val.TopBitIndex();

        std::optional<size_t> Started;
        std::optional<size_t> LastBit;

        for (size_t i = 0; i <= TopBit + 1; ++i) {
            if (Started.has_value()) {
                if (Val.GetBit(i)) {
                    LastBit = i;
                } else if (i == (*LastBit + 2) || i > TopBit) {
                    //std::cout << "Entry from bit " << *Started << ", length " << (*LastBit - *Started + 1) << "\n";

                    Res.Entries.emplace_back(
                        Cz::I::Pow(Cz::I(2), *Started),
                        *Started,
                        *LastBit - *Started + 1
                    );

                    Started.reset();
                    LastBit.reset();
                }
            } else if (Val.GetBit(i)) {
                Started = i;
                LastBit = i;
            }
        }

        return {};
    }
};