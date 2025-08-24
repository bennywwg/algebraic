// algebraic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "bignum.hpp"

using Num = BigInt<uint64_t, uint32_t>;


int main() {
    Num a = Num::FromString("5120000000000");
    Num b = Num::FromString("2550000000000");
    Num r;

    Num d = Num::Divide(r, a, b);

    std::cout << Num::ToString(a) << "\n";
    std::cout << Num::ToString(b) << "\n";

    std::cout << Num::ToString(Num::Add(a, b)) << "\n";
    std::cout << Num::ToString(Num::Multiply(a, b)) << "\n";
    std::cout << Num::ToString(Num::Divide(r, a, b)) << "\n";
    std::cout << Num::ToString(r) << "\n";

    return 0;
}
