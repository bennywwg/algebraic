// algebraic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "bignum.hpp"

using Num = BigInt<uint64_t, uint32_t>;


int main() {
    Num a = Num::FromString("805947832982430823949839");
    Num b = Num::FromString("378438373");
    Num r;

    /*
    std::cout << Num::ToHexString(a) << "\n";
    std::cout << Num::ToHexString(b) << "\n";
    std::cout << Num::ToHexString(Num::Multiply(a, b)) << "\n";

    Num d = Num::Divide(r, a, b);

    return 0;
    */

    std::cout << Num::ToString(a) << "\n";

    std::cout << Num::ToString(b) << "\n";
    std::cout << Num::ToString(Num::Multiply(a, b)) << "\n";

    std::cout << Num::ToString(Num::Add(a, b)) << "\n";
    std::cout << Num::ToString(Num::Divide(r, a, b)) << "\n";
    std::cout << Num::ToString(r) << "\n";

    return 0;
}
