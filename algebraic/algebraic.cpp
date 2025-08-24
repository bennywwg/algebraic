// algebraic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "bignum.hpp"
#include "fixed.hpp"

#include <map>

using Num = BigInt<unsigned __int128, uint64_t>;
using Rat = Rational<unsigned __int128, uint64_t>;

int main() {
    Num a = Num::FromString("80594783298243082394983980594783298243082394983983298243082394983980594783");
    Num b = Num::FromString("2430823949839805947832982430823949839");
    Num r;

    /*
    std::cout << Num::ToHexString(a) << "\n";
    std::cout << Num::ToHexString(b) << "\n";
    std::cout << Num::ToHexString(Num::Multiply(a, b)) << "\n";

    Num d = Num::Divide(r, a, b);

    return 0;
    */

    std::cout << Num::ToString(b) << "\n";
    std::cout << Num::ToString(-Num::Multiply(-a, -b)) << "\n";

    std::cout << Num::ToString(Num::Add(a, b)) << "\n";
    std::cout << Num::ToString(-Num::Divide(r, a, -b)) << "\n";
    std::cout << Num::ToString(-r) << "\n";


    std::cout << "\nRationals:\n\n";

    Rat r1 = Rat::FromString("498389340.783284934939949");
    Rat n = Rat::FromString("493498389340344044");
    Rat d = Rat::FromString("48493");
    std::cout << Rat::ToString(n / d) << "\n";

    std::map<Num, std::string> Map;
    Map[a] = "abc";

    return 0;
}