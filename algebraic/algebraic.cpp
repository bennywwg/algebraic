// algebraic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "bignum.hpp"
#include "fixed.hpp"
#include "complex.hpp"
#include "polynomial.hpp"

#include <map>

using I = BigInt<>;
using R = Rational<I>;
using Z = Complex<R>;
using P = Polynomial<Z>;

std::string Global_ToString(const I& Val) {
    return std::to_string(Val.Size());
}

int main() {
    I a = I::FromString("80594783298243082394983980594783298243082394983983298243082394983980594783");
    I b = I::FromString("2430823949839805947832982430823949839");
    I r;

    std::cout << I::ToString(b) << "\n";
    std::cout << I::ToString(-(-a / -b)) << "\n";

    std::cout << I::ToString(-(a / b)) << "\n";
    r = a % b;
    std::cout << I::ToString(-r) << "\n";

    std::cout << "\nRationals:\n\n";

    R n = R::FromString("8934034449838893403449340.783284934939949");
    R r1 = R::FromString("493498389340344044");
    R d = R::FromString("1234567899468291094980");
    std::cout << R::ToString(n / d) << "\n";

    std::cout << R::ToString(
        R::Pow(n, 10) / R::Pow(d, 10) - R::Pow(n / d, 10)
    ) << "\n";

    std::cout << "\nComplex:\n\n";

    Z complex = Z::MakeImag(n);

    std::cout << Z::ToString(complex) << "\n";

    std::cout << "\nPolynomials:\n\n";

    P num = P(2, 3) - P(3, 2) + P(4, 1) + P(5, 0);
    P denom = P(1, 1) + P(2, 0);

    std::cout << P::ToString(num) << " / " << P::ToString(denom) << " = \n";
    std::cout << P::ToString(num / denom) << " rem " << P::ToString(num % denom) << "\n";

    return 0;
}