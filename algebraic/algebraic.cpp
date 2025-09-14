// algebraic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "bignum.hpp"
#include "rational.hpp"
#include "complex.hpp"
#include "polynomial.hpp"

#include <map>

using I = BigInt<>;
using R = Rational<I>;
using Z = Complex<R>;
using P = Polynomial<Z>;

std::string RoundedString(R Val, int RoundToDecimal) {
    R scale = R::Pow(R(10), RoundToDecimal);
    
    return R::ToString(R((Val * scale).Round()) / scale);
    //return R::ToString(scale);
}

P GenerateH(unsigned int n) {
    P result(1, 0);
    for (unsigned int k = 1; k <= n; ++k) {
        unsigned int m = 2 * k - 1;
        result *= P(1, 2) - P(m * m, 0);
    }
    return result;
}

void run() {
    I a = I::FromString("80594783298243082394983980594783298243082394983983298243082394983980594783");
    I b = I::FromString("2430823949839805947832982430823949839");
    I r;

    I LHS(-100);
    I RHS(-4);

    std::cout << I::ToString(LHS - RHS) << "\n";

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

    std::cout << "\n" << R::ToString(R(0.0001)) << "\n";
    std::cout << "\n" << R::ToString(R(1e-320), 400) << "\n";

    std::cout << "\nComplex:\n\n";

    Z complex = Z::MakeImag(n);

    std::cout << Z::ToString(complex) << "\n";

    std::cout << "\nPolynomials:\n\n";

    P num = P(2, 3) - P(3, 2) + P(4, 1) + P(5, 0);
    P denom = P(1, 1) + P(2, 0);

    std::cout << "(" << P::ToString(num) << ") / (" << P::ToString(denom) << ") = \n";
    std::cout << P::ToString(num / denom) << " rem " << P::ToString(num % denom) << "\n";

    std::cout << "\nsturm\n\n";

    P roots = GenerateH(4);

    std::cout << "P = " << P::ToString(roots) << "\n";

    auto sturm = P::MakeSturmSequence(roots);

    R cauchy = P::CauchyBounds(roots);

    auto Roots = P::EvaluateRootsInRange(sturm, -cauchy, cauchy, R(0.0001));

    for (size_t i = 0; i < Roots.size(); ++i) {
        std::cout << "Root " << i << " = " << RoundedString(Roots[i], 3) << "\n";
    }
    
    std::cout << "\n\n";
}

int main() {
    try {
        run();
    } catch (std::runtime_error er) {
        std::cout << er.what() << "\n";
        return 1;
    }

    return 0;
}