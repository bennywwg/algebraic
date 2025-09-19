#pragma once

#include "polynomial.hpp"

template<typename P = Polynomial<>>
class PolyVal {
    using T = decltype(P::Term::Cof);

    P Expression;
    T Min;
    T Max;
    std::vector<P> Sturm;
public:
    static PolyVal Root(uint32_t Root, T Val) {
        if (Val < 0) {
            throw std::runtime_error("Can't take root of negative");
        } else if (Val.IsZero()) {
            Min = 0;
            Max = 0;
        } else {
            Min = 0;
            Max = Val;

            Expression = P(T(1), Root) - P(Val, 0);
        }

        
    }
}