# algebraic

P <- x^2 - 2 = pc_n*x^n + pc_(n-1)*x^(n-1) + ... + pc_1*x + pc_0
Q <- x^2 - 3 = qc_n*x^n + qc_(n-1)*x^(n-1) + ... + qc_1*x + qc_0
A <- x - s   = ac_n*s^n + ac_(n-1)*s^(n-1) + ... + ac_1*s + ac_0


Q(A(s))
    = qc_n      *(ac_n*s^n  + ac_(n-1)*s^(n-1) + ... + ac_1*s + ac_0)^n
    + qc_(n-1)  *(ac_n*s^n  + ac_(n-1)*s^(n-1) + ... + ac_1*s + ac_0)^(n-1)
    + ...
    + qc_1      *(ac_n*s^n  + ac_(n-1)*s^(n-1) + ... + ac_1*s + ac_0)
    + qc_0


res_x(P(x), Q(A(x)))
res_x(P(x), Q(x - s))
res_x(x^2 - 2, (x-s)^2 - 3)
res_x(x^2 - 2, x^2 - 2xs + s^2 - 3)
res_x(x^2 - 2, x^2 - (2s)x + (s^2 - 3))




TODO:
- Use "Minimal Polynomial Computation" w/ number field extension 𝑄(𝛼) linear algebra technique

