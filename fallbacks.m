// Sanity check of the exceptional branches: abar = 0 (x(P') = 0) in line 6 of Algorithm 3, and delta = 0
// (Proposition 4.10) in line 7 of Algorithms 4 and 6, on curves that trigger them.  For delta = 0 the steps start
// from the dual of an isogeny out of y^2 = x(x - r)(x - 2r), so that they land on a curve with 2A'^2 = 9B'.
// Run from this directory:  magma fallbacks.m
load "lib/core.m"; load "lib/testing.m"; load "lib/weierstrass.m"; load "lib/montgomery.m"; load "lib/wmodel.m";
SetSeed(1);
// 1. Algorithm 3, abar = 0 (x(P') = 0): p = 433, N = 7, y^2 = x^3 + 91x + 344, x(P) = 104
K := GF(433); N := 7;
E1 := EllipticCurve([K | 91, 344]); _, P1 := IsPoint(E1, K!104);
assert Order(P1) eq N;
Rs := []; repeat R := Random(E1); if (2*N)*R ne E1!0 then Append(~Rs, R); end if; until #Rs eq 5;
a4p, a6p, xPp, st := SWStep(K!91, K!344, K!104, K!1, N);
pr := SWPushPrecomp(st, N);
xRs2 := [];
for R in Rs do X, Z := SWPush(R[1], K!1, pr, N); Append(~xRs2, X/Z); end for;
ok, msg := CheckStep(E1, P1, EllipticCurve([K | a4p, a6p]), xPp, N, Rs, xRs2);
Report("short Weierstrass, abar = 0 (p = 433, N = 7)", ok and xPp eq 0, msg);
// 2. delta = 0 (Proposition 4.10): steps landing on y^2 = x(x - r)(x - 2r), i.e. 2A'^2 = 9B'
function DualSetup(K, A, B, N)
    Ep := EllipticCurve([K | 0, A, 0, B, 0]);
    Q := KernelPoint(Ep, N);
    E := IsogenyFromKernel(Ep, KernelPolynomial(Ep, Q, N));      // Velu: same a1, a2, a3
    xT := &+[(Ep![0, 0] + i*Q)[1] - (i*Q)[1] : i in [1..N - 1]];
    c := Coefficients(E);
    A2 := c[2] + 3*xT; B2 := c[4] + 2*c[2]*xT + 3*xT^2;
    assert c[5] + c[4]*xT + c[2]*xT^2 + xT^3 eq 0;
    Et := EllipticCurve([K | 0, -A2, 0, B2, 0]);                 // the dual kernel is rational on the twist
    Pt := KernelPoint(Et, N);
    return A2, B2, -Pt[1];
end function;
// x is the x-coordinate of a point of order N on y^2 = x(x^2 + Ax + B) or on its quadratic twist
function NTorsX(K, A, B, x, N)
    E := EllipticCurve([K | 0, A, 0, B, 0]); Et := EllipticCurve([K | 0, -A, 0, B, 0]);
    ok, P := IsPoint(E, x); if ok then return Order(P) eq N; end if;
    ok, P := IsPoint(Et, -x); return ok and Order(P) eq N;
end function;
for N in [7, 13, 31] do
    repeat
        p := RandomPrime(60);
    until p mod N eq N - 1 and p mod 8 in {1, 7} and #EllipticCurve([GF(p) | 0, 3, 0, 2, 0]) mod N eq 0;
    K := GF(p);
    A2, B2, x := DualSetup(K, K!3, K!2, N);
    Anp, Bnp, Cp, Xp, Zp := WStep(A2, B2, K!1, x, K!1, N);
    A3 := Anp/Cp; B3 := Bnp/Cp;
    Report(Sprintf("E_{A,B}, delta = 0, N = %o", N), 2*A3^2 eq 9*B3 and NTorsX(K, A3, B3, Xp/Zp, N), "");
    s2 := Sqrt(K!2);
    repeat
        A2, B2, x := DualSetup(K, 3/s2, K!1, N);
        sq, sb := IsSquare(B2);
    until sq;
    Anp, Cp, Xp, Zp := MontgomeryStep(A2/sb, K!1, x/sb, K!1, N);
    A3 := Anp/Cp;
    Report(Sprintf("Montgomery, delta = 0, N = %o", N), A3^2 eq 9/K!2 and NTorsX(K, A3, K!1, Xp/Zp, N), "");
end for;
