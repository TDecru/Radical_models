// Sanity check of Algorithm 5 (twisted Edwards curves a x^2 + y^2 = 1 + d x^2 y^2, y-coordinates).
// Since (a : d) is projective, the curve is only determined up to a quadratic twist, which x-only (here
// y-only) arithmetic does not see; the comparison therefore uses the Montgomery x-line y^2 = x^3 + A x^2 + x
// with A = 2(a+d)/(a-d) and x = (1+y)/(1-y).  Run from this directory:  magma edwards.m
load "lib/core.m"; load "lib/testing.m"; load "lib/montgomery.m"; load "lib/edwards.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
function EdwardsW(a, d)
    return EllipticCurve([Parent(a) | 0, 2*(a + d)/(a - d), 0, 1, 0]);
end function;
for N in Ns do
    repeat
        A := Random(K);
        E1 := EllipticCurve([K | 0, A, 0, 1, 0]);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    a := A + 2; d := A - 2; x := P1[1]; Y1 := x - 1; Z1 := x + 1;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        ap, dp, Yp, Zp, st := EdwardsStep(a, d, Y1, Z1, N);
        xRs2 := [];
        for R in Rs do
            Y, Z := EdwardsPush(R[1] - 1, R[1] + 1, st, N);
            Append(~xRs2, (Z + Y)/(Z - Y));
        end for;
        E2 := EdwardsW(ap, dp);
        ok, msg := CheckStep(E1, P1, E2, (Zp + Yp)/(Zp - Yp), N, Rs, xRs2);
        Report(Sprintf("twisted Edwards, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, (Zp + Yp)/(Zp - Yp));
        E1 := E2; a := ap; d := dp; Y1 := Yp; Z1 := Zp;
    end for;
end for;
