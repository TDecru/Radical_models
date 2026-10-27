// Sanity check of Algorithm 10 (Legendre curves y^2 = x(x-1)(x-lambda); Jacobi intersections use the same
// x-line, Corollary 4.24).  Run from this directory:  magma legendre.m
load "lib/core.m"; load "lib/testing.m"; load "lib/wmodel.m"; load "lib/legendre.m";
SetSeed(1);
Ns := [5, 7, 11, 13];
p := SetupPrime(Ns, 60);
K := GF(p);
function LegW(lam)
    return EllipticCurve([Parent(lam) | 0, -(1 + lam), 0, lam, 0]);
end function;
for N in Ns do
    repeat
        lam := Random(K);
        E1 := LegW(lam);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    lamn := lam; C := K!1; X1 := P1[1]; Z1 := K!1;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        lamnp, Cp, Xp, Zp, st := LegendreStep(lamn, C, X1, Z1, N);
        xRs2 := [];
        for R in Rs do
            X, Z := LegendrePush(R[1], K!1, st, N);
            Append(~xRs2, X/Z);
        end for;
        E2 := LegW(lamnp/Cp);
        ok, msg := CheckStep(E1, P1, E2, Xp/Zp, N, Rs, xRs2);
        Report(Sprintf("Legendre, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, Xp/Zp);
        E1 := E2; lamn := lamnp; C := Cp; X1 := Xp; Z1 := Zp;
    end for;
    lamnp, Cp := LegendreCodomainOnly(lamn, C, X1, Z1, N);
    Report(Sprintf("Legendre, N = %o, codomain only", N), CheckCodomain(E1, P1, LegW(lamnp/Cp), N), "");
end for;
