// Sanity check of Algorithm 10 on Jacobi intersections J_a: s^2 + c^2 = 1, a s^2 + d^2 = 1, with projective points
// (s : c : d : z) (Corollary 4.24).  Points are converted to the Legendre x-coordinate by x = (z + d)/(z - c), and
// back by c = (x^2 - 2x + a)/(x^2 - a), d = (x^2 - 2ax + a)/(x^2 - a) (Lemma 4.23).  Every point is checked on
// J_a itself, and the comparison with Velu uses the Legendre curve y^2 = x(x-1)(x-a), whose rational points have
// the same x-coordinates.  Run from this directory:  magma jacobi_intersection.m
load "lib/core.m"; load "lib/testing.m"; load "lib/wmodel.m"; load "lib/legendre.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
function JILegW(lam)
    return EllipticCurve([Parent(lam) | 0, -(1 + lam), 0, lam, 0]);
end function;
// The point (s : c : d : z) of J_a with Legendre x-coordinate x, if it is rational.
function JIPoint(a, x)
    c := (x^2 - 2*x + a)/(x^2 - a);
    d := (x^2 - 2*a*x + a)/(x^2 - a);
    ok, s := IsSquare(1 - c^2);
    if not ok then return false, []; end if;
    return true, [s, c, d, Parent(a)!1];
end function;
function JIOnCurve(a, P)
    return P[1]^2 + P[2]^2 eq P[4]^2 and a*P[1]^2 + P[3]^2 eq P[4]^2;
end function;
// (X : Z) = (z + d : z - c), the input of Algorithm 10.
function JIX(P)
    return P[4] + P[3], P[4] - P[2];
end function;
// A rational point of E1 that is not of order dividing 2N and has x^2 <> a, and the corresponding point of J_a.
function JIRandom(E1, a, N)
    repeat
        R := Random(E1);
    until R ne E1!0 and (2*N)*R ne E1!0 and R[1]^2 ne a;
    ok, RJ := JIPoint(a, R[1]);
    return R, ok and JIOnCurve(a, RJ), RJ;
end function;
for N in Ns do
    repeat
        a := Random(K);
    until a ne 0 and a ne 1 and #JILegW(a) mod N eq 0;
    E1 := JILegW(a);
    P1 := KernelPoint(E1, N);
    ok, PJ := JIPoint(a, P1[1]);
    ok := ok and JIOnCurve(a, PJ);
    if ok then X, Z := JIX(PJ); ok := X/Z eq P1[1]; end if;
    Report(Sprintf("Jacobi intersection, N = %o, kernel generator on J_a", N), ok, "");
    lamn := a; C := K!1;
    for step in [1..3] do
        Rs := []; RJs := []; okR := true;
        for j in [1..3] do
            R, okj, RJ := JIRandom(E1, a, N);
            Append(~Rs, R); Append(~RJs, RJ); okR := okR and okj;
        end for;
        X1, Z1 := JIX(PJ);
        lamnp, Cp, Xp, Zp, st := LegendreStep(lamn, C, X1, Z1, N);
        ap := lamnp/Cp;
        xRs2 := [];
        for RJ in RJs do
            XR, ZR := JIX(RJ);
            X, Z := LegendrePush(XR, ZR, st, N);
            okj, RJ2 := JIPoint(ap, X/Z);
            okR := okR and okj and JIOnCurve(ap, RJ2);
            Append(~xRs2, X/Z);
        end for;
        okP, PJ2 := JIPoint(ap, Xp/Zp);
        if okP then X2, Z2 := JIX(PJ2); okP := JIOnCurve(ap, PJ2) and X2/Z2 eq Xp/Zp; end if;
        E2 := JILegW(ap);
        ok, msg := CheckStep(E1, P1, E2, Xp/Zp, N, Rs, xRs2);
        if not okR then ok := false; msg := "pushed point not on J_a'"; end if;
        if not okP then ok := false; msg := "next generator not on J_a'"; end if;
        Report(Sprintf("Jacobi intersection, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, Xp/Zp);
        E1 := E2; a := ap; PJ := PJ2; lamn := lamnp; C := Cp;
    end for;
end for;
