// Sanity check of Algorithm 11 (twisted Hessian curves a X^3 + Y^3 + Z^3 = d XYZ, coordinate theta = X/(Y+Z)).
// The comparison uses the Weierstrass model of Lemma 4.25 with x = k theta/(3(d theta + 3)), k = d^3 - 27a.
// Run from this directory:  magma hessian.m
load "lib/core.m"; load "lib/testing.m"; load "lib/hessian.m";
SetSeed(1);
Ns := [5, 7, 11, 13];
p := SetupPrime(Ns, 60);
K := GF(p);
function HessW(a, d)
    k := d^3 - 27*a;
    return EllipticCurve([Parent(a) | d, -d^2, -k/3, k*d/3, -k^2/27]);
end function;
function HessX(a, d, Th, Xi)
    k := d^3 - 27*a; th := Th/Xi;
    return k*th/(3*(d*th + 3));
end function;
for N in Ns do
    repeat
        a := Random(K); d := Random(K);
        E1 := HessW(a, d);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    k := d^3 - 27*a; Th := 9*P1[1]; Xi := k - 3*d*P1[1];
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        ap, dp, Thp, Xip, st := HessianStep(a, d, Th, Xi, N);
        pr := HessianPushPrecomp(st, N);
        xRs2 := [];
        for R in Rs do
            T, X := HessianPush(9*R[1], k - 3*d*R[1], st, pr, N);
            Append(~xRs2, HessX(ap, dp, T, X));
        end for;
        E2 := HessW(ap, dp);
        ok, msg := CheckStep(E1, P1, E2, HessX(ap, dp, Thp, Xip), N, Rs, xRs2);
        Report(Sprintf("twisted Hessian, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, HessX(ap, dp, Thp, Xip));
        E1 := E2; a := ap; d := dp; k := d^3 - 27*a; Th := Thp; Xi := Xip;
    end for;
    ap, dp := HessianCodomainOnly(a, d, Th, Xi, N);
    Report(Sprintf("twisted Hessian, N = %o, codomain only", N), CheckCodomain(E1, P1, HessW(ap, dp), N), "");
end for;
