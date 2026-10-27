// Sanity check of Algorithm 8 (Kummer lines K_{a^2,b^2} in squared theta coordinates [X : Z]).
// The comparison uses the Montgomery curve with A = a^2/b^2 + b^2/a^2 and x = (a^2 X - b^2 Z)/(a^2 Z - b^2 X)
// (Proposition 4.16).  Run from this directory:  magma kummer.m
load "lib/core.m"; load "lib/testing.m"; load "lib/montgomery.m"; load "lib/huff.m"; load "lib/kummer.m";
SetSeed(1);
Ns := [5, 7, 11, 13];
p := SetupPrime(Ns, 60);
K := GF(p);
function KummerW(a2, b2)
    return EllipticCurve([Parent(a2) | 0, a2/b2 + b2/a2, 0, 1, 0]);
end function;
function KummerX(a2, b2, X, Z)
    return (a2*X - b2*Z)/(a2*Z - b2*X);
end function;
for N in Ns do
    repeat
        a2 := Random(K); b2 := Random(K);
        E1 := KummerW(a2, b2);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    X1 := a2*P1[1] + b2; Z1 := b2*P1[1] + a2;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        a2p, b2p, Xp, Zp, st := KummerStep(a2, b2, X1, Z1, N);
        xRs2 := [];
        for R in Rs do
            X, Z := KummerPush(a2*R[1] + b2, b2*R[1] + a2, st, N);
            Append(~xRs2, KummerX(a2p, b2p, X, Z));
        end for;
        E2 := KummerW(a2p, b2p);
        ok, msg := CheckStep(E1, P1, E2, KummerX(a2p, b2p, Xp, Zp), N, Rs, xRs2);
        Report(Sprintf("Kummer line, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, KummerX(a2p, b2p, Xp, Zp));
        E1 := E2; a2 := a2p; b2 := b2p; X1 := Xp; Z1 := Zp;
    end for;
end for;
