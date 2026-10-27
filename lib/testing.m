// Helpers for the sanity checks: a prime field, kernel points, and a comparison with Magma's own isogenies.

// A prime p of the given size with gcd(N, p - 1) = 1 for all N in Ns, so that N-th roots are unique.
function SetupPrime(Ns, bits)
    repeat
        p := RandomPrime(bits);
    until p gt 3 and &and[Gcd(N, p - 1) eq 1 : N in Ns];
    return p;
end function;

// A point of order N (prime) on E, assuming N divides #E.
function KernelPoint(E, N)
    n := #E;
    repeat
        P := (n div N) * Random(E);
    until P ne E!0;
    return P;
end function;

function KernelPolynomial(E, P, N)
    R<X> := PolynomialRing(BaseRing(E));
    return &*[X - (i*P)[1] : i in [1..(N - 1) div 2]];
end function;

// E2 and xP2 are the codomain and next generator computed by an algorithm, in Weierstrass form, and
// xRs2 are the x-coordinates computed for the images of the points Rs.  Checks that E2 is K-isomorphic
// to the Velu codomain E1/<P1>, that xP2 belongs to a point P2 of order N with dual(P2) = +-P1, and
// that the pushed x-coordinates are those of the images of the Rs.
function CheckStep(E1, P1, E2, xP2, N, Rs, xRs2)
    Ev, phi := IsogenyFromKernel(E1, KernelPolynomial(E1, P1, N));
    ok, iso := IsIsomorphic(Ev, E2);
    if not ok then return false, "codomain"; end if;
    ok, P2 := IsPoint(E2, xP2);
    if not ok or Order(P2) ne N then return false, "next generator"; end if;
    Q := DualIsogeny(phi)(Inverse(iso)(P2));
    if Q ne P1 and Q ne -P1 then return false, "not distinguished"; end if;
    for j in [1..#Rs] do
        if iso(phi(Rs[j]))[1] ne xRs2[j] then return false, "pushed point"; end if;
    end for;
    return true, "";
end function;

// Codomain only: E2 must be K-isomorphic to E1/<P1>.
function CheckCodomain(E1, P1, E2, N)
    Ev := IsogenyFromKernel(E1, KernelPolynomial(E1, P1, N));
    ok := IsIsomorphic(Ev, E2);
    return ok;
end function;

procedure Report(name, ok, msg)
    if ok then printf "%o: ok\n", name; else printf "%o: FAILED (%o)\n", name, msg; end if;
    assert ok;
end procedure;
