# Radical isogenies on every model: Magma code

Magma implementation of the algorithms in

> Thomas Decru, *Radical isogenies on every model*, 2026.

Each algorithm computes one radical isogeny step of odd degree N directly on a given curve model, using one Nth root and without inverting any field element computed during the step. A step outputs the codomain in the same model, the coordinate of a generator of the next kernel, and the images of any pushed points.

## Files

`lib/` contains the algorithms, one file per model:

| File | Contents |
|---|---|
| `core.m` | Algorithm 1 (radical core), shared by all models |
| `tate.m` | Algorithm 2: Tate normal form |
| `weierstrass.m` | Algorithm 3: short Weierstrass curves, and the codomain-only variant |
| `montgomery.m` | Algorithm 4: Montgomery curves, and the codomain-only variant |
| `edwards.m` | Algorithm 5: twisted Edwards curves |
| `wmodel.m` | Algorithm 6: curves y² = x(x² + Ax + B), and the codomain-only variant |
| `huff.m` | Algorithm 7: Huff curves |
| `kummer.m` | Algorithm 8: Kummer lines in squared theta coordinates |
| `jacobi_quartic.m` | Algorithm 9: extended Jacobi quartics |
| `legendre.m` | Algorithm 10: Legendre curves and Jacobi intersections |
| `hessian.m` | Algorithm 11: twisted Hessian curves |
| `testing.m` | helper functions for the sanity checks |

The code uses the names of the pseudocode in the paper, and a comment `// Line k` marks line k of the corresponding algorithm. Algorithm 1 is split into `CoreRadical` (lines 1–8) and `CoreCubic` (lines 9–10), since Tate normal form only needs the former.

The top-level files (`tate.m`, `weierstrass.m`, …) are sanity-check scripts, one per model, with separate scripts `jacobi_intersection.m` for Jacobi intersections and `fallbacks.m` for the exceptional branches; `run_all.m` runs all of them.

## Running the checks

The scripts load the library with relative paths, so start Magma from this directory:

```
cd magma
magma run_all.m
```

or run a single model, e.g. `magma montgomery.m`.

For each model and each N ∈ {5, 7, 9, 11, 13, 15} (N ∈ {5, 7, 11, 13, 25} for twisted Hessian curves, which need gcd(N, 6) = 1), a script runs three consecutive radical steps over a 60-bit prime field with gcd(N, p − 1) = 1, starting from a random curve with a point of order N. Each step is compared with Magma's implementation of Vélu's formulae:

- the codomain must be isomorphic over F_p to E/⟨P⟩;
- the next generator must belong to a rational point of order N that the dual isogeny maps to ±P;
- pushed points must be the images of random points under the isogeny.

The codomain-only variants are checked in the same way. The Jacobi intersection script converts every point between J_a and the Legendre x-line (Lemma 4.23) and checks it on J_a itself. The Kummer line script also uses the doubling and differential addition of Karati and Sarkar in squared theta coordinates: it checks that the conversion to the Montgomery model respects them, that the kernel generator and the next generator have order N on their Kummer lines, and that pushing commutes with doubling.

The script `fallbacks.m` runs the exceptional branches on curves that trigger them: x(P′) = 0 in Algorithm 3, and the case δ = 0 of Proposition 4.10 in Algorithms 4 and 6, i.e. steps that land on a curve y² = x(x − r)(x − 2r).

The scripts are sanity checks and are not optimized for speed. They were run with Magma V2.29.
