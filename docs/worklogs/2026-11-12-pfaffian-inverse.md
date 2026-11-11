# Shared Pfaffian and inverse for mvmc-rs

## Decisions

The downstream task [mvmc-rs#731](https://github.com/AtelierArith/mvmc-rs/issues/731)
needs the Pfaffian and inverse together. Retain L, D and the permutation during
the existing paired reduction, then use faer's unit-triangular inverse and
triangular product. A second LU or a second skew reduction duplicates cubic
work. The standalone Pfaffian specializes the factor-output flag away.

The public surface adds only the combined call and its scratch query. Retained
factors remain private; downstream code has no need to depend on their storage.
Input authority, plain-transpose complex arithmetic and input preservation
follow the existing Pfaffian. Outputs are unspecified on failure. No thread
team is created internally and no new parallel split is introduced.
The square work matrices reuse the module's established cache-line padding
at 1 KiB leading dimensions. Empty and odd dimensions need no scratch because
the combined call returns before allocating work matrices.

Numerical failure is an explicit result, including nonrepresentable intermediate
arithmetic. A native complex reciprocal can become zero after denominator-norm
overflow; publishing that as a successful zero inverse was reproduced at
1e30 for c32 and 1e200 for c64. A finite, nonzero reciprocal preflight rejects
that case without changing ordinary division order or relaxing the tests.
The API does not promise scaled extreme-range complex division.

The arithmetic follows the existing module's published paired factorization
and A[perm,perm] = L D L^T. The inverse uses existing faer primitives. C 1.4.0
was inspected for the downstream error/input contract, not translated into
this implementation. The existing source citations and licenses are retained.
No proprietary library internals were inspected.

## Verification conclusions and constraints

All four native scalar types have analytic block Pfaffian/inverse checks,
independent inverse residuals, input preservation, padded/negative-stride
input and transposed/reused output checks at panel boundaries. Retained
factors reconstruct the permuted matrix in a private child test. Empty,
odd, nonsquare, zero-pivot, nonfinite, inverse-overflow and product-underflow
cases distinguish their typed errors. Exact scratch queries are used.

The inherited library, integration and doctest suite passes on Linux x86_64
with Rust 1.96.0. Minimal std/linalg public contracts pass too. The warmed
allocation regression covers complete real/complex calls on held teams;
packing allocation behavior on other backends remains inherited, rather than
a new universal no-allocation claim. No parallel speedup is claimed.

The preliminary full-rebuild comparisons are retained in the downstream
[report](https://github.com/AtelierArith/mvmc-rs/tree/perf/faer-pfaffian/benchmark/pfaffian_faer/provider-prototype).
They compare the original provider revision 144621059bca3f0348c4ddb64d68428185528daf
plus a prototype, not Codeberg upstream. They are shared-host kernel probes,
not the faer-benchmark campaign's clean acceptance grid. This branch provides
an immutable downstream dependency candidate; provider-main adoption,
registry publication and any upstream proposal remain separate work.
