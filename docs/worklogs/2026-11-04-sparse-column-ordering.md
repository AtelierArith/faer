# Explicit sparse LU/QR column ordering

## Decisions

Keep existing COLAMD constructors unchanged and add explicit low-level symbolic
entry points taking an optional borrowed permutation. Copy both permutation
arrays into the symbolic owner before analysis; numeric factorization and
solves already consume that owner. This avoids downstream input reordering and
permutation-composition workarounds. Existing constructors delegate with None.

## Verification conclusions and constraints

Four scalar types, nonzero complex coefficients, identity/reverse/cyclic columns,
auto/simplicial/supernodal, empty/singleton/5/40 dimensions, LU reconstruction
and square/tall QR solve checks pass. Both public doctests execute. High-level LU construction, invalid permutation length tests and None/default
equivalence pass; downstream integration remains pending. The extension does not alter the
numeric kernels. The base includes PR #238 and its #237 dependency.
