# Checked generalized Schur block moves

Issue: [#250](https://github.com/tensor4all/faer/issues/250).

## Contract and design

Downstream generalized `ordschur` needs to move saved QZ blocks without
recomputing the original decomposition. The existing QZ swaps and ordering
helpers are private. The public facade moves one block between row indices,
updates optional Q/Z, supports all four native scalar types and strided
exclusive mutable views, and returns the last valid position if an adjacent
swap cannot be accepted. Either row of a real conjugate pair names its entire
block. Moving down adjusts the final row for differing block sizes.

Each adjacent window has at most four rows. S and T are independently scaled
by their largest real/imaginary component, using componentwise division so a
reciprocal cannot overflow at tiny scales. Existing private QZ arithmetic runs
on this bounded normalized window, with local identity Q/Z. Before accepting,
the facade checks both pencil reconstructions, local unitarity, homogeneous
block spectra in both normalized and actual scales, pair integrity, and finite
updated stripes. A regular infinite eigenvalue retains an exact zero beta;
a nonregular 0/0 block rejects. A rejected window does not change its input;
if updating outside stripes would overflow, the saved stripes are restored.
Earlier successful adjacent moves are retained and their current row is
reported. There are no whole-matrix copies or n-dependent dense checks.

Workspace is caller-owned, O(n), sized for the largest possible adjacent blocks
and optional vectors. All arithmetic runs sequentially even with a team token;
this follows the existing ordinary Schur reorder convention. Q/Z are accumulated
into the caller's values, not overwritten with new identities.

## Decisions and constraints

Exposing the existing unchecked private reorder directly would omit rejection
and return-position guarantees. Calling it on original extreme-scale pencils
also permits overflow/underflow in its products. The facade therefore confines
the unmodified kernel to a normalized local problem and validates before commit.
Full dense reconstruction or a new QZ factorization would add unnecessary
quadratic/cubic work and would no longer be a saved-factor block move.

Exact repeated or poorly separated blocks may reject. The contract promises a
valid last-successful state, not a successful ordering for every ill-conditioned
pencil. Standardized finite QZ factors are a precondition; this is not a
validator for arbitrary nonschur input. The numerical validation is bounded,
not a complete port of LAPACK's generalized nonsymmetric eigenproblem program.

This is an API/correctness extension. Existing arithmetic and execution paths
are unchanged; new entrypoints have no existing benchmark population. No speed
claim or performance adoption is made. A benchmark companion is a separate
recommendation, not evidence of numerical correctness.

## Validation evidence

The private-entry probe reports E0603 on efe576b; the public-entry compile probe
reports E0432 on the parent fc52fc9,
establishing the missing public contract. The new public contract tests cover
empty/singleton factors, real block sizes, regular infinity, all native types,
column/row-major views, optional outputs, failure positions, reconstruction,
unitarity, exact range endpoints, clustered/Jordan-like pencils and stripe
overflow rollback. Allocation tests measure warmed-up real/complex calls and
real conjugate-pair exchanges using exactly the queried scratch.

LAPACK validation uses DGET54's coupled residual semantics, its safe minimum
and 2*n*ulp factors, plus orthogonality, triangular structure and independent
homogeneous order checks. Deterministic fixtures are original; the translated
checker retains its source notice and the pinned OpenBLAS distribution revision.
The 1,520 configured QZ cases complete without a failed call or diagnostic,
including every real 1×1/2×2 adjacent-block combination at independent pencil
scales and both move directions. The focused library selection passes 36 tests.
The full LAPACK runner passes 18 families / 425 selected tests / 159,692 planned
and attempted cases, with no known, unexpected or unexplained failures. Its
QZ population contained 944 cases; the expanded 1,520-case QZ selection checks
the final facade independently. Existing numerical paths are unchanged, so
their validation remains applicable to the final change.

Public contract coverage passes 31 tests with threadteam and the warmed-up
allocation test passes, including real 2×2 exchanges and saved-factor eigenvalue
extraction. All six new doctests run successfully. The standalone f64 program
prints `position=0 first_eigenvalue=3 residual=1.3574770910276864e-15`.
Offline runner classification passes ten tests. The linalg implementation builds
without std; the std-only build and package listing are checked separately.

The complete ordinary unit/integration suite and nightly no_std executable are
left to the inherited hosted CI: the local routine checks, full bounded LAPACK
suite and explicit no_std linalg build cover the changed surface. This is an
explicit local-checklist waiver, not a claim that those entire suites ran locally.
Clippy 1.99 encounters inherited errors in unchanged files (including the generic
`value != value` NaN implementation in faer-traits). Its diagnostics for the new
production files are reviewed separately; no existing lint is disabled globally.

The saved-factor eigenvalue API supports mutable downstream Schur metadata:
real singleton values copy S/T directly, while conjugate pairs use the same
scaled local evaluator without allocation or re-factorization. Independent
characteristic-polynomial tests cover unequal beta entries, both real types,
reversed output strides, empty/singleton cases, regular infinity, 0/0 and
independent extreme scales. The existing migration/llms index points at the
updated Julia guide; no parallel guidance path is added.

The graded-pencil regression uses analytic row/column equivalents of a rotation
with eigenvalues `center ± i`. The initial independent-factor normalization
lost a 1e-300 entry against a 1e300 entry and returned imaginary eigenvalues
near 2.4e-224 instead of ±i. Joint row/column equilibration preserves every
nonzero component during each scaling step. The real-block evaluator also
checks its roots against the balanced characteristic polynomial; unsafe
scaling, nonfinite output, loss of conjugacy or a bad residual returns
`SchurEigenvalueError` before writing that block's output. A public error test
preserves the valid prefix and untouched suffix. Homogeneous pair scales are
free, so the polynomial test normalizes those scales before evaluation rather
than squaring arbitrary huge output representatives.
