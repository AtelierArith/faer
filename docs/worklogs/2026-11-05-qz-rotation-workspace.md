# Generalized-Hessenberg rotation range and AED workspace

Issues: [#256](https://github.com/tensor4all/faer/issues/256),
[#260](https://github.com/tensor4all/faer/issues/260).

The current base is main `e4c6d57e`. This PR owns the generalized-Hessenberg
rotation correction and AED Sylvester index reservation. The even-rounded
window query (#253, related to #252) and row-zero sentinel correction (#255)
are already in main and are not changes contributed by this PR.

As [the maintainer clarified](https://github.com/tensor4all/faer/pull/258#issuecomment-6638555959),
the earlier closure followed deletion of the #255 base branch; it was not a
review rejection or a decision to abandon this change. The PR now targets main
and its rebase conflicts have been resolved, preserving main's regressions
and this PR's additional CI registrations.

## Design

The old generalized-Hessenberg Givens helper forms a reciprocal of the norm or
cosine. A finite subnormal can make that reciprocal infinite, producing
nonfinite factors even for a finite 3x3 pencil with identity B. The public
four-type reproducer fails on unmodified 19053d8 and succeeds with the candidate.

Reuse existing arithmetic at the owning rotation layer. For real scalars,
`JacobiRotation::make_givens` computes the coefficients from bounded ratios;
adjust its sign convention to keep the cosine nonnegative. Compute the returned
eliminated entry by dividing by the dominant coefficient, whose magnitude is
at least 1/sqrt(2). This avoids both tiny reciprocals and individually rounded
subnormal products. For complex scalars, reuse the scaled complex `rotg` and
convert its sine convention explicitly. Keep the g=0 identity/returned-f branch.
No per-call allocation, dynamic dispatch or matrix-size loop is introduced.

The AED local Sylvester solve owns four usize permutation arrays. Its query
listed one. [#260](https://github.com/tensor4all/faer/issues/260) independently
reproduces this at n=18/blocking_threshold=15 after #253 plus the exact #255
sentinel correction: 2,112 queried bytes panic requesting a 24-byte LU record.
Changing only the array count gives 2,304 bytes and completes QZ; its retained
numerical assertion still fails at ratio 61.65 (tracked with #257). Reserve
the actual four arrays rather than depending on padding or
another scratch branch accidentally covering the shortfall. A child-module
regression executes a full 2x2/2x2 swap with the AED query's workspace and checks
both pencil reconstructions.

## Historical evidence

The results below describe the original candidate and its recorded baseline,
not a rerun of all checks on current main. Current-base evidence is separated
below.

Four public finite-pencil tests cover f32/f64/c32/c64, subnormals, zero/nonzero
components, complex phase, signs and large finite entries. They verify finite
factors, Hessenberg/triangular structure, Q/Z unitarity and both reconstructions.
Large inputs are normalized componentwise before reconstruction checks to avoid
vacuous overflowing tolerances.

Exact public-query QZ tests cover real single/double precision, n=18/78/96,
default/maximum valid windows, both/no vector accumulators and a dedicated
2-thread token. Outputs start as NaN sentinels, and final real blocks do not
overlap. Independent block polynomials check homogeneous values and finite
spectrum for nonsingular B. The n=258 recursive sweep is opt-in, is executed
locally, and is wired into both scheduled numerical workflows; its runtime
belongs outside the default PR test budget.

The related library selection passes 56 tests, including the bounded 1,520-case
QZ reorder family. Existing public generalized Schur contract/allocation tests
also pass. Final compile/format and stacked-base checks are recorded in the PR.
The original n=96 workspace MWE completes with 40,960 queried bytes after the
query correction, versus its original 40,448-byte allocation panic. Its debug
execution took minutes; that is a local runtime observation, not a speed claim.

## Remaining limits

This is a range/bounds correction, not a proof of all QZ numerical accuracy.
[#257](https://github.com/tensor4all/faer/issues/257) remains open: the n=78
maximum-window fixture has a DGET54-style coupled residual ratio 436.6 on
the recorded baseline `19053d8` versus 0.7724 from LAPACK. The local range-safety candidate
reduces the ratio to 178.6 but does not meet the diagnostic threshold 20.
No tolerance is raised and no full-window accuracy fix is claimed here.

There is no performance adoption or benchmark claim. Existing Jacobi arithmetic
is reused; real/complex routes specialize at compile time. Hosted CI still owns
the full feature/MSRV/no_std matrix and coverage. The local checks are scoped to
the changed generalized-Hessenberg/QZ paths; unrelated numerical families are
unchanged. Colmat retains its current wrapper correction until merged upstream
fixes and the original reproductions plus wrapper regressions pass.


## Parent regression integration

At the cc317ac checkpoint, #255 added its six-dimension public window-32 regression and
runs both sentinel/query tests in CI. Integration preserves those checks and
this candidate's rotation/workspace checks in the same CI step. Production
`faer/src` and `faer-traits/src` trees remain identical to 12bd4bc; the parent
update changes only tests/workflow. No gate, tolerance or timeout is weakened.

The integrated focused run passes eight tests: four rotation tests, the row-zero
regression, the six-dimension window-32 query regression and both real-precision
exact-query suites. Recursive cases retain their opt-in scheduled coverage; the
unchanged production code was already validated for those cases and on both
aarch64 numerical jobs. The new head still runs its required hosted checks.


## Parent test-runtime adjustment

The later #255 head 68ac9eb reduces its window-32 public regression from
31/32/95/96/97/191 to 31/32/96/97. The original reported n=96 case remains;
n=95 and n=191 are no longer covered by that sweep. This is a reduction in
coverage, inherited from the parent's test-runtime adjustment, rather than a
numerical repair. This candidate's separate maximum/default-window suites at
n=18/78/96 and opt-in recursive n=258 suites remain unchanged; they do not
replace the removed window-32 cases.

The merged parent changes only that test. Production sources, manifests and
Cargo.lock are identical to cc317ac and the previously validated 12bd4bc.
The integrated focused run passes all eight ordinary tests (four scalar-range,
row-zero, four-dimension window-32 and both real-precision exact-query suites).
The two recursive tests retain their scheduled opt-in coverage. Required hosted
CI runs again for this combined head; no performance claim is added.


## Main integration after #255

Maintainers merged #255 as e6bb7be and retargeted this PR to main after its
base branch was removed. The only merge conflict was the CI test command;
the resolved command retains every main regression and both rotation/workspace
suites. The merge's complete tree is identical to 470360d, so the eight local
focused tests and existing numerical evidence apply unchanged. This record is
the only subsequent content change; numerical sources, manifests and Cargo.lock
remain identical to the validated candidate. Required checks run for the new
ancestry/head; no main merge or registry publication is performed here.


## Current main compatibility

The candidate now builds on main `e4c6d57e`. The sentinel correction and
public QZ regressions are supplied by main; this PR retains only its rotation
and four-array workspace corrections, tests and validation registrations.
The rotation helper and the new public regression sources are unchanged from
`e475baf`; main's existing even-rounded AED workspace assertion is retained.
Historical validation above describes its recorded revisions; current-head
hosted validation must complete again before adoption.

Current-base focused validation passes eight tests: four scalar rotation tests,
the public row-zero regression, the public window-query regression, and both
real-precision exact-query suites. The two opt-in recursive tests were not rerun
locally on this base. Hosted checks on `99bda9dc` are tracked separately; earlier
full-suite and aarch64 results do not establish success of the rebased head.
