# Iterative callback control

## Decisions

- Extend the canonical CG/BiCGSTAB/LSMR callbacks through `IterationControl`,
  accepting both existing unit-returning closures and standard `ControlFlow<()>`.
  A separate solver variant or panic-based cancellation would duplicate the
  contract or conflate normal stopping with numerical failure.
- A requested stop returns the current iterate and native residual estimates,
  `iter_count`, and `info.stopped = true`. It does not assert convergence and
  takes priority over the post-update convergence test. Initial convergence
  requires no callback. Explicit turbofish calls now provide two type parameters
  (`::<T, _>`); inferred legacy calls continue to compile.
- BiCGSTAB's successful half-step reaches the observer exactly once (#242).
  No solver recurrence, workspace requirement or scalar kernel changes.

## Verification conclusions and constraints

- A standalone Rust identity example on main `1284351` reports one iteration,
  zero callbacks and zero residual; the final-observation assertion fails.
- Five focused tests pass: real/complex single/double precision (nonzero
  imaginary coefficients and RHS), early stop before convergence, unchanged
  continuing unit/ControlFlow callbacks, terminal BiCGSTAB observation, and
  initial convergence with empty/block RHS. One executable API example passes.
- The changed path is matrix-free callback dispatch only. It adds no allocation,
  scratch requirement or per-element branch. No timing or performance adoption
  claim is made. Full feature/coverage validation remains with hosted CI.
