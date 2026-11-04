# Saved sparse LU components (#235)

## Decisions

Expose `LuRef::factors` as owned CSC factors, and borrow the matching numeric and
symbolic owner through `solvers::Lu::as_ref`. This preserves saved-factor reuse
and permutation conventions without exposing private storage. Supernodal export
splits the combined diagonal panel and appends lower and transposed upper blocks;
finalized row indices already contain the global pivot permutation. Extraction
allocates at the explicit owned-output boundary and leaves the source unchanged.

## Verification conclusions and constraints

Reconstruction and solve residuals pass for f32/f64/c32/c64, three algorithm
thresholds, and dimensions 0,1,3,12,40. Complex fixtures include nonzero imaginary
coefficients. Both new public examples execute. An unrelated existing Cholesky
doctest requires additional features; the minimal-feature broad doctest run is
not globally green. The new API depends on the separate no-Rayon fix in #237.
No numerical factorization or solve kernel is changed by this API extension.
Downstream Colmat adoption, broader shapes and lifecycle validation remain with
AtelierArith/Colmat.jl#22.
