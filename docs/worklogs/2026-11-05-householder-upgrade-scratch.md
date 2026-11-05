# Caller-owned Householder factor-upgrade workspace

The compact-WY deep reduction currently allocates a contiguous `threads*n*n`
partial array on each warmed team call (faer #93, comment 6648218975; also #168
item 5). The four-scalar standalone MWE measures one allocation with width 2 and
none with width 1. Its independent Gram checks pass. This change supplies that
same working set through the canonical upgrade API and a matching scratch query;
it does not change the partition, cost gates, reduction order, or arithmetic.

The query must preserve the dispatch costs embedded in ParShape and model the
actual recursive join budgets. Concurrent diagonal children get disjoint stack
regions; scratch from sequential phases is reused. Generic initialized partials
retain the previous live-T and drop contract. No zero-allocation or performance
claim is established until the final tests and required acceptance evidence run.

The API change reaches unpivoted/pivoted QR, Hessenberg/tridiagonal reductions,
and bidiagonalization (hence EVD/SVD). Block widths must be inputs to the
tridiagonalization/bidiagonalization scratch queries so they can describe the
actual factor outputs, rather than over-reserving for every possible width.
QR look-ahead must account for simultaneous panel and trailing-update scratch.

Predeclared acceptance scope: qr_factor/qr and pivoted QR; Hessenberg/general and
self-adjoint EVD; bidiagonalization/SVD, all four scalars and declared layouts at
1T/4T/8T, with untouched LU/solve controls. Canonical library tests, a focused
four-scalar exact-query/Gram regression, no_step_alloc and no_std checks are
required. This note does not waive campaign, hosted, or hardware requirements.
No registry publication or automatic adoption is authorized by this record.

Remaining evidence and decisions will be recorded against the tested candidate.

## Initial functional evidence

The source candidate builds with std/linalg/threadteam. The complete library run
passes 680 tests (162 pre-existing ignored), the four-scalar public exact-query/
layout test passes, and no_step_alloc passes with added deep upgrades at all
four scalar types and its 4/8/70-width token budgets. The independent release MWE
reports zero allocation at widths 1/2 for every scalar, with the same parallel
section counts and Gram errors as the pinned parent. This evidence precedes the
separate #265 correction and is not the final adoption evidence.

The allocation work exposed another valid public branch: multiple target blocks
split at block_count/2 columns instead of that count times block_size. The pinned
main 9-column/block-2/previous-1 standalone MWE exits 101 with three wrong
coefficients. It is reported separately as tensor4all/faer#265 and will be fixed
in its own commit, with the matching scratch partition corrected too.

Campaign acceptance, final no_std/doctests, final regression and hosted checks
remain unverified. The change is an experiment, not adopted.
