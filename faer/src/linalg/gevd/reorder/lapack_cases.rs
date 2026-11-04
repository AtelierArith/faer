// DGET54 checker adapted from LAPACK 3.12.0 TESTING/EIG/dget54.f,
// distributed by OpenBLAS 791b02a1a00c996b66a894f6df93421c53b62d80:
// https://github.com/OpenMathLib/OpenBLAS/blob/791b02a1a00c996b66a894f6df93421c53b62d80/lapack-netlib/TESTING/EIG/dget54.f
// Copyright (c) 1992-2023 The University of Tennessee and The University
// of Tennessee Research Foundation. All rights reserved.
// Copyright (c) 2000-2023 The University of California Berkeley.
// Copyright (c) 2006-2023 The University of Colorado Denver.
// SPDX-License-Identifier: BSD-3-Clause
// See COPYING.LAPACK.BSD.
// Changes: native real/complex arithmetic (adjoints for complex), a common
// component scale prevents overflowing the checker itself, and deterministic
// analytic pencils replace LAPACK's generator. The 2*n*ulp denominator and
// safe-minimum floor are retained. Ordering is an independent homogeneous
// cross-product oracle, not the production block eigenvalue helper.
fn get54(
	a: MatRef<'_, Scalar>,
	b: MatRef<'_, Scalar>,
	s: MatRef<'_, Scalar>,
	t: MatRef<'_, Scalar>,
	q: MatRef<'_, Scalar>,
	v: MatRef<'_, Scalar>,
) -> Real {
	let n = a.nrows();
	if n == 0 {
		return 0.0;
	}
	let mut scale = 0.0 as Real;
	for matrix in [a, b, s, t] {
		for j in 0..n {
			for i in 0..n {
				scale = scale
					.max(real(&matrix[(i, j)]).abs())
					.max(imag(&matrix[(i, j)]).abs());
			}
		}
	}
	if !scale.is_finite() {
		return Real::INFINITY;
	}
	if scale == 0.0 {
		return 0.0;
	}
	let normalized = |matrix: MatRef<'_, Scalar>| {
		Mat::from_fn(n, n, |i, j| {
			scalar(real(&matrix[(i, j)]) / scale, imag(&matrix[(i, j)]) / scale)
		})
	};
	let an = normalized(a);
	let bn = normalized(b);
	let sn = normalized(s);
	let tn = normalized(t);
	let za = adjoint(v);
	let ar = product(product(q, sn.as_ref()).as_ref(), za.as_ref());
	let br = product(product(q, tn.as_ref()).as_ref(), za.as_ref());
	let ae = Mat::from_fn(n, n, |i, j| an[(i, j)] - ar[(i, j)]);
	let be = Mat::from_fn(n, n, |i, j| bn[(i, j)] - br[(i, j)]);
	let abnorm = norm_one(an.as_ref())
		.max(norm_one(bn.as_ref()))
		.max(Real::MIN_POSITIVE / scale);
	let wnorm = norm_one(ae.as_ref()).max(norm_one(be.as_ref()));
	let dimension = (2 * n) as Real;
	if abnorm > wnorm {
		(wnorm / abnorm) / (dimension * Real::EPSILON)
	} else if abnorm * scale < 1.0 {
		(wnorm.min(dimension * abnorm) / abnorm) / (dimension * Real::EPSILON)
	} else {
		(wnorm / abnorm).min(dimension) / (dimension * Real::EPSILON)
	}
}

fn qz_case(
	run: &mut Run,
	n: usize,
	kind: usize,
	sa: Real,
	sb: Real,
	layout: usize,
	par: Par,
) {
	let ctx = format!("n={n} kind={kind} scale={sa:e}/{sb:e} layout={layout}");
	let s0 = Mat::from_fn(n, n, |i, j| {
		if i == j {
			let value = match kind {
				1 => 1.0,
				2 => 1.0 + Real::EPSILON.sqrt() * i as Real,
				3 => Real::EPSILON.powf(i as Real / n.max(1) as Real),
				_ => 1.0 + i as Real,
			};
			scalar(
				value,
				if <Scalar as ComplexField>::IS_REAL {
					0.0
				} else {
					0.25 * value
				},
			) * sa
		} else if i < j && kind != 4 {
			sample(i, j) * sa * 0.02
		} else {
			z(0.0)
		}
	});
	let t0 = Mat::from_fn(n, n, |i, j| {
		if i == j {
			z(if i + 1 == n { 0.0 } else { sb })
		} else if i < j && kind != 4 {
			sample(j, i) * sb * 0.01
		} else {
			z(0.0)
		}
	});
	let mut sp = guarded(n, n, layout);
	let mut tp = guarded(n, n, layout);
	interior_mut(sp.as_mut(), layout).copy_from(&s0);
	interior_mut(tp.as_mut(), layout).copy_from(&t0);
	let mut q = Mat::<Scalar>::identity(n, n);
	let mut v = q.clone();
	let mut memory = MemBuffer::new(reorder_schur_scratch::<Scalar>(n));
	let mut result = Ok(0);
	if !run.call(&ctx, || {
		result = reorder(
			interior_mut(sp.as_mut(), layout),
			interior_mut(tp.as_mut(), layout),
			Some(q.as_mut()),
			Some(v.as_mut()),
			n.saturating_sub(1),
			0,
			par,
			MemStack::new(&mut memory),
		);
		// Ill-separated finite blocks may reject. Their last-valid state is
		// checked by the same reconstruction and ordering diagnostics below.
		Ok::<_, ()>(())
	}) {
		return;
	}
	let s = interior(sp.as_ref(), layout);
	let t = interior(tp.as_ref(), layout);
	let position = match result {
		Ok(position) => position,
		Err(error) => error.position,
	};
	if kind == 0 || kind == 4 {
		run.require(&ctx, "well-separated move rejected", result.is_ok());
	}
	run.ratio(
		&ctx,
		"DGET54",
		get54(s0.as_ref(), t0.as_ref(), s, t, q.as_ref(), v.as_ref()),
		20.0,
	);
	run.ratio(
		&ctx,
		"Q unitary",
		unitary_ratio(q.as_ref(), Real::EPSILON),
		20.0,
	);
	run.ratio(
		&ctx,
		"Z unitary",
		unitary_ratio(v.as_ref(), Real::EPSILON),
		20.0,
	);
	run.require(
		&ctx,
		"guard overwritten",
		guards_ok(sp.as_ref()) && guards_ok(tp.as_ref()),
	);
	for j in 0..n {
		for i in j + 1..n {
			run.require(
				&ctx,
				"triangular structure",
				s[(i, j)] == z(0.0) && t[(i, j)] == z(0.0),
			);
		}
	}
	for i in 0..n {
		let original = if i < position {
			i
		} else if i == position {
			n - 1
		} else {
			i - 1
		};
		let a = scalar(
			real(&s0[(original, original)]) / sa,
			imag(&s0[(original, original)]) / sa,
		);
		let b = scalar(
			real(&t0[(original, original)]) / sb,
			imag(&t0[(original, original)]) / sb,
		);
		let c = scalar(real(&s[(i, i)]) / sa, imag(&s[(i, i)]) / sa);
		let d = scalar(real(&t[(i, i)]) / sb, imag(&t[(i, i)]) / sb);
		run.ratio(
			&ctx,
			"homogeneous order",
			magnitude(a * d - c * b) / Real::EPSILON / (n.max(1) as Real),
			20.0,
		);
		run.require(
			&ctx,
			"infinite eigenvalue beta",
			(t[(i, i)] == z(0.0)) == (original + 1 == n),
		);
	}
}

#[test]
fn lapack_qz_reorder_checkers() {
	let a = Mat::<Scalar>::identity(2, 2);
	let mut bad = a.clone();
	bad[(0, 1)] = z(1.0);
	assert_eq!(
		get54(
			a.as_ref(),
			a.as_ref(),
			a.as_ref(),
			a.as_ref(),
			a.as_ref(),
			a.as_ref()
		),
		0.0
	);
	crate::assert!(
		get54(
			a.as_ref(),
			a.as_ref(),
			bad.as_ref(),
			a.as_ref(),
			a.as_ref(),
			a.as_ref()
		) > 20.0
	);
}
#[test]
#[ignore = "opt-in LAPACK numerical validation (#31)"]
fn lapack_qz_reorder_cases() {
	let mut run = Run::new("qz_reorder");
	for n in [0, 1, 2, 3, 8, 33] {
		for kind in 0..5 {
			for layout in 0..4 {
				qz_case(&mut run, n, kind, 1.0, 1.0, layout, Par::Seq);
			}
		}
	}
	let mut planned = 6 * 5 * 4;
	if <Scalar as ComplexField>::IS_REAL {
		let scales = [1.0, Real::MIN_POSITIVE.sqrt(), Real::MAX.sqrt()];
		for n1 in [1, 2] {
			for n2 in [1, 2] {
				for sa in scales {
					for sb in scales {
						for layout in 0..4 {
							for down in [false, true] {
								real_block_case(
									&mut run, n1, n2, sa, sb, layout, down,
								);
							}
						}
					}
				}
			}
		}
		planned += 2 * 2 * 3 * 3 * 4 * 2;
	}
	run.finish(planned);
}
#[test]
#[ignore = "opt-in LAPACK numerical validation (#31)"]
fn lapack_qz_reorder_range() {
	let mut run = Run::new("qz_reorder");
	let scales = [
		Real::MIN_POSITIVE,
		Real::MIN_POSITIVE.sqrt(),
		1.0,
		Real::MAX.sqrt(),
		Real::MAX / 16.0,
	];
	for sa in scales {
		for sb in scales {
			for layout in 0..4 {
				qz_case(&mut run, 2, 0, sa, sb, layout, Par::Seq);
			}
		}
	}
	run.finish(5 * 5 * 4);
}
#[test]
#[ignore = "opt-in LAPACK numerical validation (#31)"]
fn lapack_qz_reorder_subnormal() {
	let mut run = Run::new("qz_reorder");
	for scale in [Real::from_bits(1), Real::MIN_POSITIVE / 2.0] {
		for layout in 0..4 {
			qz_case(&mut run, 2, 4, scale, scale, layout, Par::Seq);
		}
	}
	run.finish(2 * 4);
}
#[cfg(feature = "threadteam")]
#[test]
#[ignore = "opt-in LAPACK numerical validation (#31)"]
fn lapack_qz_reorder_parallel() {
	let mut run = Run::new("qz_reorder");
	for threads in [2, 4] {
		crate::team::with_team(threads, |mut par| {
			for layout in 0..4 {
				qz_case(&mut run, 8, 0, 1.0, 1.0, layout, par.reborrow());
			}
		});
	}
	run.finish(2 * 4);
}

// Independent determinant-polynomial coefficients, normalized per pencil
// factor. This avoids taking alpha/beta ratios, including at infinity.
fn block_polynomial(
	s: MatRef<'_, Scalar>,
	t: MatRef<'_, Scalar>,
	sa: Real,
	sb: Real,
) -> [Real; 3] {
	if s.nrows() == 1 {
		return [0.0, -real(&t[(0, 0)]) / sb, real(&s[(0, 0)]) / sa];
	}
	let a = |i, j| real(&s[(i, j)]) / sa;
	let b = |i, j| real(&t[(i, j)]) / sb;
	[
		b(0, 0) * b(1, 1) - b(0, 1) * b(1, 0),
		-(a(0, 0) * b(1, 1) + a(1, 1) * b(0, 0)
			- a(0, 1) * b(1, 0)
			- a(1, 0) * b(0, 1)),
		a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0),
	]
}
fn real_block_case(
	run: &mut Run,
	n1: usize,
	n2: usize,
	sa: Real,
	sb: Real,
	layout: usize,
	down: bool,
) {
	let n = n1 + n2;
	let ctx = format!(
		"blocks={n1}/{n2} scale={sa:e}/{sb:e} layout={layout} down={down}"
	);
	let mut a = Mat::<Scalar>::zeros(n, n);
	let mut b = a.clone();
	for j in 0..n {
		for i in 0..j {
			a[(i, j)] = z(0.02) * sample(i, j);
			b[(i, j)] = z(0.01) * sample(j, i);
		}
	}
	let mut offset = 0;
	for (index, size) in [n1, n2].into_iter().enumerate() {
		if size == 2 {
			a[(offset, offset)] = z(0.1 + 0.4 * index as Real);
			a[(offset + 1, offset + 1)] = a[(offset, offset)];
			a[(offset, offset + 1)] = z(1.0);
			a[(offset + 1, offset)] = z(-1.0);
			b[(offset, offset)] = z(1.2);
			b[(offset + 1, offset + 1)] = z(0.7);
		} else {
			a[(offset, offset)] = z(2.0 + index as Real);
			b[(offset, offset)] = z(if index == 1 { 0.0 } else { 1.0 });
		}
		offset += size;
	}
	for j in 0..n {
		for i in 0..n {
			a[(i, j)] *= sa;
			b[(i, j)] *= sb;
		}
	}
	let mut sp = guarded(n, n, layout);
	let mut tp = guarded(n, n, layout);
	interior_mut(sp.as_mut(), layout).copy_from(&a);
	interior_mut(tp.as_mut(), layout).copy_from(&b);
	let mut q = Mat::<Scalar>::identity(n, n);
	let mut v = q.clone();
	let mut mem = MemBuffer::new(reorder_schur_scratch::<Scalar>(n));
	let (from, to, position) = if down { (0, n1, n2) } else { (n1, 0, 0) };
	if !run.call(&ctx, || {
		reorder(
			interior_mut(sp.as_mut(), layout),
			interior_mut(tp.as_mut(), layout),
			Some(q.as_mut()),
			Some(v.as_mut()),
			from,
			to,
			Par::Seq,
			MemStack::new(&mut mem),
		)
		.map(|got| assert_eq!(got, position))
	}) {
		return;
	}
	let s = interior(sp.as_ref(), layout);
	let t = interior(tp.as_ref(), layout);
	run.ratio(
		&ctx,
		"DGET54 blocks",
		get54(a.as_ref(), b.as_ref(), s, t, q.as_ref(), v.as_ref()),
		20.0,
	);
	run.ratio(
		&ctx,
		"Q unitary",
		unitary_ratio(q.as_ref(), Real::EPSILON),
		20.0,
	);
	run.ratio(
		&ctx,
		"Z unitary",
		unitary_ratio(v.as_ref(), Real::EPSILON),
		20.0,
	);
	run.require(
		&ctx,
		"guards",
		guards_ok(sp.as_ref()) && guards_ok(tp.as_ref()),
	);
	for (old, new, size) in [(0, n2, n1), (n1, 0, n2)] {
		let pa = block_polynomial(
			a.submatrix(old, old, size, size),
			b.submatrix(old, old, size, size),
			sa,
			sb,
		);
		let pb = block_polynomial(
			s.submatrix(new, new, size, size),
			t.submatrix(new, new, size, size),
			sa,
			sb,
		);
		let an = pa.iter().fold(0.0 as Real, |x, y| x.max(y.abs()));
		let bn = pb.iter().fold(0.0 as Real, |x, y| x.max(y.abs()));
		for i in 0..3 {
			for j in 0..3 {
				run.ratio(
					&ctx,
					"block polynomial",
					((pa[i] / an) * (pb[j] / bn) - (pa[j] / an) * (pb[i] / bn))
						.abs() / Real::EPSILON
						/ n as Real,
					20.0,
				);
			}
		}
		if size == 2 {
			run.require(
				&ctx,
				"conjugate pair split",
				s[(new + 1, new)] != z(0.0),
			);
		}
		if size == 1 {
			run.require(
				&ctx,
				"infinite beta",
				(t[(new, new)] == z(0.0)) == (b[(old, old)] == z(0.0)),
			);
		}
	}
	for j in 0..n {
		for i in j + 1..n {
			run.require(&ctx, "T triangular", t[(i, j)] == z(0.0));
			if i > j + 1 || (i == n2 && j == n2 - 1) {
				run.require(&ctx, "S quasi-triangular", s[(i, j)] == z(0.0));
			}
		}
	}
}
