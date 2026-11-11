//! Public shared-factor Pfaffian/inverse contract; mathematical oracles only.
mod common;
use common::*;
use faer::linalg::pfaffian::{
	PfaffianError, pfaffian, pfaffian_inverse, pfaffian_inverse_scratch,
	pfaffian_scratch,
};
use reborrow::*;

fn residual<T: ComplexField>(a: MatRef<'_, T>, inverse: MatRef<'_, T>) {
	let infinity_norm = |a: MatRef<'_, T>| {
		let mut maximum = zero::<T::Real>();
		for i in 0..a.nrows() {
			let mut sum = zero::<T::Real>();
			for j in 0..a.ncols() {
				sum = &sum + a[(i, j)].abs();
			}
			maximum = maximum.fmax(sum);
		}
		maximum
	};
	let identity = Mat::<T>::identity(a.nrows(), a.nrows());
	let error = (&a * &inverse - identity).norm_max();
	// n products per residual entry, normalized by the operand norms.
	let bound = eps::<T::Real>()
		* from_f64::<T::Real>(64.0 * a.nrows().max(1) as f64)
		* infinity_norm(a)
		* infinity_norm(inverse);
	assert!(
		error <= bound,
		"inverse residual {error:?}, bound {bound:?}"
	);
	let skew_error = (inverse + inverse.transpose()).norm_max();
	assert!(
		skew_error <= bound,
		"plain-transpose skew error {skew_error:?}"
	);
}

fn all_contracts<T: ComplexField>(
	entry: impl Fn(f64, f64) -> T,
	complex: bool,
	tiny: f64,
	huge: f64,
) {
	// Analytic disconnected blocks: Pf=ab and each inverse upper entry=-1/a.
	let a = Mat::from_fn(4, 4, |i, j| {
		let value = match (i.min(j), i.max(j)) {
			(0, 1) => entry(2.0, 1.0),
			(2, 3) => entry(3.0, -0.5),
			_ => entry(0.0, 0.0),
		};
		if i > j { -value } else { value }
	});
	let mut inverse = Mat::<T>::zeros(4, 4);
	let mut memory = MemBuffer::new(pfaffian_inverse_scratch::<T>(4));
	let pf = pfaffian_inverse(
		a.as_ref(),
		inverse.as_mut(),
		Par::Seq,
		MemStack::new(&mut memory),
	)
	.unwrap();
	let expected_pf = if complex {
		entry(6.5, 2.0)
	} else {
		entry(6.0, 0.0)
	};
	let expected = Mat::from_fn(4, 4, |i, j| {
		let value = match (i.min(j), i.max(j)) {
			(0, 1) if complex => entry(-0.4, 0.2),
			(0, 1) => entry(-0.5, 0.0),
			(2, 3) if complex => entry(-12.0 / 37.0, -2.0 / 37.0),
			(2, 3) => entry(-1.0 / 3.0, 0.0),
			_ => entry(0.0, 0.0),
		};
		if i > j { -value } else { value }
	});
	let tolerance = eps::<T::Real>() * from_f64::<T::Real>(64.0 * 4.0);
	assert!((&pf - &expected_pf).abs() <= &tolerance * expected_pf.abs());
	assert!(
		(&inverse - &expected).norm_l2() <= &tolerance * expected.norm_l2()
	);
	residual(a.as_ref(), inverse.as_ref());

	for n in [2, 6, 34, 64, 66, 128, 130, 256] {
		let base = Mat::from_fn(n, n, |i, j| {
			let (row, col) = if i > j { (i, j) } else { (j, i) };
			let value = if i == j {
				entry(0.0, 0.0)
			} else if col % 2 == 0 && row == col + 1 {
				entry(1.25, 0.125)
			} else {
				entry(0.02 * sample(row, col), 0.01 * sample(col, row))
			};
			if i < j { -value } else { value }
		});
		let a =
			Mat::from_fn(n, n, |i, j| base[((i + 1) % n, (j + 1) % n)].copy());
		let mut padded =
			Mat::from_fn(n + 4, n + 4, |_, _| entry(f64::NAN, f64::NAN));
		padded
			.as_mut()
			.subrows_mut(1, n)
			.subcols_mut(2, n)
			.copy_from(&a);
		let canary = entry(17.0, 0.0);
		let mut output = Mat::from_fn(n + 4, n + 4, |_, _| canary.copy());
		let mut memory = MemBuffer::new(pfaffian_inverse_scratch::<T>(n));
		let mut pf_memory = MemBuffer::new(pfaffian_scratch::<T>(n));
		let reference =
			pfaffian(a.as_ref(), Par::Seq, MemStack::new(&mut pf_memory))
				.unwrap();
		let mut run = |par: Par<'_>, reversed: bool| {
			// Negative strides and a transposed output window are accepted;
			// the exact scratch query and output buffer are reused across calls.
			let input = padded.as_ref().subrows(1, n).subcols(2, n);
			let input = if reversed {
				input.reverse_rows().reverse_cols()
			} else {
				input
			};
			let mut view = output
				.as_mut()
				.subrows_mut(1, n)
				.subcols_mut(2, n)
				.transpose_mut();
			let value = pfaffian_inverse(
				input,
				view.rb_mut(),
				par,
				MemStack::new(&mut memory),
			)
			.unwrap();
			residual(input, view.rb());
			// Reversal has n(n-1)/2 transpositions, odd iff n/2 is odd.
			let expected = if reversed && n / 2 % 2 == 1 {
				-&reference
			} else {
				reference.copy()
			};
			assert!(
				(&value - &expected).abs()
					<= eps::<T::Real>()
						* from_f64::<T::Real>(64.0 * n as f64)
						* reference.abs()
			);
		};
		run(Par::Seq, true);
		run(Par::Seq, false);
		#[cfg(feature = "threadteam")]
		faer::team::with_team(3, |par| run(par, true));
		for i in 0..n {
			for j in 0..n {
				// Immutable input preservation is an exact storage contract.
				assert_eq!(padded[(i + 1, j + 2)], a[(i, j)]);
			}
		}
		for i in 0..n + 4 {
			for j in 0..n + 4 {
				if !(1..n + 1).contains(&i) || !(2..n + 2).contains(&j) {
					// These are untouched canaries, not computed numerical results.
					assert_eq!(output[(i, j)], canary);
				}
			}
		}
	}

	for n in [0, 1, 3] {
		let a = Mat::<T>::zeros(n, n);
		let mut out = Mat::<T>::zeros(n, n);
		let mut memory = MemBuffer::new(pfaffian_inverse_scratch::<T>(n));
		let value = pfaffian_inverse(
			a.as_ref(),
			out.as_mut(),
			Par::Seq,
			MemStack::new(&mut memory),
		);
		if n == 0 {
			assert!(
				(value.unwrap() - entry(1.0, 0.0)).abs() <= eps::<T::Real>()
			);
		} else {
			assert!(matches!(value, Err(PfaffianError::OddDimension)));
		}
	}
	let nonsquare = Mat::<T>::zeros(2, 3);
	let mut out = Mat::<T>::zeros(2, 2);
	let mut memory = MemBuffer::new(pfaffian_inverse_scratch::<T>(2));
	assert!(matches!(
		pfaffian_inverse(
			nonsquare.as_ref(),
			out.as_mut(),
			Par::Seq,
			MemStack::new(&mut memory)
		),
		Err(PfaffianError::NonSquare)
	));
	for value in [0.0, f64::NAN, f64::INFINITY, tiny] {
		let a = Mat::from_fn(2, 2, |i, j| {
			if i > j {
				entry(value, 0.0)
			} else {
				entry(f64::NAN, f64::NAN)
			}
		});
		let result = pfaffian_inverse(
			a.as_ref(),
			out.as_mut(),
			Par::Seq,
			MemStack::new(&mut memory),
		);
		if value == 0.0 {
			assert!(matches!(result, Err(PfaffianError::ZeroPivot { .. })));
		} else {
			assert!(matches!(result, Err(PfaffianError::NumericalFailure)));
		}
	}
	// Ignored diagonal/upper NaNs do not contaminate a valid lower triangle.
	let a = Mat::from_fn(2, 2, |i, j| {
		if i > j {
			entry(2.0, 0.0)
		} else {
			entry(f64::NAN, f64::NAN)
		}
	});
	let value = pfaffian_inverse(
		a.as_ref(),
		out.as_mut(),
		Par::Seq,
		MemStack::new(&mut memory),
	)
	.unwrap();
	assert!(
		(value + entry(2.0, 0.0)).abs()
			<= eps::<T::Real>() * from_f64::<T::Real>(8.0)
	);
	assert!(
		(out[(0, 1)].copy() - entry(0.5, 0.0)).abs()
			<= eps::<T::Real>() * from_f64::<T::Real>(8.0)
	);
	for value in [2.0f64.powi(-20), 2.0f64.powi(20), huge] {
		let a = Mat::from_fn(2, 2, |i, j| {
			if i == j {
				entry(0.0, 0.0)
			} else if i > j {
				entry(value, 0.0)
			} else {
				entry(-value, 0.0)
			}
		});
		match pfaffian_inverse(
			a.as_ref(),
			out.as_mut(),
			Par::Seq,
			MemStack::new(&mut memory),
		) {
			Ok(pf) => {
				// Ratios/products stay O(1); no overflowing absolute norm can
				// turn this extreme-scale inverse check into a vacuous bound.
				let scale = from_f64::<T::Real>(value);
				let pf_error = (pf.real() / &scale + one::<T::Real>())
					.abs()
					.fmax((pf.imag() / &scale).abs());
				let inverse_error = (out[(0, 1)].real() * &scale
					- one::<T::Real>())
				.abs()
				.fmax((out[(0, 1)].imag() * &scale).abs());
				let bound = eps::<T::Real>() * from_f64::<T::Real>(64.0);
				assert!(
					pf_error <= bound && inverse_error <= bound,
					"scale={value:e} Pf error={pf_error:?}, inverse error={inverse_error:?}"
				);
			},
			Err(error) => {
				assert!(matches!(error, PfaffianError::NumericalFailure))
			},
		}
	}
	let small = if tiny > 1e-100 { 1e-25 } else { 1e-200 };
	let underflow = Mat::from_fn(4, 4, |i, j| {
		if i / 2 == j / 2 && i != j {
			entry(if i > j { small } else { -small }, 0.0)
		} else {
			entry(0.0, 0.0)
		}
	});
	let mut out = Mat::<T>::zeros(4, 4);
	let mut memory = MemBuffer::new(pfaffian_inverse_scratch::<T>(4));
	assert!(matches!(
		pfaffian_inverse(
			underflow.as_ref(),
			out.as_mut(),
			Par::Seq,
			MemStack::new(&mut memory)
		),
		Err(PfaffianError::NumericalFailure)
	));
}

#[test]
fn real_single_precision_contract() {
	all_contracts(|re, _| re as f32, false, 1e-44, 1e30);
}
#[test]
fn real_double_precision_contract() {
	all_contracts(|re, _| re, false, 1e-320, 1e200);
}
#[test]
fn complex_single_precision_contract() {
	all_contracts(|re, im| c32::new(re as f32, im as f32), true, 1e-44, 1e30);
}
#[test]
fn complex_double_precision_contract() {
	all_contracts(c64::new, true, 1e-320, 1e200);
}
