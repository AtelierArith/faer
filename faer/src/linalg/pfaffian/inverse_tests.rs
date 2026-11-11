//! Retained-factor reconstruction, independently of the inverse calculation.
use super::*;
use crate::assert;

fn reconstruction<T: ComplexField>(entry: impl Fn(f64, f64) -> T) {
	for n in [2, 6, 34, 66] {
		// Rotate the paired blocks so the first pivot is at the far corner.
		let base = Mat::from_fn(n, n, |i, j| {
			let (row, col) = if i > j { (i, j) } else { (j, i) };
			let value = if row == col {
				entry(0.0, 0.0)
			} else if col % 2 == 0 && row == col + 1 {
				entry(1.25, 0.125)
			} else {
				entry(
					((row * 17 + col * 31) % 23) as f64 / 1024.0,
					((row * 11 + col * 7) % 19) as f64 / 2048.0,
				)
			};
			if i < j { -value } else { value }
		});
		let a =
			Mat::from_fn(n, n, |i, j| base[((i + 1) % n, (j + 1) % n)].copy());
		let mut lower = Mat::<T>::zeros(n, n);
		let mut diagonal = vec![zero::<T>(); n / 2];
		let mut permutation = vec![0; n];
		let mut memory = dyn_stack::MemBuffer::new(pfaffian_scratch::<T>(n));
		let mut check = |par: Par<'_>| {
			pfaffian_factor(
				a.as_ref(),
				lower.as_mut(),
				&mut permutation,
				&mut diagonal,
				par,
				dyn_stack::MemStack::new(&mut memory),
			)
			.unwrap();
			let d = Mat::from_fn(n, n, |i, j| {
				if i % 2 == 1 && j + 1 == i {
					diagonal[j / 2].copy()
				} else if j % 2 == 1 && i + 1 == j {
					-&diagonal[i / 2]
				} else {
					zero::<T>()
				}
			});
			let expected = Mat::from_fn(n, n, |i, j| {
				a[(permutation[i], permutation[j])].copy()
			});
			let actual = &lower * &d * lower.transpose();
			let error = (actual - &expected).norm_l2();
			let bound = eps::<T::Real>()
				* from_f64::<T::Real>(64.0 * n as f64)
				* expected.norm_l2();
			assert!(
				error <= bound,
				"n={n} reconstruction {error:?}, bound {bound:?}"
			);
		};
		check(Par::Seq);
		#[cfg(feature = "threadteam")]
		crate::team::with_team(3, |par| check(par));
	}
}

#[test]
fn retained_factors_reconstruct_all_native_scalar_types() {
	reconstruction(|re, _| re as f32);
	reconstruction(|re, _| re);
	reconstruction(|re, im| c32::new(re as f32, im as f32));
	reconstruction(c64::new);
}
