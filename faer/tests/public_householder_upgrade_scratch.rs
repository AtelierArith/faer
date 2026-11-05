//! Caller-sized upgrade scratch preserves scalar/layout/previous-factor
//! contracts.
mod common;
use common::*;
use faer::linalg::householder::{
	upgrade_householder_factor, upgrade_householder_factor_scratch,
};
use reborrow::*;
trait Value: ComplexField {
	fn value(re: f64, im: f64) -> Self;
}
impl Value for f32 {
	fn value(re: f64, _: f64) -> Self {
		re as f32
	}
}
impl Value for f64 {
	fn value(re: f64, _: f64) -> Self {
		re
	}
}
impl Value for c32 {
	fn value(re: f64, im: f64) -> Self {
		Self::new(re as f32, im as f32)
	}
}
impl Value for c64 {
	fn value(re: f64, im: f64) -> Self {
		Self::new(re, im)
	}
}
fn gram<T: ComplexField>(a: MatRef<'_, T>, i: usize, j: usize) -> T {
	let mut sum = zero::<T>();
	for row in 0..a.nrows() {
		sum += conj(&a[(row, i)]) * &a[(row, j)];
	}
	sum
}
fn check<T: Value>(mut par: Par<'_>) {
	for (m, n) in [
		(0, 0),
		(1, 1),
		(3, 2),
		(96, 9),
		(384, 32),
		(2048, 16),
		(96, 64),
	] {
		for row_major in [false, true] {
			for prev in [1, 8] {
				let b = n.max(1);
				let basis = Mat::<T>::from_fn(m, n, |i, j| {
					if i < j {
						zero()
					} else if i == j {
						one()
					} else {
						T::value(0.02 * sample(i, j), 0.01 * sample(j, i))
					}
				});
				let mut input = if row_major {
					basis.transpose().to_owned()
				} else {
					basis.clone()
				};
				let mut out = Mat::<T>::full(n, n, T::value(19.0, 0.0));
				let logical_basis = if row_major {
					input.as_mut().transpose_mut()
				} else {
					input.as_mut()
				};
				let mut logical_out = if row_major {
					out.as_mut().transpose_mut()
				} else {
					out.as_mut()
				};
				for j in 0..n {
					for i in 0..=j {
						if i == j {
							logical_out[(i, j)] =
								gram(basis.as_ref(), i, j) * from_f64::<T>(0.5);
						} else if i / prev == j / prev {
							logical_out[(i, j)] = gram(basis.as_ref(), i, j);
						} else {
							logical_out[(i, j)] = T::value(23.0, 0.0);
						}
					}
				}
				let req = upgrade_householder_factor_scratch::<T>(
					m,
					n,
					b,
					prev,
					par.shape(),
				);
				let mut memory = MemBuffer::new(req);
				upgrade_householder_factor(
					logical_out.rb_mut(),
					logical_basis.rb(),
					b,
					prev,
					par.reborrow(),
					MemStack::new(&mut memory),
				);
				let bound = eps::<T::Real>()
					* from_f64::<T::Real>(32.0 * m.max(1) as f64);
				for j in 0..n {
					for i in 0..n {
						if i > j {
							assert_eq!(
								logical_out[(i, j)],
								T::value(19.0, 0.0)
							);
							continue;
						}
						let expected = gram(basis.as_ref(), i, j)
							* if i == j { from_f64::<T>(0.5) } else { one() };
						assert!(
							abs(&(logical_out[(i, j)].clone() - expected))
								<= bound,
							"({m},{n}) ({i},{j}) prev={prev}"
						);
					}
				}
				if matches!(par.shape(), ParShape::Seq) || prev >= 8 || n <= 8 {
					assert_eq!(req.size_bytes(), 0);
				}
			}
		}
	}
}
#[test]
fn caller_scratch_matches_four_scalar_layouts() {
	check::<f32>(Par::Seq);
	check::<f64>(Par::Seq);
	check::<c32>(Par::Seq);
	check::<c64>(Par::Seq);
	#[cfg(feature = "threadteam")]
	faer::team::with_team(2, |mut par| {
		check::<f32>(par.reborrow());
		check::<f64>(par.reborrow());
		check::<c32>(par.reborrow());
		check::<c64>(par);
	});
}
