use super::*;
use dyn_stack::MemBuffer;

#[test]
fn scratch_shortage_panics_before_mutation() {
	let mut s = crate::mat![[1.0, 0.25], [0.0, 2.0]];
	let mut t = Mat::<f64>::identity(2, 2);
	let original = (s.clone(), t.clone());
	let mut memory = [];
	let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
		let _ = reorder_schur_real(
			s.as_mut(),
			t.as_mut(),
			None,
			None,
			1,
			0,
			Par::Seq,
			MemStack::new(&mut memory),
		);
	}));
	crate::assert!(panic.is_err());
	assert_eq!((s, t), original);
}

#[test]
fn homogeneous_pairs_distinguish_infinite_and_nonregular_blocks() {
	let infinite = pair(1.0f64, 0.0, 0.0, 0.0);
	let same = pair(-2.0f64, 0.0, 0.0, 0.0);
	let nonregular = pair(0.0f64, 0.0, 0.0, 0.0);
	crate::assert!(same_pair(&infinite, &same, &f64::EPSILON));
	crate::assert!(!same_pair(&infinite, &nonregular, &f64::EPSILON));
	let finite = pair(1.0f64, 0.0, f64::from_bits(1), 0.0);
	crate::assert!(!same_pair(&infinite, &finite, &f64::EPSILON));
}

#[test]
fn no_op_requires_no_scratch() {
	let mut s = crate::mat![[1.0f64]];
	let mut t = crate::mat![[2.0f64]];
	let mut mem = MemBuffer::new(reorder_schur_scratch::<f64>(1));
	assert_eq!(
		reorder_schur_real(
			s.as_mut(),
			t.as_mut(),
			None,
			None,
			0,
			0,
			Par::Seq,
			MemStack::new(&mut mem)
		),
		Ok(0)
	);
	assert_eq!(s[(0, 0)] / t[(0, 0)], 0.5);
}

#[test]
fn invalid_shapes_and_indices_panic_before_mutation() {
	for kind in 0..3 {
		let mut s = crate::mat![[1.0f64, 0.25], [0.0, 2.0]];
		let mut t = Mat::<f64>::zeros(2, if kind == 0 { 3 } else { 2 });
		let mut q = Mat::<f64>::zeros(2, 1);
		let originals = (s.clone(), t.clone(), q.clone());
		let mut memory = [];
		let panic =
			std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
				let _ = reorder_schur_real(
					s.as_mut(),
					t.as_mut(),
					if kind == 1 { Some(q.as_mut()) } else { None },
					None,
					if kind == 2 { 2 } else { 1 },
					0,
					Par::Seq,
					MemStack::new(&mut memory),
				);
			}));
		crate::assert!(panic.is_err());
		assert_eq!((s, t, q), originals);
	}
}
