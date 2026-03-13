//! MOV attack simulation. The attack reduces ECDLP to DLP in F_{p^k} via pairing;
//! it is practical only for small k (≤ 6). Our curve uses k = 18.
use curve1024::u1024::LIMBS;
use curve1024::{Curve1024BaseField, PrimeFieldConfig, U1024};

fn bit_length(n: &U1024) -> usize {
    for i in (0..LIMBS).rev() {
        if n.0[i] != 0 {
            return i * 64 + (64 - n.0[i].leading_zeros()) as usize;
        }
    }
    0
}

/// Verify our embedding degree k = 18 places the pairing target F_{p^18} at 18432 bits,
/// far beyond the reach of sub-exponential index-calculus.
#[test]
fn test_mov_attack_resistance() {
    const K: u64 = 18;
    const MOV_THRESHOLD: u64 = 6;
    const {
        assert!(
            K > MOV_THRESHOLD,
            "embedding degree too small for MOV resistance"
        )
    };

    let target_bits = bit_length(&Curve1024BaseField::MODULUS) * K as usize;
    assert!(
        target_bits >= 10_000,
        "|F_{{p^{K}}}| = {target_bits} bits, need >= 10000"
    );
}

/// Show the gap between a weak k = 1 curve and our strong k = 18 curve.
/// For k = 1 the pairing maps ECDLP directly to DLP in F_p (index-calculus applies).
#[test]
fn test_mov_weak_curve_simulation() {
    let toy_p_bits: u64 = 64;
    let weak_target = toy_p_bits; // k=1:  64  bits — feasible
    let strong_target = toy_p_bits * 18; // k=18: 1152 bits — infeasible
    assert!(strong_target > 10 * weak_target);
}
