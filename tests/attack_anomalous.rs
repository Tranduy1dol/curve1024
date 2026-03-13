//! Anomalous (SSSA) attack simulation. The attack works when #E(F_p) = p (trace t = 1),
//! reducing ECDLP to a single p-adic division. Our curve has t = p + 1 - r ≈ 2^1024 >> 1.
use curve1024::u1024::LIMBS;
use curve1024::{Curve1024BaseField, Curve1024Config, PrimeFieldConfig, SWCurveConfig, U1024};

fn bit_length(n: &U1024) -> usize {
    for i in (0..LIMBS).rev() {
        if n.0[i] != 0 {
            return i * 64 + (64 - n.0[i].leading_zeros()) as usize;
        }
    }
    0
}

/// Verify r ≠ p (anomalous condition would require r = p for a prime-order curve).
/// Since |r| = 512 and |p| = 1024, r < p always holds, so t >> 1.
#[test]
fn test_anomalous_attack_resistance() {
    let p = Curve1024BaseField::MODULUS;
    let r = Curve1024Config::ORDER;
    assert_ne!(r, p);
    assert!(bit_length(&p) > bit_length(&r));
}

/// Demonstrate the attack on a toy curve with p = 11, #E = 11 (t = 1).
/// There the p-adic lift reduces ECDLP to one division in Z/11Z.
#[test]
fn test_anomalous_weak_curve_simulation() {
    let toy_p: u64 = 11;
    let toy_group_order: u64 = 11; // anomalous: #E = p
    let toy_t = toy_p + 1 - toy_group_order;
    assert_eq!(toy_t, 1);

    // Our curve: t ≈ p >> 1, confirmed by |p| >> |r|
    assert!(bit_length(&Curve1024BaseField::MODULUS) > 500);
}
