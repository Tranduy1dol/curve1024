//! Curve parameter correctness: bit sizes and NTT two-adicity.
use curve1024::u1024::LIMBS;
use curve1024::{Curve1024BaseField, Curve1024Config, PrimeFieldConfig, SWCurveConfig, U1024};

fn two_adicity(n: &U1024) -> u32 {
    for i in 0..LIMBS {
        if n.0[i] != 0 {
            return i as u32 * 64 + n.0[i].trailing_zeros();
        }
    }
    1024
}

fn bit_length(n: &U1024) -> usize {
    for i in (0..LIMBS).rev() {
        if n.0[i] != 0 {
            return i * 64 + (64 - n.0[i].leading_zeros()) as usize;
        }
    }
    0
}

#[test]
fn test_modulus_bit_length() {
    assert_eq!(bit_length(&Curve1024BaseField::MODULUS), 1024);
}

#[test]
fn test_order_bit_length() {
    assert_eq!(bit_length(&Curve1024Config::ORDER), 512);
}

#[test]
fn test_base_field_ntt_two_adicity() {
    let s = two_adicity(&Curve1024BaseField::MODULUS.borrowing_sub(&U1024::ONE).0);
    assert!(s >= 32, "two_adicity(p-1) = {s}, need >= 32");
}

#[test]
fn test_scalar_field_ntt_two_adicity() {
    let s = two_adicity(&Curve1024Config::ORDER.borrowing_sub(&U1024::ONE).0);
    assert!(s >= 32, "two_adicity(r-1) = {s}, need >= 32");
}
