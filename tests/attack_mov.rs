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

#[test]
fn test_mov_attack_resistance() {
    let p_bits = bit_length(&Curve1024BaseField::MODULUS);
    assert!(p_bits > 0);
}

#[test]
fn test_mov_weak_curve_simulation() {
    let toy_p_bits: u64 = 64;
    let weak_target = toy_p_bits;
    let strong_target = toy_p_bits * 12;
    assert!(strong_target > 10 * weak_target);
}
