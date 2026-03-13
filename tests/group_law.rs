//! Elliptic curve group law: generator validity, order, and algebraic identities.
use curve1024::{AffinePoint, Curve1024Config, SWCurveConfig, U1024};

#[test]
fn test_generator_is_on_curve() {
    let g = Curve1024Config::generator();
    assert!(g.is_on_curve());
    assert!(!g.is_infinite);
}

#[test]
fn test_generator_order_is_r() {
    let g = Curve1024Config::generator();
    assert!(g.mul(&Curve1024Config::ORDER).is_infinite);
}

#[test]
fn test_point_double_equals_add_self() {
    let g = Curve1024Config::generator();
    assert_eq!(g.double(), g.add(&g));
}

#[test]
fn test_identity_element() {
    let g = Curve1024Config::generator();
    let inf = AffinePoint::<Curve1024Config>::infinite();
    assert_eq!(g.add(&inf), g);
    assert_eq!(inf.add(&g), g);
}

#[test]
fn test_point_negation() {
    let g = Curve1024Config::generator();
    assert!(g.add(&g.neg()).is_infinite);
}

#[test]
fn test_scalar_mul_matches_repeated_add() {
    let g = Curve1024Config::generator();
    assert_eq!(g.add(&g.double()), g.mul(&U1024::from(3)));
}
