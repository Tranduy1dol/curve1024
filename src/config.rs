use crate::{AffinePoint, PrimeFieldConfig, PrimeFieldElement, SWCurveConfig, U1024};

include!(concat!(env!("OUT_DIR"), "/constants.rs"));

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve1024BaseField;

impl PrimeFieldConfig for Curve1024BaseField {
    const MODULUS: U1024 = CURVE_MODULUS;
    const R2: U1024 = CURVE_R2;
    const N_PRIME: U1024 = CURVE_N_PRIME;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve1024ScalarField;

impl PrimeFieldConfig for Curve1024ScalarField {
    const MODULUS: U1024 = CURVE_ORDER;
    const R2: U1024 = CURVE_ORDER_R2;
    const N_PRIME: U1024 = CURVE_ORDER_N_PRIME;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve1024Config;

impl SWCurveConfig for Curve1024Config {
    type BaseField = Curve1024BaseField;
    type ScalarField = Curve1024ScalarField;

    const COEFF_A: PrimeFieldElement<Curve1024BaseField> =
        PrimeFieldElement::from_montgomery(CURVE_A);
    const COEFF_B: PrimeFieldElement<Curve1024BaseField> =
        PrimeFieldElement::from_montgomery(CURVE_B);
    const ORDER: U1024 = CURVE_ORDER;

    fn generator() -> AffinePoint<Self> {
        AffinePoint::new(
            PrimeFieldElement::from_montgomery(CURVE_GENERATOR_X),
            PrimeFieldElement::from_montgomery(CURVE_GENERATOR_Y),
        )
    }
}
