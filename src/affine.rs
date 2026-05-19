use std::marker::PhantomData;

use subtle::ConditionallySelectable;

use crate::{
    U1024,
    prime_field::{PrimeFieldConfig as FieldConfig, PrimeFieldElement as FieldElement},
};

pub trait SWCurveConfig: 'static + Copy + Clone + Eq + PartialEq + std::fmt::Debug {
    type BaseField: FieldConfig + std::fmt::Debug;

    type ScalarField: FieldConfig;

    const COEFF_A: FieldElement<Self::BaseField>;

    const COEFF_B: FieldElement<Self::BaseField>;

    const ORDER: U1024;

    #[inline]
    fn mul_by_a(element: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        if Self::COEFF_A.is_zero() {
            FieldElement::<Self::BaseField>::zero()
        } else {
            element * Self::COEFF_A
        }
    }

    #[inline]
    fn add_b(element: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        element + Self::COEFF_B
    }

    fn generator() -> AffinePoint<Self>;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AffinePoint<C: SWCurveConfig> {
    pub x: FieldElement<C::BaseField>,
    pub y: FieldElement<C::BaseField>,
    pub is_infinite: bool,
    _config: PhantomData<C>,
}

impl<C: SWCurveConfig> AffinePoint<C> {
    pub fn new(x: FieldElement<C::BaseField>, y: FieldElement<C::BaseField>) -> Self {
        let point = Self {
            x,
            y,
            is_infinite: false,
            _config: PhantomData,
        };
        assert!(point.is_on_curve());
        point
    }

    pub fn infinite() -> Self {
        let zero = FieldElement::<C::BaseField>::zero();
        Self {
            x: zero,
            y: zero,
            is_infinite: true,
            _config: PhantomData,
        }
    }

    pub fn is_on_curve(&self) -> bool {
        if self.is_infinite {
            return true;
        }

        let y2 = self.y.square();
        let x3 = self.x.square() * self.x;
        let ax = C::mul_by_a(self.x);
        let rhs = C::add_b(x3 + ax);

        y2 == rhs
    }

    pub fn neg(&self) -> Self {
        if self.is_infinite {
            return *self;
        }

        Self {
            x: self.x,
            y: -self.y,
            is_infinite: false,
            _config: PhantomData,
        }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_infinite {
            return *rhs;
        }
        if rhs.is_infinite {
            return *self;
        }
        if self.neg() == *rhs {
            return Self::infinite();
        }

        if *self == *rhs {
            return self.double();
        }

        let num = rhs.y - self.y;
        let den = rhs.x - self.x;
        let lambda = num * den.inv();

        let x3 = lambda.square() - self.x - rhs.x;
        let y3 = lambda * (self.x - x3) - self.y;

        Self::new(x3, y3)
    }

    pub fn double(&self) -> Self {
        if self.is_infinite || self.y.is_zero() {
            return Self::infinite();
        }

        let three = FieldElement::new(U1024::from(3));
        let two = FieldElement::new(U1024::from(2));

        let num = three * self.x.square() + C::COEFF_A;
        let den = two * self.y;
        let lambda = num * den.inv();

        let x3 = lambda.square() - self.x - self.x;
        let y3 = lambda * (self.x - x3) - self.y;

        Self::new(x3, y3)
    }

    pub fn mul(&self, scalar: &U1024) -> Self {
        let mut r0 = Self::infinite();
        let mut r1 = *self;

        for i in (0..1024).rev() {
            let bit = scalar.bit(i) as u8;
            Self::conditional_swap(&mut r0, &mut r1, bit.into());
            r1 = r0.add(&r1);
            r0 = r0.double();
            Self::conditional_swap(&mut r0, &mut r1, bit.into());
        }

        r0
    }
}

impl<C: SWCurveConfig> ConditionallySelectable for AffinePoint<C> {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        let a_inf = a.is_infinite as u8;
        let b_inf = b.is_infinite as u8;
        Self {
            x: FieldElement::conditional_select(&a.x, &b.x, choice),
            y: FieldElement::conditional_select(&a.y, &b.y, choice),
            is_infinite: u8::conditional_select(&a_inf, &b_inf, choice) != 0,
            _config: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime_field::{PrimeFieldConfig, PrimeFieldElement};

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestBase;

    // p = 17. R = 1 mod 17.
    impl PrimeFieldConfig for TestBase {
        const MODULUS: U1024 = U1024([17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const R2: U1024 = U1024::ONE;
        const N_PRIME: U1024 = U1024([0x0f0f0f0f0f0f0f0f; 16]);
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestScalar;
    // Curve order = 19
    impl PrimeFieldConfig for TestScalar {
        const MODULUS: U1024 = U1024([19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const R2: U1024 = U1024::ONE;
        const N_PRIME: U1024 = U1024([0x0f0f0f0f0f0f0f0f; 16]);
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestCurve;

    // y^2 = x^3 + 2x + 2 mod 17
    impl SWCurveConfig for TestCurve {
        type BaseField = TestBase;
        type ScalarField = TestScalar;

        const COEFF_A: PrimeFieldElement<Self::BaseField> =
            PrimeFieldElement::from_montgomery(U1024([
                2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]));
        const COEFF_B: PrimeFieldElement<Self::BaseField> =
            PrimeFieldElement::from_montgomery(U1024([
                2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]));
        const ORDER: U1024 = U1024([19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        fn generator() -> AffinePoint<Self> {
            AffinePoint::new(
                PrimeFieldElement::from_montgomery(U1024([
                    5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ])),
                PrimeFieldElement::from_montgomery(U1024([
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ])),
            )
        }
    }

    type P = AffinePoint<TestCurve>;
    type F = PrimeFieldElement<TestBase>;

    #[test]
    fn test_curve_point() {
        let generator = TestCurve::generator();
        assert!(generator.is_on_curve());
        assert!(!generator.is_infinite);

        let inf = P::infinite();
        assert!(inf.is_infinite);
        assert!(inf.is_on_curve());
    }

    #[test]
    fn test_point_addition() {
        let p = TestCurve::generator();

        let p2_add = p.add(&p);
        let p2_double = p.double();
        assert_eq!(p2_add, p2_double, "P + P should equal 2P");

        // 2P = (6, 3) analytically computed
        let expected_p2 = P::new(F::new(U1024::from_u64(6)), F::new(U1024::from_u64(3)));
        assert_eq!(p2_double, expected_p2);

        let inf = P::infinite();
        assert_eq!(p.add(&inf), p);
        assert_eq!(inf.add(&p), p);
        assert_eq!(p.add(&p.neg()), inf);
    }

    #[test]
    fn test_point_multiplication() {
        let p = TestCurve::generator();

        // P * 2 = 2P
        let p2_mul = p.mul(&U1024::from_u64(2));
        assert_eq!(p2_mul, p.double());

        // P * order = Inf
        let p_order = p.mul(&TestCurve::ORDER);
        assert!(p_order.is_infinite);

        // P * 1 = P
        assert_eq!(p.mul(&U1024::ONE), p);

        // P * 0 = Inf
        assert!(p.mul(&U1024::ZERO).is_infinite);
    }

    #[test]
    fn test_negation() {
        let p = TestCurve::generator();
        let neg_p = p.neg();
        assert_eq!(neg_p.x, p.x);
        assert_eq!(neg_p.y, -p.y);
        assert!(neg_p.is_on_curve());

        let inf = P::infinite();
        assert_eq!(inf.neg(), inf);
    }
}
