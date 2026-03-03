use std::marker::PhantomData;

use crate::{
    U1024,
    prime_field::{PrimeFieldConfig as FieldConfig, PrimeFieldElement as FieldElement},
};

pub trait SWCurveConfig: 'static + Copy + Clone + Eq + PartialEq {
    type BaseField: FieldConfig;

    const COEFF_A: FieldElement<Self::BaseField>;

    const COEFF_B: FieldElement<Self::BaseField>;

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
}

#[derive(Clone, Copy, PartialEq, Eq)]
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
        let ax = C::mul_by_a(x3);
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
        let mut result = Self::infinite();
        let mut base = *self;

        for i in 0..1024 {
            let limb_idx = i / 64;
            let bit_idx = i % 64;
            if (scalar.0[limb_idx] >> bit_idx) & 1 == 1 {
                result = result.add(&base);
            }
            base = base.double();
        }

        result
    }
}
