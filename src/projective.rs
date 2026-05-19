use std::marker::PhantomData;

use subtle::ConditionallySelectable;

use crate::{AffinePoint, PrimeFieldElement, SWCurveConfig, U1024};

#[derive(Clone, Copy)]
pub struct ProjectivePoint<C: SWCurveConfig> {
    x: PrimeFieldElement<C::BaseField>,
    y: PrimeFieldElement<C::BaseField>,
    z: PrimeFieldElement<C::BaseField>,
    _config: PhantomData<C>,
}

impl<C: SWCurveConfig> ProjectivePoint<C> {
    pub fn from_affine(value: &AffinePoint<C>) -> Self {
        if value.is_infinite {
            return Self::infinity();
        }

        Self {
            x: value.x,
            y: value.y,
            z: PrimeFieldElement::<C::BaseField>::one(),
            _config: PhantomData,
        }
    }

    pub fn to_affine(&self) -> AffinePoint<C> {
        if self.z.is_zero() {
            return AffinePoint::<C>::infinity();
        }

        let z_inv = self.z.inv();
        let z_inv_square = z_inv.square();
        let z_inv_cube = z_inv_square * z_inv;

        let x = self.x * z_inv_square;
        let y = self.y * z_inv_cube;

        AffinePoint::<C>::new(x, y)
    }
}

impl<C: SWCurveConfig> ProjectivePoint<C> {
    pub fn infinity() -> Self {
        Self {
            x: PrimeFieldElement::<C::BaseField>::one(),
            y: PrimeFieldElement::<C::BaseField>::one(),
            z: PrimeFieldElement::<C::BaseField>::zero(),
            _config: PhantomData,
        }
    }

    pub fn double(&self) -> Self {
        if self.z.is_zero() {
            return *self;
        }

        let two = PrimeFieldElement::<C::BaseField>::new(U1024::from(2));
        let three = PrimeFieldElement::<C::BaseField>::new(U1024::from(3));
        let four = PrimeFieldElement::<C::BaseField>::new(U1024::from(4));
        let eight = PrimeFieldElement::<C::BaseField>::new(U1024::from(8));

        let a = self.y.square();
        let b = four * self.x * a;
        let c = eight * a.square();
        let d = three * self.x.square() + C::mul_by_a(self.z.square().square());

        let x = d.square() - two * b;
        Self {
            x,
            y: d * (b - x) - c,
            z: two * self.y * self.z,
            _config: PhantomData,
        }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        if self.z.is_zero() {
            return *rhs;
        }
        if rhs.z.is_zero() {
            return *self;
        }

        let u1 = self.x * rhs.z.square();
        let u2 = rhs.x * self.z.square();

        let s1 = self.y * rhs.z.square() * rhs.z;
        let s2 = rhs.y * self.z.square() * self.z;

        let h = u2 - u1;
        let r = s2 - s1;

        if h.is_zero() {
            if r.is_zero() {
                return self.double();
            } else {
                return Self::infinity();
            }
        }

        let h2 = h.square();
        let h3 = h2 * h;
        let two = PrimeFieldElement::<C::BaseField>::new(U1024::from(2));

        let x = r.square() - h3 - two * u1 * h2;
        Self {
            x,
            y: r * (u1 * h2 - x) - s1 * h3,
            z: h * self.z * rhs.z,
            _config: PhantomData,
        }
    }
}

impl<C: SWCurveConfig> ConditionallySelectable for ProjectivePoint<C> {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        Self {
            x: PrimeFieldElement::conditional_select(&a.x, &b.x, choice),
            y: PrimeFieldElement::conditional_select(&a.y, &b.y, choice),
            z: PrimeFieldElement::conditional_select(&a.z, &b.z, choice),
            _config: PhantomData,
        }
    }
}
