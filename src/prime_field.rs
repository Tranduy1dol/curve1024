use std::{
    marker::PhantomData,
    ops::{Add, Mul, Neg, Sub},
};

use subtle::ConditionallySelectable;

use crate::U1024;

pub trait PrimeFieldConfig {
    const MODULUS: U1024;
    const R2: U1024;
    const N_PRIME: U1024;
}

pub struct PrimeFieldElement<C: PrimeFieldConfig> {
    value: U1024,
    _config: PhantomData<C>,
}

impl<C: PrimeFieldConfig> PrimeFieldElement<C> {
    pub fn new(value: U1024) -> Self {
        let (lo, hi) = value.widening_mul(&C::R2);
        let value = Self::reduce(&lo, &hi);

        Self {
            value,
            _config: PhantomData,
        }
    }

    pub fn zero() -> Self {
        Self {
            value: U1024::ZERO,
            _config: PhantomData,
        }
    }

    pub fn one() -> Self {
        Self {
            value: U1024::ONE,
            _config: PhantomData,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.value == U1024::ZERO
    }

    fn reduce(lo: &U1024, hi: &U1024) -> U1024 {
        let (m, _) = lo.widening_mul(&C::N_PRIME);

        let (mn_lo, mn_hi) = m.widening_mul(&C::MODULUS);

        let (_, carry_lo) = lo.carrying_add(&mn_lo);

        let (res_hi, carry_hi) = hi.carrying_add(&mn_hi);

        let (t, carry_final) = if carry_lo {
            res_hi.carrying_add(&U1024::ONE)
        } else {
            (res_hi, false)
        };

        if carry_hi || carry_final {
            let (sub_res, _) = t.borrowing_sub(&C::MODULUS);
            return sub_res;
        }

        let (sub_res, borrow) = t.borrowing_sub(&C::MODULUS);
        if !borrow { sub_res } else { t }
    }
}

impl<C: PrimeFieldConfig> Add for PrimeFieldElement<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let (sum, carry) = self.value.carrying_add(&rhs.value);
        let (sub_res, borrow) = sum.borrowing_sub(&C::MODULUS);
        let use_sub = (carry || !borrow) as u8;

        Self {
            value: U1024::conditional_select(&sub_res, &sum, use_sub.into()),
            _config: PhantomData,
        }
    }
}

impl<C: PrimeFieldConfig> Sub for PrimeFieldElement<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let (diff, borrow) = self.value.borrowing_sub(&rhs.value);
        let (corr, _) = diff.carrying_add(&C::MODULUS);

        Self {
            value: U1024::conditional_select(&corr, &diff, (borrow as u8).into()),
            _config: PhantomData,
        }
    }
}

impl<C: PrimeFieldConfig> Mul for PrimeFieldElement<C> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let (lo, hi) = self.value.widening_mul(&rhs.value);
        let res = Self::reduce(&lo, &hi);

        Self {
            value: res,
            _config: PhantomData,
        }
    }
}

impl<C: PrimeFieldConfig> Neg for PrimeFieldElement<C> {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_zero() {
            self
        } else {
            let (value, _) = C::MODULUS.borrowing_sub(&self.value);
            Self {
                value,
                _config: PhantomData,
            }
        }
    }
}
