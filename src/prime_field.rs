use std::{
    marker::PhantomData,
    ops::{Add, Mul, Neg, Sub},
};

use subtle::ConditionallySelectable;

use crate::U1024;

pub trait PrimeFieldConfig: 'static + Copy + Clone + Eq + PartialEq {
    const MODULUS: U1024;
    const R2: U1024;
    const N_PRIME: U1024;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrimeFieldElement<C: PrimeFieldConfig> {
    value: U1024,
    _config: PhantomData<C>,
}

impl<C: PrimeFieldConfig> PrimeFieldElement<C> {
    pub const ZERO: Self = Self {
        value: U1024::ZERO,
        _config: PhantomData,
    };

    pub const fn from_montgomery(value: U1024) -> Self {
        Self {
            value,
            _config: PhantomData,
        }
    }

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
        Self::new(U1024::ONE)
    }

    pub fn is_zero(&self) -> bool {
        self.value == U1024::ZERO
    }

    pub fn inv(&self) -> Self {
        let two = U1024::from(2);
        let (p_minus_2, _) = C::MODULUS.borrowing_sub(&two);
        self.pow(p_minus_2)
    }

    pub fn pow(&self, exp: U1024) -> Self {
        let mut res = Self::one();
        let mut base = *self;

        for i in 0..16 {
            let mut limb = exp.0[i];
            for _ in 0..64 {
                let bit = ((limb & 1) == 1) as u8;
                let product = res * base;

                res = Self::conditional_select(&res, &product, bit.into());
                base = base.square();

                limb >>= 1;
            }
        }
        res
    }

    pub fn square(&self) -> Self {
        *self * *self
    }

    pub fn to_u1024(&self) -> U1024 {
        Self::reduce(&self.value, &U1024::ZERO)
    }

    pub fn to_bytes(&self) -> [u8; 128] {
        let canonical = self.to_u1024();
        canonical.to_be_bytes()
    }

    pub fn from_bytes(bytes: &[u8; 128]) -> Self {
        let value = U1024::from_be_bytes(bytes);
        Self::new(value)
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
            value: U1024::conditional_select(&sum, &sub_res, use_sub.into()),
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
            value: U1024::conditional_select(&diff, &corr, (borrow as u8).into()),
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

impl<C: PrimeFieldConfig> ConditionallySelectable for PrimeFieldElement<C> {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        Self {
            value: U1024::conditional_select(&a.value, &b.value, choice),
            _config: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestConfig;

    impl PrimeFieldConfig for TestConfig {
        const MODULUS: U1024 = U1024([17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const R2: U1024 = U1024::ONE;
        const N_PRIME: U1024 = U1024([0x0f0f0f0f0f0f0f0f; 16]);
    }

    type F = PrimeFieldElement<TestConfig>;

    #[test]
    fn test_basics() {
        let a = F::new(U1024::from_u64(5));
        assert_eq!(a.to_u1024(), U1024::from_u64(5));

        let b = F::new(U1024::from_u64(12));
        assert_eq!(b.to_u1024(), U1024::from_u64(12));

        assert!(a != b);
    }

    #[test]
    fn test_arithmetic() {
        let a = F::new(U1024::from_u64(10));
        let b = F::new(U1024::from_u64(12));

        // 10 + 12 = 22 = 5 mod 17
        let add = a + b;
        assert_eq!(add.to_u1024(), U1024::from_u64(5));

        // 10 - 12 = -2 = 15 mod 17
        let sub = a - b;
        assert_eq!(sub.to_u1024(), U1024::from_u64(15));

        // 10 * 12 = 120 = 1 mod 17
        let mul = a * b;
        assert_eq!(mul.to_u1024(), U1024::from_u64(1));

        // 10^-1 = 12 mod 17
        assert_eq!(a.inv(), b);
        assert_eq!(b.inv(), a);

        let neg_a = -a;
        assert_eq!(neg_a.to_u1024(), U1024::from_u64(7));

        // Square: 10^2 = 100 = 15 mod 17
        let sq = a.square();
        assert_eq!(sq.to_u1024(), U1024::from_u64(15));
    }

    #[test]
    fn test_pow() {
        let a = F::new(U1024::from_u64(2));
        // 2^4 = 16 mod 17
        assert_eq!(a.pow(U1024::from_u64(4)).to_u1024(), U1024::from_u64(16));
        // Fermat's little theorem: 2^16 = 1 mod 17
        assert_eq!(a.pow(U1024::from_u64(16)).to_u1024(), U1024::from_u64(1));
    }

    #[test]
    fn test_zero_one() {
        let zero = F::zero();
        let one = F::one();

        let a = F::new(U1024::from_u64(7));

        assert_eq!(a * zero, zero);
        assert_eq!(a * one, a);
        assert_eq!(a + zero, a);

        assert!(zero.is_zero());
        assert!(!one.is_zero());

        assert_eq!(F::ZERO, F::zero());
    }

    #[test]
    fn test_from_montgomery() {
        let a = F::new(U1024::from_u64(10));
        // Since R = 1 mod 17, the internal Montgomery form inherently equals the standard literal value.
        let raw = F::from_montgomery(U1024::from_u64(10));
        assert_eq!(a, raw);
    }

    #[test]
    fn test_bytes() {
        let a = F::new(U1024::from_u64(14));
        let bytes = a.to_bytes();
        let b = F::from_bytes(&bytes);
        assert_eq!(a, b);
    }

    #[test]
    fn test_conditionally_selectable() {
        let a = F::new(U1024::from_u64(5));
        let b = F::new(U1024::from_u64(10));

        let choice_0 = subtle::Choice::from(0);
        let choice_1 = subtle::Choice::from(1);

        assert_eq!(F::conditional_select(&a, &b, choice_0), a);
        assert_eq!(F::conditional_select(&a, &b, choice_1), b);
    }
}
