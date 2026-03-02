use std::fmt;

use subtle::ConditionallySelectable;

const LIMBS: usize = 16;

#[repr(align(64))]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct U1024(pub [u64; LIMBS]);

impl U1024 {
    pub const NUM_LIMBS: usize = LIMBS;

    pub const ZERO: Self = Self([0; LIMBS]);

    pub const ONE: Self = Self([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    pub fn from_hex(hex: &str) -> Self {
        let hex = hex.trim_start_matches("0x");
        assert!(hex.len() <= 256, "Hex string too long for U1024");

        let mut res = Self::ZERO;
        let mut limb_idx = 0;
        let mut char_idx = hex.len();

        while char_idx > 0 {
            let start = char_idx.saturating_sub(16);
            let chunk = &hex[start..char_idx];

            let val = u64::from_str_radix(chunk, 16).expect("Invalid hex character");

            if limb_idx < LIMBS {
                res.0[limb_idx] = val;
            }
            limb_idx += 1;
            char_idx = start;
        }
        res
    }

    pub fn from_u64(v: u64) -> Self {
        let mut arr = [0; LIMBS];
        arr[0] = v;
        U1024(arr)
    }

    pub fn widening_mul(&self, rhs: &Self) -> (Self, Self) {
        let mut res = [0u64; LIMBS * 2];
        let mut i = 0;
        while i < LIMBS {
            let mut carry = 0u64;
            let mut j = 0;
            while j < LIMBS {
                let k = i + j;
                let val = res[k] as u128 + (self.0[i] as u128 * rhs.0[j] as u128) + carry as u128;
                res[k] = val as u64;
                carry = (val >> 64) as u64;
                j += 1;
            }
            let mut k = i + LIMBS;
            while carry > 0 && k < LIMBS * 2 {
                let val = res[k] as u128 + carry as u128;
                res[k] = val as u64;
                carry = (val >> 64) as u64;
                k += 1;
            }
            i += 1;
        }

        let mut low = [0u64; LIMBS];
        let mut high = [0u64; LIMBS];

        let mut k = 0;
        while k < LIMBS {
            low[k] = res[k];
            high[k] = res[k + LIMBS];
            k += 1;
        }

        (U1024(low), U1024(high))
    }

    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        if *divisor == Self::ZERO {
            panic!("Division by zero");
        }

        if *self < *divisor {
            return (Self::ZERO, *self);
        }

        if *self == *divisor {
            return (Self::ONE, Self::ZERO);
        }

        let mut quotient = Self::ZERO;
        let mut remainder = Self::ZERO;

        for i in (0..1024).rev() {
            remainder = remainder.shl(1);

            let bit = self.bit(i) as u8;
            let rem_plus_one = remainder.carrying_add(&Self::ONE).0;
            remainder = Self::conditional_select(&rem_plus_one, &remainder, bit.into());

            let (sub_result, borrow) = remainder.borrowing_sub(divisor);
            let should_subtract = (!borrow) as u8;

            remainder = Self::conditional_select(&sub_result, &remainder, should_subtract.into());

            let quotient_with_bit = quotient.with_bit(i);
            quotient =
                Self::conditional_select(&quotient_with_bit, &quotient, should_subtract.into());
        }

        (quotient, remainder)
    }

    pub fn borrowing_sub(&self, rhs: &Self) -> (Self, bool) {
        let mut ret = [0u64; LIMBS];
        let mut borrow = 0u64;
        let mut i = 0;

        while i < LIMBS {
            let (diff1, b1) = self.0[i].overflowing_sub(rhs.0[i]);
            let (diff2, b2) = diff1.overflowing_sub(borrow);

            ret[i] = diff2;
            borrow = (b1 as u64) + (b2 as u64);
            i += 1;
        }
        (U1024(ret), borrow != 0)
    }

    pub fn carrying_add(&self, rhs: &Self) -> (Self, bool) {
        let mut ret = U1024::ZERO;
        let mut carry = 0u64;

        for i in 0..LIMBS {
            let (sum1, c1) = self.0[i].overflowing_add(rhs.0[i]);
            let (sum2, c2) = sum1.overflowing_add(carry);
            ret.0[i] = sum2;
            carry = (c1 as u64) + (c2 as u64);
        }
        (ret, carry != 0)
    }

    pub fn shr(&self, n: usize) -> Self {
        if n >= 1024 {
            return Self::ZERO;
        }

        if n == 0 {
            return *self;
        }

        let limb_shift = n / 64;
        let bit_shift = n % 64;

        let mut result = [0u64; LIMBS];

        if bit_shift == 0 {
            result[limb_shift..LIMBS].copy_from_slice(&self.0[..(LIMBS - limb_shift)]);
        } else {
            for (i, result_limb) in result.iter_mut().enumerate().skip(limb_shift) {
                let src_idx = i - limb_shift;
                *result_limb = self.0[src_idx] >> bit_shift;
                if src_idx > 0 {
                    *result_limb |= self.0[src_idx - 1] << (64 - bit_shift);
                }
            }
        }

        Self(result)
    }

    pub fn shl(&self, n: usize) -> Self {
        if n >= 1024 {
            return Self::ZERO;
        }

        if n == 0 {
            return *self;
        }

        let limb_shift = n / 64;
        let bit_shift = n % 64;

        let mut result = [0u64; LIMBS];

        if bit_shift == 0 {
            result[limb_shift..LIMBS].copy_from_slice(&self.0[..(LIMBS - limb_shift)]);
        } else {
            for (i, result_limb) in result.iter_mut().enumerate().skip(limb_shift) {
                let src_idx = i - limb_shift;
                *result_limb = self.0[src_idx] << bit_shift;
                if src_idx > 0 {
                    *result_limb |= self.0[src_idx - 1] >> (64 - bit_shift);
                }
            }
        }

        Self(result)
    }

    fn with_bit(&self, index: usize) -> Self {
        if index >= 1024 {
            return *self;
        }

        let limb_idx = index / 64;
        let bit_idx = index % 64;

        let mut result = self.0;
        result[limb_idx] |= 1u64 << bit_idx;

        Self(result)
    }

    fn bit(&self, index: usize) -> bool {
        if index >= 1024 {
            return false;
        }

        let limb_idx = index / 64;
        let bit_idx = index % 64;

        (self.0[limb_idx] >> bit_idx) & 1 == 1
    }
}

impl fmt::Display for U1024 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0x")?;
        for limb in self.0.iter().rev() {
            write!(f, "{:016x}", limb)?;
        }
        Ok(())
    }
}

impl From<u64> for U1024 {
    fn from(value: u64) -> Self {
        U1024::from_u64(value)
    }
}

impl ConditionallySelectable for U1024 {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        let mut res = U1024([0; LIMBS]);
        for i in 0..LIMBS {
            let a_val = core::hint::black_box(a.0[i]);
            let b_val = core::hint::black_box(b.0[i]);

            res.0[i] = u64::conditional_select(&a_val, &b_val, choice);
        }
        res
    }
}
