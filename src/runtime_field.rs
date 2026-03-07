use crate::U1024;
use crate::u1024::LIMBS;

pub struct RuntimeField {
    pub modulus: U1024,
    n_prime: U1024,
    r2: U1024,
}

impl RuntimeField {
    pub fn new(modulus: U1024) -> Self {
        let n_prime = Self::compute_n_prime(&modulus);
        let r2 = Self::compute_r2(&modulus);
        Self {
            modulus,
            n_prime,
            r2,
        }
    }

    pub fn reduce(&self, a: &U1024) -> U1024 {
        a.div_rem(&self.modulus).1
    }

    pub fn add(&self, a: &U1024, b: &U1024) -> U1024 {
        let (sum, carry) = a.carrying_add(b);
        let (sub_res, borrow) = sum.borrowing_sub(&self.modulus);
        let use_sub = carry || !borrow;
        if use_sub { sub_res } else { sum }
    }

    pub fn mul(&self, a: &U1024, b: &U1024) -> U1024 {
        let a_mont = self.encode_mont(a);
        let b_mont = self.encode_mont(b);
        let r_mont = self.mont_mul(&a_mont, &b_mont);
        self.decode_mont(&r_mont)
    }

    pub fn pow(&self, base: &U1024, exp: &U1024) -> U1024 {
        if self.modulus == U1024::ONE {
            return U1024::ZERO;
        }

        let mut result = self.encode_mont(&U1024::ONE);
        let mut b = self.encode_mont(&self.reduce(base));

        let top_limb = (0..LIMBS)
            .rev()
            .find(|&i| exp.0[i] != 0)
            .map_or(0, |i| i + 1);

        for i in 0..top_limb {
            let mut limb = exp.0[i];
            for _ in 0..64 {
                if limb & 1 == 1 {
                    result = self.mont_mul(&result, &b);
                }
                b = self.mont_mul(&b, &b);
                limb >>= 1;
            }
        }

        self.decode_mont(&result)
    }

    pub fn inv(&self, a: &U1024) -> U1024 {
        assert!(!a.is_zero());
        let two = U1024::from_u64(2);
        let (exp, _) = self.modulus.borrowing_sub(&two);
        self.pow(a, &exp)
    }

    pub fn encode_mont(&self, a: &U1024) -> U1024 {
        self.mont_mul(a, &self.r2)
    }

    pub fn decode_mont(&self, a: &U1024) -> U1024 {
        Self::mont_reduce(a, &U1024::ZERO, &self.modulus, &self.n_prime)
    }

    fn mont_mul(&self, a: &U1024, b: &U1024) -> U1024 {
        let (lo, hi) = a.widening_mul(b);
        Self::mont_reduce(&lo, &hi, &self.modulus, &self.n_prime)
    }

    fn mont_reduce(lo: &U1024, hi: &U1024, n: &U1024, n_prime: &U1024) -> U1024 {
        let (m, _) = lo.widening_mul(n_prime);
        let (mn_lo, mn_hi) = m.widening_mul(n);
        let (_, carry_lo) = lo.carrying_add(&mn_lo);
        let (res_hi, carry_hi) = hi.carrying_add(&mn_hi);
        let (t, carry_final) = if carry_lo {
            res_hi.carrying_add(&U1024::ONE)
        } else {
            (res_hi, false)
        };
        if carry_hi || carry_final {
            let (sub_res, _) = t.borrowing_sub(n);
            return sub_res;
        }
        let (sub_res, borrow) = t.borrowing_sub(n);
        if !borrow { sub_res } else { t }
    }

    fn compute_n_prime(n: &U1024) -> U1024 {
        let mut inv = U1024::ONE;
        let two = U1024::from_u64(2);
        for _ in 0..10 {
            let ni = n.widening_mul(&inv).0;
            let diff = two.borrowing_sub(&ni).0;
            inv = inv.widening_mul(&diff).0;
        }
        U1024::ZERO.borrowing_sub(&inv).0
    }

    fn compute_r2(n: &U1024) -> U1024 {
        let (neg_n, _) = U1024::ZERO.borrowing_sub(n);
        let r_mod_n = neg_n.div_rem(n).1;
        r_mod_n.mod_mul(&r_mod_n, n)
    }
}
