use crate::U1024;
use crate::affine::{AffinePoint, SWCurveConfig};

use super::hash_message;

#[derive(Clone, Debug)]
pub struct EcdsaSignature {
    pub r: U1024,
    pub s: U1024,
}

impl EcdsaSignature {
    /// 1. k = random in [1, n-1]
    /// 2. R = k * G
    /// 3. r = R.x mod n (retry if r == 0)
    /// 4. e = hash_message(message)
    /// 5. s = k^{-1} * (e + r * private_key) mod n (retry if s == 0)
    pub fn sign<C: SWCurveConfig>(private_key: &U1024, message: &[u8]) -> Self {
        let e = hash_message(message);

        let k = U1024::random_below(&C::ORDER);
        let r_point = C::generator().mul(&k);

        let r = r_point.x.to_u1024().mod_reduce(&C::ORDER);

        let k_inv = k.mod_inverse(&C::ORDER);
        let s =
            (e.mod_add(&r.mod_mul(private_key, &C::ORDER), &C::ORDER)).mod_mul(&k_inv, &C::ORDER);
        Self { r, s }
    }

    /// 1. Check r, s in [1, n-1]
    /// 2. e = hash_message(message)
    /// 3. w = s^{-1} mod n
    /// 4. u1 = (e * w) mod n, u2 = (r * w) mod n
    /// 5. P = u1 * G + u2 * public_key
    /// 6. Valid if P is not infinity and P.x mod n == r
    pub fn verify<C: SWCurveConfig>(&self, public_key: &AffinePoint<C>, message: &[u8]) -> bool {
        if self.r.is_zero() || self.s.is_zero() {
            return false;
        }

        if self.r >= C::ORDER || self.s >= C::ORDER {
            return false;
        }

        let e = hash_message(message);

        let w = self.s.mod_inverse(&C::ORDER);
        let u1 = e.mod_mul(&w, &C::ORDER);
        let u2 = self.r.mod_mul(&w, &C::ORDER);

        let u1_g = C::generator().mul(&u1);
        let u2_pub = public_key.mul(&u2);
        let p = u1_g.add(&u2_pub);

        if p.is_infinite {
            return false;
        }

        p.x.to_u1024().mod_reduce(&C::ORDER) == self.r 
    }
}
