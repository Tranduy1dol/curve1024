use crate::{AffinePoint, PrimeFieldElement, SWCurveConfig, U1024};

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
        let k = U1024::rand(&C::ORDER);
        let r_point = C::generator().mul(&k);

        let e = PrimeFieldElement::<C::ScalarField>::new(hash_message(message));
        let r = PrimeFieldElement::<C::ScalarField>::new(r_point.x.to_u1024());
        let k_inv = PrimeFieldElement::<C::ScalarField>::new(k).inv();
        let private_key = PrimeFieldElement::<C::ScalarField>::new(*private_key);

        let s = k_inv * (e + r * private_key);
        Self {
            r: r.value,
            s: s.value,
        }
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

        let e = PrimeFieldElement::<C::ScalarField>::new(hash_message(message));
        let w = PrimeFieldElement::<C::ScalarField>::new(self.s).inv();
        let r = PrimeFieldElement::<C::ScalarField>::new(self.r);

        let u1 = e * w;
        let u2 = r * w;

        let u1_g = C::generator().mul(&u1.value);
        let u2_pub = public_key.mul(&u2.value);
        let p = u1_g.add(&u2_pub);

        if p.is_infinite {
            return false;
        }

        PrimeFieldElement::<C::ScalarField>::new(p.x.value).value == self.r
    }
}
