use crate::U1024;
use crate::affine::{AffinePoint, SWCurveConfig};

use super::hash_message;

#[derive(Clone, Debug)]
pub struct SchnorrSignature<C: SWCurveConfig> {
    pub r_point: AffinePoint<C>,
    pub s: U1024,
}

impl<C: SWCurveConfig> SchnorrSignature<C> {
    fn challenge_hash(r: &AffinePoint<C>, message: &[u8], order: &U1024) -> U1024 {
        let r_x = r.x.to_bytes();
        let r_y = r.y.to_bytes();

        let mut combined = Vec::with_capacity(256 + message.len());
        combined.extend_from_slice(&r_x);
        combined.extend_from_slice(&r_y);
        combined.extend_from_slice(message);

        hash_message(&combined).mod_reduce(order)
    }

    /// 1. k = random in [1, n-1]
    /// 2. R = k * G
    /// 3. e = challenge_hash(R, message, n)
    /// 4. s = (k + e * private_key) mod n
    pub fn sign(private_key: &U1024, message: &[u8]) -> Self {
        let k = U1024::random_below(&C::ORDER);
        let r = C::generator().mul(&k);
        let e = Self::challenge_hash(&r, message, &C::ORDER);
        let s = e.mod_mul(private_key, &C::ORDER).mod_add(&k, &C::ORDER);
        Self { r_point: r, s }
    }

    /// 1. e = challenge_hash(R, message, n)
    /// 2. V1 = s * G
    /// 3. V2 = R + e * public_key
    /// 4. Valid if V1 == V2
    pub fn verify(&self, public_key: &AffinePoint<C>, message: &[u8]) -> bool {
        let e = Self::challenge_hash(&self.r_point, message, &C::ORDER);
        let v1 = C::generator().mul(&self.s);
        let v2 = self.r_point.add(&public_key.mul(&e));
        v1 == v2
    }
}
