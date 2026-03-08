use crate::{AffinePoint, PrimeFieldElement, SWCurveConfig, U1024};

use super::hash_message;

#[derive(Clone, Debug)]
pub struct SchnorrSignature<C: SWCurveConfig> {
    pub r_point: AffinePoint<C>,
    pub s: U1024,
}

impl<C: SWCurveConfig> SchnorrSignature<C> {
    fn challenge_hash(r: &AffinePoint<C>, message: &[u8]) -> PrimeFieldElement<C::ScalarField> {
        let r_x = r.x.to_bytes();
        let r_y = r.y.to_bytes();

        let mut combined = Vec::with_capacity(256 + message.len());
        combined.extend_from_slice(&r_x);
        combined.extend_from_slice(&r_y);
        combined.extend_from_slice(message);

        PrimeFieldElement::new(hash_message(&combined))
    }

    /// 1. k = random in [1, n-1]
    /// 2. R = k * G
    /// 3. e = challenge_hash(R, message, n)
    /// 4. s = (k + e * private_key) mod n
    pub fn sign(private_key: &U1024, message: &[u8]) -> Self {
        let k = U1024::rand(&C::ORDER);
        let r_point = C::generator().mul(&k);

        let k = PrimeFieldElement::<C::ScalarField>::new(k);
        let e = Self::challenge_hash(&r_point, message);
        let private_key = PrimeFieldElement::<C::ScalarField>::new(*private_key);
        let s = k + e * private_key;

        Self {
            r_point,
            s: s.value,
        }
    }

    /// 1. e = challenge_hash(R, message, n)
    /// 2. V1 = s * G
    /// 3. V2 = R + e * public_key
    /// 4. Valid if V1 == V2
    pub fn verify(&self, public_key: &AffinePoint<C>, message: &[u8]) -> bool {
        if self.s >= C::ORDER {
            return false;
        }

        let e = Self::challenge_hash(&self.r_point, message);
        let v1 = C::generator().mul(&self.s);
        let v2 = self.r_point.add(&public_key.mul(&e.value));

        v1 == v2
    }
}
