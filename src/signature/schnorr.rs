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
            s: s.to_u1024(),
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
        let v2 = self.r_point.add(&public_key.mul(&e.to_u1024()));

        v1 == v2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime_field::{PrimeFieldConfig, PrimeFieldElement};

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestBase;

    // p = 17. R = 1 mod 17.
    impl PrimeFieldConfig for TestBase {
        const MODULUS: U1024 = U1024([17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const R2: U1024 = U1024::ONE;
        const N_PRIME: U1024 = U1024([0x0f0f0f0f0f0f0f0f; 16]);
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestScalar;

    // order = 19
    // R^2 mod 19 = (2^1024)^2 mod 19 = 6
    // N_PRIME = -19^-1 mod 2^1024
    impl PrimeFieldConfig for TestScalar {
        const MODULUS: U1024 = U1024([19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const R2: U1024 = U1024([6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        const N_PRIME: U1024 = U1024([
            0x79435e50d79435e5,
            0x435e50d79435e50d,
            0x5e50d79435e50d79,
            0x50d79435e50d7943,
            0xd79435e50d79435e,
            0x9435e50d79435e50,
            0x35e50d79435e50d7,
            0xe50d79435e50d794,
            0x0d79435e50d79435,
            0x79435e50d79435e5,
            0x435e50d79435e50d,
            0x5e50d79435e50d79,
            0x50d79435e50d7943,
            0xd79435e50d79435e,
            0x9435e50d79435e50,
            0x35e50d79435e50d7,
        ]);
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct TestCurve;

    // y^2 = x^3 + 2x + 2 mod 17
    // ORDER = 19
    // Generator = (5, 1)
    impl SWCurveConfig for TestCurve {
        type BaseField = TestBase;
        type ScalarField = TestScalar;

        const COEFF_A: PrimeFieldElement<Self::BaseField> =
            PrimeFieldElement::from_montgomery(U1024([
                2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]));
        const COEFF_B: PrimeFieldElement<Self::BaseField> =
            PrimeFieldElement::from_montgomery(U1024([
                2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]));
        const ORDER: U1024 = U1024([19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        fn generator() -> AffinePoint<Self> {
            AffinePoint::new(
                PrimeFieldElement::from_montgomery(U1024([
                    5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ])),
                PrimeFieldElement::from_montgomery(U1024([
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ])),
            )
        }
    }

    #[test]
    fn test_valid_schnorr_signature() {
        let private_key = U1024::from_u64(10);
        let public_key = TestCurve::generator().mul(&private_key);

        let message = b"Hello, Schnorr!";

        let sig = SchnorrSignature::sign(&private_key, message);
        assert!(
            sig.verify(&public_key, message),
            "Signature failed verification"
        );

        // Verify tampering breaks validation
        let mut bad_sig = sig.clone();
        bad_sig.s = bad_sig.s.carrying_add(&U1024::ONE).0;
        assert!(
            !bad_sig.verify(&public_key, message),
            "Tampered signature should fail"
        );
    }
}
