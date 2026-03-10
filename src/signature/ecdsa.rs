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
            r: r.to_u1024(),
            s: s.to_u1024(),
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

        let u1_g = C::generator().mul(&u1.to_u1024());
        let u2_pub = public_key.mul(&u2.to_u1024());
        let p = u1_g.add(&u2_pub);

        if p.is_infinite {
            return false;
        }

        let px_mod_n = PrimeFieldElement::<C::ScalarField>::new(p.x.to_u1024());
        px_mod_n.to_u1024() == self.r
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
    fn test_valid_ecdsa_signature() {
        let private_key = U1024::from_u64(10);
        let public_key = TestCurve::generator().mul(&private_key);

        let message = b"Hello, ECDSA!";

        let mut valid_signature_found = false;

        for _ in 0..10 {
            let sig = EcdsaSignature::sign::<TestCurve>(&private_key, message);

            if sig.r.is_zero() || sig.s.is_zero() {
                continue;
            }

            assert!(
                sig.verify::<TestCurve>(&public_key, message),
                "Signature failed verification"
            );

            // Verify tampering breaks validation
            let mut bad_sig = sig.clone();
            bad_sig.s = bad_sig.s.carrying_add(&U1024::ONE).0;
            assert!(
                !bad_sig.verify::<TestCurve>(&public_key, message),
                "Tampered signature should fail"
            );

            valid_signature_found = true;
            break;
        }

        assert!(
            valid_signature_found,
            "Failed to generate a valid non-zero signature in 10 attempts"
        );
    }
}
