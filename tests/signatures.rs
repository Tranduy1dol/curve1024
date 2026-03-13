//! Schnorr and ECDSA sign/verify round-trips.
use curve1024::{Curve1024Config, EcdsaSignature, KeyPair, SWCurveConfig, SchnorrSignature};

fn ecdsa_sign(key: &curve1024::U1024, msg: &[u8]) -> EcdsaSignature {
    let mut sig = EcdsaSignature::sign::<Curve1024Config>(key, msg);
    while sig.r.is_zero() || sig.s.is_zero() {
        sig = EcdsaSignature::sign::<Curve1024Config>(key, msg);
    }
    sig
}

#[test]
fn test_schnorr_sign_verify_roundtrip() {
    let kp = KeyPair::<Curve1024Config>::generate();
    let pk = Curve1024Config::generator().mul(&kp.private_key);
    let sig = SchnorrSignature::<Curve1024Config>::sign(&kp.private_key, b"hello schnorr");
    assert!(sig.verify(&pk, b"hello schnorr"));
}

#[test]
fn test_schnorr_rejects_wrong_message() {
    let kp = KeyPair::<Curve1024Config>::generate();
    let pk = Curve1024Config::generator().mul(&kp.private_key);
    let sig = SchnorrSignature::<Curve1024Config>::sign(&kp.private_key, b"correct");
    assert!(!sig.verify(&pk, b"tampered"));
}

#[test]
fn test_schnorr_rejects_wrong_key() {
    let kp = KeyPair::<Curve1024Config>::generate();
    let wrong_kp = KeyPair::<Curve1024Config>::generate();
    let wrong_pk = Curve1024Config::generator().mul(&wrong_kp.private_key);
    let sig = SchnorrSignature::<Curve1024Config>::sign(&kp.private_key, b"message");
    assert!(!sig.verify(&wrong_pk, b"message"));
}

#[test]
fn test_ecdsa_sign_verify_roundtrip() {
    let kp = KeyPair::<Curve1024Config>::generate();
    let pk = Curve1024Config::generator().mul(&kp.private_key);
    let sig = ecdsa_sign(&kp.private_key, b"hello ecdsa");
    assert!(sig.verify::<Curve1024Config>(&pk, b"hello ecdsa"));
}

#[test]
fn test_ecdsa_rejects_wrong_message() {
    let kp = KeyPair::<Curve1024Config>::generate();
    let pk = Curve1024Config::generator().mul(&kp.private_key);
    let sig = ecdsa_sign(&kp.private_key, b"correct");
    assert!(!sig.verify::<Curve1024Config>(&pk, b"tampered"));
}
