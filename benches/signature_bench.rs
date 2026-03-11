use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use curve1024::{Curve1024Config, EcdsaSignature, KeyPair, SWCurveConfig, SchnorrSignature, U1024};
use std::{hint::black_box, time::Duration};

fn bench_keygen(c: &mut Criterion) {
    let mut group = c.benchmark_group("KeyGen");
    group
        .sample_size(10)
        .measurement_time(Duration::from_secs(30));
    group.bench_function("KeyPair generate", |b| {
        b.iter(|| KeyPair::<Curve1024Config>::generate())
    });
    group.finish();
}

fn bench_signatures(c: &mut Criterion) {
    let keypair = KeyPair::<Curve1024Config>::generate();
    let private_key = keypair.private_key;
    let public_key = keypair.public_key;
    let message = b"Hello, World - we are benchmarking curve signatures!";

    // Pre-generate a valid ECDSA signature (r, s both non-zero)
    let ecdsa_sig = {
        let mut sig = EcdsaSignature::sign::<Curve1024Config>(&private_key, message);
        while sig.r.is_zero() || sig.s.is_zero() {
            sig = EcdsaSignature::sign::<Curve1024Config>(&private_key, message);
        }
        sig
    };
    let schnorr_sig = SchnorrSignature::<Curve1024Config>::sign(&private_key, message);

    let mut group = c.benchmark_group("Signatures");
    group
        .sample_size(10)
        .measurement_time(Duration::from_secs(30));

    group.bench_function("ECDSA/sign", |b| {
        b.iter(|| {
            EcdsaSignature::sign::<Curve1024Config>(black_box(&private_key), black_box(message))
        })
    });

    group.bench_function("ECDSA/verify", |b| {
        b.iter(|| ecdsa_sig.verify::<Curve1024Config>(black_box(&public_key), black_box(message)))
    });

    group.bench_function("Schnorr/sign", |b| {
        b.iter(|| {
            SchnorrSignature::<Curve1024Config>::sign(black_box(&private_key), black_box(message))
        })
    });

    group.bench_function("Schnorr/verify", |b| {
        b.iter(|| schnorr_sig.verify(black_box(&public_key), black_box(message)))
    });

    group.finish();
}

fn bench_scalar_mul(c: &mut Criterion) {
    let g = Curve1024Config::generator();
    let scalars: Vec<U1024> = (1u64..=4).map(U1024::from).collect();

    let mut group = c.benchmark_group("ScalarMul");
    group
        .sample_size(10)
        .measurement_time(Duration::from_secs(30));
    for scalar in &scalars {
        group.bench_with_input(BenchmarkId::new("G * k", scalar), scalar, |b, s| {
            b.iter(|| g.mul(black_box(s)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_keygen, bench_signatures, bench_scalar_mul);
criterion_main!(benches);
