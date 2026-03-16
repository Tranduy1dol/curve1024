# Curve1024

[![Rust](https://img.shields.io/badge/rust-1.70%2B-blue.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A complete, from-scratch **1024-bit Elliptic Curve Cryptography (ECC)** library written purely in Rust. Designed with cryptographic sovereignty and memory safety in mind, Curve1024 relies on zero external cryptographic dependencies. It provides a massive 256-bit security margin, targeting high-security, long-term, and offline signing scenarios.

## 🌟 Key Features

- **Pure-Rust 1024-bit BigNum Arithmetic (`U1024`)**: Custom, memory-safe large integer operations over 1024-bit fields, circumventing vulnerabilities inherent in legacy C/C++ bignum implementations.
- **Novel KSS18 Pairing-Friendly Curve**: Synthesized via an improved Cocks-Pinch algorithm. Features an embedding degree $k=18$, group order $r \approx 2^{512}$, and base field $p \approx 2^{1024}$.
- **Optimized Field Arithmetic**: Heavily utilizes Montgomery multiplication to bypass expensive large-integer divisions context, ensuring a viable performance-to-security trade-off.
- **Built-in Signature Schemes**: Complete, thoroughly tested implementations of **Schnorr** and **ECDSA** digital signature algorithms.
- **Robust Attack Resistance**:
  - **Pollard's rho / ECDLP**: Requires $O(2^{256})$ operations (absolute margin against classical computing limits).
  - **MOV Attack**: The extension field $\mathbb{F}_{p^{18}}$ is a massive 18,432 bits, making index-calculus DLP completely infeasible.
  - **Anomalous (SSSA) Attack**: Curve cardinality is rigorously checked ($\#E(\mathbb{F}_p) \neq p$).
  - **TNFS Attack Resistance**: Ensures NTT-friendly base and scalar fields.

## 📐 Mathematical Specifications

- **Curve Form**: Short Weierstrass ($y^2 \equiv x^3 + b \pmod p$)
- **Base Field ($p$)**: ~1024 bits
- **Scalar Field ($r$)**: ~512 bits
- **Embedding Degree ($k$)**: 18
- **Security Level**: 256-bit symmetric equivalent (exceeds NIST 128-bit post-2030 standards)

## 🏗️ Architecture

Curve1024 is structured into distinct, decoupled architectural layers:

1. **`U1024` (Big Integer Layer)**: Underlying memory-safe 1024-bit arithmetic operations.
2. **`PrimeField` (Field Layer)**: Modulo mathematics utilizing Montgomery space for optimized computations.
3. **`AffinePoint` (Curve Layer)**: Geometric operations, curve point verification, and scalar multiplication (`Double-and-Add`).
4. **Signatures**: Mathematical construction of Schnorr and ECDSA schemas.
5. **`curve1024-sig`**: Integrated Command Line Interface (CLI) for key generation, document signing, and verification.

## 🚀 Getting Started

### Prerequisites
- [Rust toolchain](https://rustup.rs/) (`cargo`, `rustc`)

### Build the Project
```bash
git clone https://github.com/Tranduy1dol/lumen-math.git
cd lumen-math
cargo build --release
```

### Run the Test Suite
Curve1024 is bundled with an extensive suite of 36 mathematically rigorous test cases, including base arithmetic edge cases and simulated attack analysis (e.g., MOV and Anomalous context checks).

```bash
cargo test
```

*Tip: To evaluate the computationally-heavy benchmarking, you can run `cargo bench`.*

## 🔒 Security Disclaimer

**Warning & Academic Transparency:**
While this library achieves high mathematical correctness and effectively neutralizes standard algorithmic attacks, the scalar multiplication algorithm currently leverages a *Double-and-Add* approach. Thus, it **may be vulnerable to timing side-channel attacks**, as execution time correlates with the scalar's Hamming weight. 

Future developments aim to replace this with a constant-time *Montgomery Ladder* and optimize through *Jacobian/Projective constraints*. **Curve1024 is currently best suited for academic demonstrations, secure protocol designs, and offline signing operations.**

## 📜 License

This project is open source and strictly licensed under the MIT License.
