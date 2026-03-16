# Curve1024

[![Rust](https://img.shields.io/badge/rust-1.70%2B-blue.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduce

**Curve1024** is a complete, from-scratch **1024-bit Elliptic Curve Cryptography (ECC)** library written purely in Rust. Designed with cryptographic sovereignty and memory safety at its core, Curve1024 relies on zero external cryptographic dependencies. 

By utilizing a monstrous 1024-bit base field ($p \approx 2^{1024}$) and a 512-bit scalar field ($r \approx 2^{512}$), this library provides an absolute 256-bit symmetric security margin. It is specifically targeted at high-security, long-term, and offline signing environments where extreme attack resistance is prioritized.

## Key Features

### The Library
- **Pure-Rust 1024-bit BigNum Arithmetic (`U1024`)**: Custom, memory-safe large integer operations over 1024-bit fields, circumventing vulnerabilities inherent in legacy C/C++ implementations.
- **Novel KSS18 Pairing-Friendly Curve**: Synthesized via an improved Cocks-Pinch algorithm. Features an embedding degree $k=18$.
- **Optimized Field Arithmetic**: Heavily utilizes Montgomery multiplication to bypass expensive large-integer divisions, ensuring a viable performance-to-security trade-off.
- **Robust Attack Resistance**:
  - *Pollard's rho*: Requires $O(2^{256})$ operations.
  - *MOV Attack*: The extension field $\mathbb{F}_{p^{18}}$ is a massive 18,432 bits, making index-calculus DLP unfeasible.
  - *Anomalous (SSSA)*: Curve cardinality is rigorously mathematically checked ($\#E(\mathbb{F}_p) \neq p$).
  - *TNFS*: Ensures NTT-friendly base and scalar fields.

### The CLI Tool (`curve1024-sig`)
- **GPG-like Interface**: An intuitive command-line interface for managing keys and signing/verifying files.
- **Multiple Signature Schemes**: Offers drop-in support for both **Schnorr** and **ECDSA** digital signature algorithms out of the box.
- **File Signatures**: Directly signs files, appending the signature securely to the file with a specific magic footer identifying the schema.

---

## Usage

### As a Library

Add `curve1024` to your `Cargo.toml` dependencies.

```rust
use curve1024::{
    AffinePoint, Curve1024BaseField, Curve1024Config, EcdsaSignature, KeyPair, 
    PrimeFieldElement, SchnorrSignature, U1024,
};

fn main() {
    // 1. Generate a new 1024-bit keypair
    let keypair = KeyPair::<Curve1024Config>::generate();
    let message = b"Confidential system payload";

    // 2. Sign using Schnorr
    let schnorr_sig = SchnorrSignature::<Curve1024Config>::sign(&keypair.private_key, message);
    let schnorr_valid = schnorr_sig.verify(&keypair.public_key, message);
    assert!(schnorr_valid);

    // 3. Or sign using ECDSA
    let ecdsa_sig = EcdsaSignature::sign::<Curve1024Config>(&keypair.private_key, message);
    let ecdsa_valid = ecdsa_sig.verify::<Curve1024Config>(&keypair.public_key, message);
    assert!(ecdsa_valid);
}
```

### The CLI Tool (`curve1024-sig`)

You can install the CLI tool directly from the repository using `cargo install`.

```bash
cargo install --git https://github.com/Tranduy1dol/curve1024.git curve1024-sig
```

Alternatively, you can build it from source:
```bash
git clone https://github.com/Tranduy1dol/curve1024.git
cd curve1024
cargo build --release
```

**CLI Commands Overview:**

```text
A GPG-like tool for elliptic curve key management and digital signatures

Usage: curve1024-sig <COMMAND>

Commands:
  keygen      
  sign        
  verify      
  export-pub  
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

---

## 🔒 Security Disclaimer

*Warning & Academic Transparency:* While this library achieves high mathematical correctness and effectively neutralizes standard algorithmic attacks, the scalar multiplication algorithm currently leverages a *Double-and-Add* approach. Thus, it **may be vulnerable to timing side-channel attacks**. Curve1024 is currently best suited for academic demonstrations, secure protocol designs, and offline signing operations.

## 📜 License
This project is open-source and strictly licensed under the MIT License.
