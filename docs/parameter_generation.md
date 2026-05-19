# Parameter Generation

This document describes the methodology used to generate the elliptic curve parameters for Curve1024, including the Cocks-Pinch construction, the NTT-friendliness constraints, and the readability-aware prime selection.

## Construction Method

Curve1024 uses the **Cocks-Pinch method** to construct a pairing-friendly elliptic curve with embedding degree $k = 18$ over a 1024-bit prime field. The construction proceeds in two phases:

1. **`cocks_pinch.rs`** — Finds the field primes $(p, r)$ and CM parameters $(t, y)$
2. **`complex_multiplication.rs`** — Derives the curve equation $(a, b)$ and a generator point $G$

### Phase 1: Prime Generation (Cocks-Pinch)

#### Step 1 — Scalar Field Prime $r$

A random 86-bit seed $T$ is chosen such that:

$$r = \Phi_{18}(T) = T^6 - T^3 + 1$$

is a 512-bit prime, where $\Phi_{18}$ is the 18th cyclotomic polynomial. This guarantees that $r \mid \Phi_{18}(T)$, which is a necessary condition for the curve to have embedding degree $k = 18$.

To ensure NTT-friendliness, $T$ is constrained to be a multiple of $2^{\lceil s/3 \rceil}$, where $s$ is the desired minimum two-adicity. Since $r - 1 = T^3(T^3 - 1)$, if $T \equiv 0 \pmod{2^{12}}$, then $2^{36} \mid (r - 1)$.

#### Step 2 — Base Field Prime $p$

Given a valid $r$, the algorithm computes:

1. A square root $\sqrt{-D} \pmod{r}$ (with $D = 3$)
2. Initial values $t_0 \equiv T^i + 1 \pmod{r}$ and $y_0 \equiv (t_0 - 2) \cdot (\sqrt{-D})^{-1} \pmod{r}$
3. Lifted values $t = t_0 + h_t \cdot r$ and $y = y_0 + h_y \cdot r$

The base field prime is then:

$$p = \frac{t^2 + D \cdot y^2}{4}$$

The algorithm searches over lift parameters $(h_t, h_y)$ until $p$ is a 1024-bit prime with the required two-adicity.

### Phase 2: Curve Equation (Complex Multiplication)

Given the CM discriminant $D = 3$, the $j$-invariant is $j = 0$ (from the Hilbert class polynomial $H_{-3}(x) = x$). This yields an initial curve family $y^2 = x^3 + b$.

The correct twist $b$ is identified by checking which of the six possible CM twist orders is divisible by $r$, and a generator $G$ of prime order $r$ is found by cofactor multiplication.

## NTT-Friendly Prime Selection

Both primes are constructed to satisfy:

$$p \equiv 1 \pmod{2^{36}}, \quad r \equiv 1 \pmod{2^{36}}$$

This two-adicity of 36 means both $\mathbb{F}_p$ and $\mathbb{F}_r$ contain primitive $2^{36}$-th roots of unity, enabling Number Theoretic Transforms (NTT) of length up to $2^{36} \approx 6.8 \times 10^{10}$.

The NTT-friendliness arises naturally from the construction:
- **For $r$**: Since $T \equiv 0 \pmod{2^{12}}$, we get $r - 1 = T^3(T^3 - 1)$, and $T^3$ contributes a factor of $2^{36}$.
- **For $p$**: The lift search filters candidates to ensure $2^{36} \mid (p - 1)$.

## Readability-Aware Search

### Motivation

Cryptographic primes are typically opaque hex blobs, making auditing and testing difficult. To improve auditability, the parameter generator includes a **readability scoring system** that evaluates candidate primes across multiple structural criteria and selects the most human-readable valid parameters.

### Scoring Criteria

Each candidate prime is scored on eight independent criteria:

| Criterion | What it rewards | Weight |
|---|---|---|
| **Binary sparseness** | Low Hamming weight (fewer set bits) | `(512 − hw) × 2` |
| **Two-adicity** | Higher power of 2 dividing $n - 1$ | `adicity × 3` |
| **Zero limbs** | Entire 64-bit blocks being zero | `count × 60` |
| **Full limbs** | Limbs equal to `0xFFFF...FFFF` (near power of 2) | `count × 50` |
| **Longest zero run** | Contiguous stretch of zero limbs | `run × 40` |
| **Hex zero density** | Many `0` digits in hex representation | `(zeros − 16) × 2` |
| **Simple top limb** | Few set bits in the leading 64-bit limb | `(32 − hw) × 3` |
| **Repeating limbs** | Adjacent identical limb values | `pairs × 25` |

The combined score weights the base field prime $p$ at 1.5× the scalar field prime $r$, since $p$ appears more frequently in implementations.

### Usage

```bash
# Fast mode — first valid prime (original behavior, ~10 minutes)
cargo run --release --example cocks_pinch

# Readable mode — best of 50 successes, up to 8 hours
cargo run --release --example cocks_pinch -- --readable

# Custom: best of 100 successes, up to 12 hours
cargo run --release --example cocks_pinch -- --readable --target 100 --hours 12
```

In readable mode, the generator:
1. Runs the normal Cocks-Pinch search
2. For each valid $r$, searches the **full lift range** for the best-scoring $p$ (instead of taking the first valid one)
3. Tracks the global best $(p, r)$ pair across all successful attempts
4. Writes the best result to `config/curve1024.toml` at the end

## Current Parameters

**Seed:** $T = \texttt{0x26704d2ace0facdd539} \cdot 2^{12}$

**Scalar field order** ($r = T^6 - T^3 + 1$, 512 bits):
```
r = 0x c042c338 9e72044d 9ec8078a ea7bc954
      ebf209d9 18f0ef27 425259e2 47711927
      97463bf2 3832d0de 5e46e47f c61cab50
      4bcfd36f c2a945c5 1055d970 00000001
```

**Base field prime** ($p = (t^2 + 3y^2)/4$, 1024 bits):
```
p = 0x c0859da8 3a500988 10eea35c 61fea054
      c8669025 99220d30 c93f3826 5ed6b830
      5785b8a0 73e482e8 685cfeeb 388dc45e
      bca16862 45b326fe 824ca74a c0844160
      f4e28980 6084e50e 6c02f32c b94d786b
      f80fb080 c07d6a07 73f2b459 1eb70df3
      9367aa32 ebb42449 e0a648f6 1a20e876
      f971f3bf 976a7d6a dd20d970 00000001
```

**Curve equation:** $E: y^2 = x^3 + 41$

### Readability Metrics

| Metric | $p$ | $r$ |
|---|---|---|
| Hamming weight | 464/1024 | 235/512 |
| Two-adicity | 36 | 36 |
| Zero 64-bit limbs | 0/16 | 8/16 |
| Hex zero digits | 33/256 | 16/128 |

Both primes share the tail pattern `...d97·2³⁶ + 1`, a consequence of $T \equiv 0 \pmod{2^{12}}$ propagating through the construction.

## Comparison with Standard Curve Presentations

The Cocks-Pinch method provides a polynomial structure for $r$ but not for $p$. This is a fundamental difference from parametric families like BN or BLS:

| Family | $r$ | $p$ | $\rho = \log p / \log r$ |
|---|---|---|---|
| **BN** | $r(x) = 36x^4 + 36x^3 + 18x^2 + 6x + 1$ | $p(x) = 36x^4 + 36x^3 + 24x^2 + 6x + 1$ | 1.0 |
| **BLS12** | $r(x) = x^4 - x^2 + 1$ | $p(x) = (x - 1)^2(x^4 - x^2 + 1)/3 + x$ | 1.5 |
| **KSS-18** | $r(x) = x^6 + 37x^3 + 343$ | $p(x) = (x^8 + 5x^7 + \cdots)/21$ | 1.33 |
| **Curve1024** (Cocks-Pinch) | $r = \Phi_{18}(T)$ | $p = (t^2 + 3y^2)/4$ | 2.0 |

The higher $\rho$-value of 2.0 means Curve1024 uses a larger base field relative to the scalar field compared to parametric families. This is the standard tradeoff of the Cocks-Pinch method: greater flexibility in parameter selection at the cost of a larger $p$.

## References

- Cocks, C. and Pinch, R. (2001). "ID-based cryptosystems based on the Weil pairing."
- Freeman, D., Scott, M., and Teske, E. (2010). "A Taxonomy of Pairing-Friendly Elliptic Curves." *Journal of Cryptology*, 23(2), 224–280.
- Kachisa, E., Schaefer, E., and Scott, M. (2008). "Constructing Brezing-Weng Pairing-Friendly Elliptic Curves Using Elements in the Cyclotomic Field." *Pairing-Based Cryptography*, LNCS 5209.
