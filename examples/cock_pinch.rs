use curve1024::{U1024, u1024::LIMBS};
use rand::Rng;

struct CurveParams {
    p: U1024, // base field prime
    r: U1024, // group order (prime)
    t: U1024, // trace of Frobenius
    y: U1024, // CM parameter
    k: u64,   // embedding degree
    d: U1024, // CM discriminant |D|
}

fn cyclotomic_eval(k: u64, t_val: &U1024) -> U1024 {
    assert!(k == 18, "Only k=18 is supported");

    let t2 = t_val.widening_mul(t_val).0;
    let t3 = t2.widening_mul(t_val).0;
    let t6 = t3.widening_mul(&t3).0;

    let (diff, _) = t6.borrowing_sub(&t3);
    diff.carrying_add(&U1024::ONE).0
}

/// Miller-Rabin primality test for U1024 values.
fn is_prime(n: &U1024) -> bool {
    if *n == U1024::ONE {
        return false;
    }

    if *n == U1024::from(2) || *n == U1024::from(3) {
        return true;
    }

    if n.0[0] & 1 == 0 {
        return false;
    }

    let mut d = n.borrowing_sub(&U1024::ONE).0;
    let mut s = 0u64;
    while d.0[0] & 1 == 0 {
        d = d.shr(1);
        s += 1;
    }

    for i in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let a = U1024::from(i);
        if a >= *n {
            continue;
        }

        let mut x = a.mod_pow(&d, n);
        if x == U1024::ONE || x == n.borrowing_sub(&U1024::ONE).0 {
            continue;
        }

        let mut composite = true;
        for _ in 0..s.saturating_sub(1) {
            x = x.mod_mul(&x, n);
            if x == n.borrowing_sub(&U1024::ONE).0 {
                composite = false;
                break;
            }
        }

        if composite {
            return false;
        }
    }
    true
}

/// Compute sqrt(-D) mod r, returns None if -D is not a QR mod r.
fn sqrt_mod(a: &U1024, p: &U1024) -> Option<U1024> {
    let half_p_minus_1 = p.borrowing_sub(&U1024::ONE).0.div_rem(&U1024::from(2)).0;
    if a.mod_pow(&half_p_minus_1, p) != U1024::ONE {
        return None;
    }

    let mut q = a.borrowing_sub(&U1024::ONE).0;
    let mut s = 0;
    while q.0[0] & 1 == 0 {
        q = q.shr(1);
        s += 1;
    }

    let mut z = U1024::from(2);
    while z.mod_pow(&half_p_minus_1, p) != a.borrowing_sub(&U1024::ONE).0 {
        z = z.carrying_add(&U1024::ONE).0;
    }

    let mut m = s;
    let mut c = z.mod_pow(&q, p);
    let mut t = a.mod_pow(&q, p);
    let half_q_plus_1 = q.carrying_add(&U1024::ONE).0.div_rem(&U1024::from(2)).0;
    let mut r = a.mod_pow(&half_q_plus_1, p);

    loop {
        if t == U1024::ONE {
            return Some(r);
        }

        let mut i = 1;
        let mut temp = t.mod_mul(&t, p);
        while temp != U1024::ONE {
            temp = temp.mod_mul(&temp, p);
            i += 1;
        }

        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = b.mod_mul(&b, p);
        }

        m = i;
        c = b.mod_mul(&b, p);
        t = t.mod_mul(&c, p);
        r = r.mod_mul(&b, p);
    }
}

/// Compute gcd(a, b) using binary GCD.
fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        return a;
    } else {
        return gcd(b, a % b);
    }
}

fn bit_length(n: &U1024) -> usize {
    for i in (0..LIMBS).rev() {
        if n.0[i] != 0 {
            return i * 64 + (64 - n.0[i].leading_zeros()) as usize;
        }
    }

    0
}

/// Search for pairing-friendly curve parameters.
///
/// Inputs:
///   k        — embedding degree
///   d        — CM discriminant |D|  (we use -D internally)
///   t_start  — start of T search range
///   t_max    — end of T search range
///   lambda_r — target bit-length for r
///   lambda_p — target bit-length for p
fn cocks_pinch(
    k: u64,
    d: &U1024,
    t_start: u64,
    t_max: u64,
    lambda_r: usize,
    lambda_p: usize,
) -> Option<CurveParams> {
    for t_val in t_start..=t_max {
        let t_u1024 = U1024::from(t_val);

        // 6a. r = Φ_k(T)
        let r = cyclotomic_eval(k, &t_u1024);

        // 6b. Skip if r is not prime
        if !is_prime(&r) {
            continue;
        }

        // 6c. Skip if bit_length(r) != lambda_r
        if bit_length(&r) != lambda_r {
            continue;
        }

        // 6d. Compute sqrt(-D) mod r, skip if not a QR
        let sqrt_neg_d = match sqrt_mod(d, &r) {
            Some(v) => v,
            None => continue,
        };

        // 6e. Try each i with gcd(i, k) == 1
        for i in 1..k {
            if gcd(i, k) != 1 {
                continue;
            }

            // t_0 = T^i + 1 (mod r)
            let t0 = compute_t0(&t_u1024, i, &r);

            // y_0 = (t_0 - 2) / sqrt(-D) (mod r)
            let y0 = compute_y0(&t0, &sqrt_neg_d, &r);

            // Try random h_t, h_y to inflate p to lambda_p bits
            if let Some(params) = lift_to_prime(k, d, &r, &t0, &y0, lambda_p) {
                return Some(params);
            }
        }
    }

    None
}

/// Compute t_0 = T^i + 1 (mod r)
fn compute_t0(t_val: &U1024, i: u64, r: &U1024) -> U1024 {
    t_val.mod_pow(&U1024::from(i), r).mod_add(&U1024::ONE, r)
}

/// Compute y_0 = (t_0 - 2) * (sqrt(-D))^{-1} (mod r)
fn compute_y0(t0: &U1024, sqrt_neg_d: &U1024, r: &U1024) -> U1024 {
    let t_minus_2 = t0.borrowing_sub(&U1024::from(2)).0;
    t_minus_2.mod_mul(&sqrt_neg_d.mod_inverse(&r), &r)
}

/// Lift (t_0, y_0) with random multipliers h_t, h_y to find
/// t = t_0 + h_t * r,  y = y_0 + h_y * r
/// such that p = (t^2 + D*y^2)/4 is prime and has lambda_p bits.
fn lift_to_prime(
    k: u64,
    d: &U1024,
    r: &U1024,
    t0: &U1024,
    y0: &U1024,
    lambda_p: usize,
) -> Option<CurveParams> {
    let target_half = lambda_p / 2 + 1;

    let two_target = U1024::ONE.shl(target_half);
    let h_base = two_target.div_rem(r).0;

    for _ in 0..500 {
        let h_t = h_base
            .carrying_add(&U1024::from(rand::rng().random_range(1..500)))
            .0;
        let h_y = h_base
            .carrying_add(&U1024::from(rand::rng().random_range(1..500)))
            .0;

        let t = t0.carrying_add(&h_t.widening_mul(r).0).0;
        let y = y0.carrying_add(&h_y.widening_mul(r).0).0;

        let t_sq = t.widening_mul(&t).0;
        let y_sq = y.widening_mul(&y).0;
        let dy_sq = d.widening_mul(&y_sq).0;
        let numerator = t_sq.carrying_add(&dy_sq).0;

        if numerator.0[0] & 3 != 0 {
            continue;
        }

        let p = numerator.shr(2);
        if bit_length(&p) != lambda_p {
            continue;
        }

        if !is_prime(&p) {
            continue;
        }

        return Some(CurveParams {
            p,
            r: *r,
            t,
            y,
            k,
            d: *d,
        });
    }

    None
}

fn main() {
    println!("=== Cocks-Pinch Curve Parameter Generator ===\n");

    let k = 18;
    let d = U1024::from(3);
    let t_start = 4_000_000_000_000u64;
    let t_max = 5_000_000_000_000u64;
    let lambda_r = 256; // target r bit-length
    let lambda_p = 1024; // target p bit-length

    println!("Searching with k={k}, D={d}, T in [{t_start}, {t_max}]");
    println!("Target: r ~ {lambda_r} bits, p ~ {lambda_p} bits\n");

    match cocks_pinch(k, &d, t_start, t_max, lambda_r, lambda_p) {
        Some(params) => {
            println!("Found pairing-friendly curve!");
            println!("  p = {}", params.p); println!("  r = {}", params.r);
            println!("  t = {}", params.t);
            println!("  y = {}", params.y);
            println!("  k = {}", params.k);
            println!("  D = {}", params.d);
            println!("  p bits = {}", bit_length(&params.p));
            println!("  r bits = {}", bit_length(&params.r));

            // Write to curve.toml
            let config = CurveConfig {
                modulus: format!("{}", params.p),
                order: format!("{}", params.r),
                trace: format!("{}", params.t),
                cm_y: format!("{}", params.y),
                embedding_degree: params.k,
                discriminant: format!("{}", params.d),
            };

            let toml_str = toml::to_string_pretty(&config).expect("Failed to serialize config");
            std::fs::write("curve.toml", &toml_str).expect("Failed to write curve.toml");
            println!("\nConfig written to curve.toml");
        }
        None => {
            println!("No curve found in the given search range.");
        }
    }
}

#[derive(serde::Serialize, serde::Deserialize)]
struct CurveConfig {
    modulus: String,       // p as hex
    order: String,         // r as hex (group order)
    trace: String,         // t as hex (trace of Frobenius)
    cm_y: String,          // y as hex
    embedding_degree: u64, // k
    discriminant: String,  // D as hex
}
