use curve1024::runtime_field::RuntimeField;
use curve1024::{U1024, u1024::LIMBS};

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

const SMALL_PRIMES: [u64; 100] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
];

/// Miller-Rabin primality test with small-prime sieve.
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

    for &p in &SMALL_PRIMES {
        let sp = U1024::from(p);
        if *n == sp {
            return true;
        }
        if n.div_rem(&sp).1.is_zero() {
            return false;
        }
    }

    let field = RuntimeField::new(*n);
    let mut d = n.borrowing_sub(&U1024::ONE).0;
    let mut s = 0u64;
    while d.0[0] & 1 == 0 {
        d = d.shr(1);
        s += 1;
    }

    let n_minus_1 = n.borrowing_sub(&U1024::ONE).0;

    for i in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let a = U1024::from(i);
        if a >= *n {
            continue;
        }

        let mut x = field.pow(&a, &d);
        if x == U1024::ONE || x == n_minus_1 {
            continue;
        }

        let mut composite = true;
        for _ in 0..s.saturating_sub(1) {
            x = field.mul(&x, &x);
            if x == n_minus_1 {
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
    let field = RuntimeField::new(*p);
    let half_p_minus_1 = p.borrowing_sub(&U1024::ONE).0.div_rem(&U1024::from(2)).0;
    if field.pow(a, &half_p_minus_1) != U1024::ONE {
        return None;
    }

    let mut q = p.borrowing_sub(&U1024::ONE).0;
    let mut s = 0;
    while q.0[0] & 1 == 0 {
        q = q.shr(1);
        s += 1;
    }

    let mut z = U1024::from(2);
    while field.pow(&z, &half_p_minus_1) != p.borrowing_sub(&U1024::ONE).0 {
        z = z.carrying_add(&U1024::ONE).0;
    }

    let mut m = s;
    let mut c = field.pow(&z, &q);
    let mut t = field.pow(a, &q);
    let half_q_plus_1 = q.carrying_add(&U1024::ONE).0.div_rem(&U1024::from(2)).0;
    let mut r = field.pow(a, &half_q_plus_1);

    loop {
        if t == U1024::ONE {
            return Some(r);
        }

        let mut i = 1;
        let mut temp = field.mul(&t, &t);
        while temp != U1024::ONE && i < m {
            temp = field.mul(&temp, &temp);
            i += 1;
        }

        if i == m {
            return None; // Not a quadratic residue
        }

        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = field.mul(&b, &b);
        }

        m = i;
        c = field.mul(&b, &b);
        t = field.mul(&t, &c);
        r = field.mul(&r, &b);
    }
}

/// Compute gcd(a, b) using binary GCD.
fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
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
///   k           — embedding degree
///   d           — CM discriminant |D|
///   lambda_r    — target bit-length for r
///   lambda_p    — target bit-length for p
///   max_attempts — max random T values to try
fn cocks_pinch(
    k: u64,
    d: &U1024,
    lambda_r: usize,
    lambda_p: usize,
    max_attempts: u64,
) -> Option<CurveParams> {
    eprintln!("Computing T range via binary search...");
    let setup_time = std::time::Instant::now();

    // Find exact T range where Φ_k(T) has exactly lambda_r bits via binary search
    let t_min = {
        let mut lo = U1024::ONE;
        let mut hi = U1024::ONE.shl(lambda_r / 6 + 2);
        while lo < hi {
            let mid = lo.carrying_add(&hi).0.shr(1);
            if bit_length(&cyclotomic_eval(k, &mid)) < lambda_r {
                lo = mid.carrying_add(&U1024::ONE).0;
            } else {
                hi = mid;
            }
        }
        lo
    };
    eprintln!("  T_min found in {:.2?}", setup_time.elapsed());

    let t_max = {
        let mut lo = t_min;
        let mut hi = U1024::ONE.shl(lambda_r / 6 + 2);
        while lo < hi {
            let mid = lo.carrying_add(&hi).0.shr(1);
            if bit_length(&cyclotomic_eval(k, &mid)) <= lambda_r {
                lo = mid.carrying_add(&U1024::ONE).0;
            } else {
                hi = mid;
            }
        }
        lo
    };
    eprintln!("  T_max found in {:.2?}", setup_time.elapsed());

    println!("  T_min = {}", t_min);
    println!("  T_max = {}", t_max);
    println!(
        "  T range width: ~{} bits",
        bit_length(&t_max.borrowing_sub(&t_min).0)
    );

    for attempt in 0..max_attempts {
        let t0_time = std::time::Instant::now();

        let t_range = t_max.borrowing_sub(&t_min).0;
        let t_u1024 = t_min.carrying_add(&U1024::random_below(&t_range)).0;

        let r = cyclotomic_eval(k, &t_u1024);

        let bl = bit_length(&r);
        if bl != lambda_r {
            continue;
        }

        // Quick small-prime sieve before full is_prime
        let mut sieve_fail = false;
        for &p in &SMALL_PRIMES {
            let sp = U1024::from(p);
            if r != sp && r.div_rem(&sp).1.is_zero() {
                sieve_fail = true;
                break;
            }
        }
        if sieve_fail {
            continue;
        }

        if !is_prime(&r) {
            continue;
        }
        eprintln!("PRIME! {:.2?}", t0_time.elapsed());

        println!(
            "[attempt {attempt}] Found prime r, r bits={}, took {:.2?}",
            bit_length(&r),
            t0_time.elapsed()
        );

        // Compute sqrt(-D) mod r, skip if not a QR
        let sqrt_neg_d = match sqrt_mod(d, &r) {
            Some(v) => v,
            None => {
                println!(
                    "[attempt {attempt}] -D is not a QR mod r, took {:.2?}",
                    t0_time.elapsed()
                );
                continue;
            }
        };

        let r_field = RuntimeField::new(r);

        // Try each i with gcd(i, k) == 1
        for i in 1..k {
            if gcd(i, k) != 1 {
                continue;
            }

            // t_0 = T^i + 1 (mod r)
            let t0 = compute_t0(&t_u1024, i, &r_field);

            // y_0 = (t_0 - 2) / sqrt(-D) (mod r)
            let y0 = compute_y0(&t0, &sqrt_neg_d, &r_field);

            // Try random h_t, h_y to inflate p to lambda_p bits
            if let Some(params) = lift_to_prime(k, d, &r, &t0, &y0, lambda_p) {
                println!(
                    "[attempt {attempt}] SUCCESS, total attempt took {:.2?}",
                    t0_time.elapsed()
                );
                return Some(params);
            }
        }

        println!(
            "[attempt {attempt}] lift_to_prime failed, took {:.2?}",
            t0_time.elapsed()
        );
    }

    None
}

/// Compute t_0 = T^i + 1 (mod r)
fn compute_t0(t_val: &U1024, i: u64, field: &RuntimeField) -> U1024 {
    field.add(&field.pow(t_val, &U1024::from(i)), &U1024::ONE)
}

/// Compute y_0 = (t_0 - 2) * (sqrt(-D))^{-1} (mod r)
fn compute_y0(t0: &U1024, sqrt_neg_d: &U1024, field: &RuntimeField) -> U1024 {
    let t_minus_2 = t0.borrowing_sub(&U1024::from(2)).0;
    field.mul(&t_minus_2, &field.inv(sqrt_neg_d))
}

fn add_2048(a_lo: &U1024, a_hi: &U1024, b_lo: &U1024, b_hi: &U1024) -> (U1024, U1024) {
    let (lo, carry) = a_lo.carrying_add(b_lo);
    let mut hi = a_hi.carrying_add(b_hi).0;
    if carry {
        hi = hi.carrying_add(&U1024::ONE).0;
    }
    (lo, hi)
}

fn shr_2048(lo: &U1024, hi: &U1024, shift: usize) -> (U1024, U1024) {
    let new_hi = hi.shr(shift);
    let mut new_lo = lo.shr(shift);
    let bottom_bits_of_hi = hi.0[0] & ((1 << shift) - 1);
    new_lo.0[15] |= bottom_bits_of_hi << (64 - shift);
    (new_lo, new_hi)
}

fn get_lift(val0: &U1024, r: &U1024, h: i32) -> U1024 {
    if h >= 0 {
        let hr = U1024::from(h as u64).widening_mul(r).0;
        hr.carrying_add(val0).0
    } else {
        let hr = U1024::from((-h) as u64).widening_mul(r).0;
        let (diff, borrow) = val0.borrowing_sub(&hr);
        if borrow {
            hr.borrowing_sub(val0).0
        } else {
            diff
        }
    }
}

/// Lift (t_0, y_0) with small signed multipliers h_t, h_y to find
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
    for ht in -20..=20 {
        for hy in -20..=20 {
            let t = get_lift(t0, r, ht);
            let y = get_lift(y0, r, hy);

            let (tsq_lo, tsq_hi) = t.widening_mul(&t);
            let (ysq_lo, ysq_hi) = y.widening_mul(&y);

            let mut dy_sq_lo = U1024::ZERO;
            let mut dy_sq_hi = U1024::ZERO;
            for _ in 0..d.0[0] {
                let (lo, hi) = add_2048(&dy_sq_lo, &dy_sq_hi, &ysq_lo, &ysq_hi);
                dy_sq_lo = lo;
                dy_sq_hi = hi;
            }

            let (num_lo, num_hi) = add_2048(&tsq_lo, &tsq_hi, &dy_sq_lo, &dy_sq_hi);

            if num_lo.0[0] & 3 != 0 {
                continue;
            }

            let (p_lo, p_hi) = shr_2048(&num_lo, &num_hi, 2);

            if p_hi != U1024::ZERO {
                continue; // Overflowed U1024 entirely
            }

            let p = p_lo;
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
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqrt_mod() {
        let p = U1024::from(65537u64); // prime, p ≡ 1 mod 4
        let a = U1024::from(3u64); // not QR (3 is a non-residue mod 65537)
        let b = U1024::from(4u64); // QR, sqrt is 2
        let c = U1024::from(9u64); // QR, sqrt is 3

        // 3^32768 mod 65537 == 65536 (-1), so 3 is not a QR.
        assert_eq!(sqrt_mod(&a, &p), None);

        let sqrt_b = sqrt_mod(&b, &p).unwrap();
        assert!(sqrt_b.0[0] == 2 || p.borrowing_sub(&sqrt_b).0.0[0] == 2);

        let sqrt_c = sqrt_mod(&c, &p).unwrap();
        assert!(sqrt_c.0[0] == 3 || p.borrowing_sub(&sqrt_c).0.0[0] == 3);

        let p7 = U1024::from(7u64);
        let d = U1024::from(2u64); // QR mod 7, sqrt is 3 or 4
        let sqrt_d = sqrt_mod(&d, &p7).unwrap();
        assert!(sqrt_d.0[0] == 3 || sqrt_d.0[0] == 4);
    }
}

fn main() {
    println!("=== Cocks-Pinch Curve Parameter Generator ===\n");

    let k = 18;
    let d = U1024::from(3);
    let lambda_r = 512; // target r bit-length
    let lambda_p = 1024; // target p bit-length
    let max_attempts = 100_000; // max random T values to try

    println!("k={k}, D={d}, target: r~{lambda_r} bits, p~{lambda_p} bits");
    println!("Max attempts: {max_attempts}\n");

    let start = std::time::Instant::now();

    match cocks_pinch(k, &d, lambda_r, lambda_p, max_attempts) {
        Some(params) => {
            println!("Found pairing-friendly curve!");
            println!("  p = {}", params.p);
            println!("  r = {}", params.r);
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
            println!("Total time: {:.2?}", start.elapsed());
        }
        None => {
            println!("No curve found after {max_attempts} attempts.");
            println!("Total time: {:.2?}", start.elapsed());
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
