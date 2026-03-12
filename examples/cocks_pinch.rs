use curve1024::{U1024, u1024::LIMBS};

const CONFIG_PATH: &str = "config/curve1024.toml";

const SMALL_PRIMES: [u64; 100] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
];

const MILLER_RABIN_WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

struct CurveParams {
    p: U1024,
    r: U1024,
    t: U1024,
    y: U1024,
    k: u64,
    d: U1024,
}

impl CurveParams {
    fn to_toml(&self) -> String {
        format!(
            "[field]\nmodulus = \"{p}\"\norder = \"{r}\"\n\n[cm]\ntrace = \"{t}\"\ny = \"{y}\"\nembedding_degree = {k}\ndiscriminant = \"{d}\"\n",
            p = self.p,
            r = self.r,
            t = self.t,
            y = self.y,
            k = self.k,
            d = self.d
        )
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

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn passes_small_prime_sieve(n: &U1024) -> bool {
    for &p in &SMALL_PRIMES {
        let sp = U1024::from(p);
        if *n == sp {
            return true;
        }
        if n.div_rem(&sp).1.is_zero() {
            return false;
        }
    }
    true
}

fn is_prime(n: &U1024) -> bool {
    if *n == U1024::ONE || n.0[0] & 1 == 0 {
        return *n == U1024::from(2);
    }
    if !passes_small_prime_sieve(n) {
        return false;
    }

    let field = RuntimeField::new(*n);
    let n_minus_1 = n.borrowing_sub(&U1024::ONE).0;

    // Factor n-1 = d * 2^s
    let mut d = n_minus_1;
    let mut s = 0u64;
    while d.0[0] & 1 == 0 {
        d = d.shr(1);
        s += 1;
    }

    for &witness in &MILLER_RABIN_WITNESSES {
        let a = U1024::from(witness);
        if a >= *n {
            continue;
        }

        let mut x = field.pow(&a, &d);
        if x == U1024::ONE || x == n_minus_1 {
            continue;
        }

        let mut found_minus_one = false;
        for _ in 0..s.saturating_sub(1) {
            x = field.mul(&x, &x);
            if x == n_minus_1 {
                found_minus_one = true;
                break;
            }
        }

        if !found_minus_one {
            return false;
        }
    }
    true
}

fn sqrt_mod(a: &U1024, p: &U1024) -> Option<U1024> {
    let field = RuntimeField::new(*p);
    let p_minus_1 = p.borrowing_sub(&U1024::ONE).0;
    let half_p_minus_1 = p_minus_1.shr(1);

    // Euler criterion
    if field.pow(a, &half_p_minus_1) != U1024::ONE {
        return None;
    }

    // Factor p-1 = q * 2^s
    let mut q = p_minus_1;
    let mut s = 0u32;
    while q.0[0] & 1 == 0 {
        q = q.shr(1);
        s += 1;
    }

    // Find a quadratic non-residue
    let mut z = U1024::from(2);
    while field.pow(&z, &half_p_minus_1) != p_minus_1 {
        z = z.carrying_add(&U1024::ONE).0;
    }

    let mut m = s;
    let mut c = field.pow(&z, &q);
    let mut t = field.pow(a, &q);
    let q_plus_1_half = q.carrying_add(&U1024::ONE).0.shr(1);
    let mut root = field.pow(a, &q_plus_1_half);

    loop {
        if t == U1024::ONE {
            return Some(root);
        }

        // Find least i such that t^(2^i) ≡ 1
        let mut i = 1u32;
        let mut temp = field.mul(&t, &t);
        while temp != U1024::ONE && i < m {
            temp = field.mul(&temp, &temp);
            i += 1;
        }

        if i == m {
            return None;
        }

        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = field.mul(&b, &b);
        }

        m = i;
        c = field.mul(&b, &b);
        t = field.mul(&t, &c);
        root = field.mul(&root, &b);
    }
}

fn cyclotomic_phi18(t: &U1024) -> U1024 {
    let t2 = t.widening_mul(t).0;
    let t3 = t2.widening_mul(t).0;
    let t6 = t3.widening_mul(&t3).0;
    t6.borrowing_sub(&t3).0.carrying_add(&U1024::ONE).0
}

fn add_2048(a: (&U1024, &U1024), b: (&U1024, &U1024)) -> (U1024, U1024) {
    let (lo, carry) = a.0.carrying_add(b.0);
    let mut hi = a.1.carrying_add(b.1).0;
    if carry {
        hi = hi.carrying_add(&U1024::ONE).0;
    }
    (lo, hi)
}

fn shr_2048(lo: &U1024, hi: &U1024, shift: usize) -> (U1024, U1024) {
    let new_hi = hi.shr(shift);
    let mut new_lo = lo.shr(shift);
    let spill = hi.0[0] & ((1 << shift) - 1);
    new_lo.0[LIMBS - 1] |= spill << (64 - shift);
    (new_lo, new_hi)
}

fn apply_lift(base: &U1024, r: &U1024, h: i32) -> U1024 {
    if h >= 0 {
        let offset = U1024::from(h as u64).widening_mul(r).0;
        offset.carrying_add(base).0
    } else {
        let offset = U1024::from((-h) as u64).widening_mul(r).0;
        let (diff, borrow) = base.borrowing_sub(&offset);
        if borrow {
            offset.borrowing_sub(base).0
        } else {
            diff
        }
    }
}

/// Returns the largest s such that 2^s divides n.
fn two_adicity(n: &U1024) -> u32 {
    for limb_idx in 0..LIMBS {
        if n.0[limb_idx] != 0 {
            return (limb_idx as u32) * 64 + n.0[limb_idx].trailing_zeros();
        }
    }
    1024
}

struct LiftParams<'a> {
    t0: &'a U1024,
    y0: &'a U1024,
    target_p_bits: usize,
    /// Minimum two-adicity of p - 1  (p = d·2^x + 1 form)
    min_base_two_adicity: u32,
    /// Actual two-adicity of r - 1; controls how far ht·r shifts bits of p.
    r_two_adicity: u32,
}

fn try_lift_to_prime(k: u64, d: &U1024, r: &U1024, lp: &LiftParams<'_>) -> Option<CurveParams> {
    let d_small = d.0[0];
    let extra = lp
        .min_base_two_adicity
        .saturating_sub(lp.r_two_adicity)
        .min(10);
    let half_range = 20i64 + (1i64 << extra);

    for ht in -half_range..=half_range {
        for hy in -half_range..=half_range {
            let t = apply_lift(lp.t0, r, ht as i32);
            let y = apply_lift(lp.y0, r, hy as i32);

            let t_sq = t.widening_mul(&t);
            let y_sq = y.widening_mul(&y);

            let mut d_y_sq = (U1024::ZERO, U1024::ZERO);
            for _ in 0..d_small {
                d_y_sq = add_2048((&d_y_sq.0, &d_y_sq.1), (&y_sq.0, &y_sq.1));
            }

            let numerator = add_2048((&t_sq.0, &t_sq.1), (&d_y_sq.0, &d_y_sq.1));

            if numerator.0.0[0] & 3 != 0 {
                continue;
            }

            let (p, p_hi) = shr_2048(&numerator.0, &numerator.1, 2);
            if p_hi != U1024::ZERO || bit_length(&p) != lp.target_p_bits {
                continue;
            }

            if two_adicity(&p.borrowing_sub(&U1024::ONE).0) < lp.min_base_two_adicity {
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

fn find_t_range(target_bits: usize) -> (U1024, U1024) {
    let upper_bound = U1024::ONE.shl(target_bits / 6 + 2);

    let t_min = {
        let (mut lo, mut hi) = (U1024::ONE, upper_bound);
        while lo < hi {
            let mid = lo.carrying_add(&hi).0.shr(1);
            if bit_length(&cyclotomic_phi18(&mid)) < target_bits {
                lo = mid.carrying_add(&U1024::ONE).0;
            } else {
                hi = mid;
            }
        }
        lo
    };

    let t_max = {
        let (mut lo, mut hi) = (t_min, upper_bound);
        while lo < hi {
            let mid = lo.carrying_add(&hi).0.shr(1);
            if bit_length(&cyclotomic_phi18(&mid)) <= target_bits {
                lo = mid.carrying_add(&U1024::ONE).0;
            } else {
                hi = mid;
            }
        }
        lo
    };

    (t_min, t_max)
}

fn cocks_pinch(
    k: u64,
    d: &U1024,
    target_r_bits: usize,
    target_p_bits: usize,
    min_scalar_two_adicity: u32,
    min_base_two_adicity: u32,
    max_attempts: u64,
) -> Option<CurveParams> {
    let (t_min, t_max) = find_t_range(target_r_bits);

    // r = T^6 - T^3 + 1, so r - 1 = T^3(T^3 - 1).
    // If T ≡ 0 (mod 2^k), then two_adicity(r-1) = 3k.
    // => align T to multiples of 2^ceil(s/3) to guarantee two_adicity(r-1) >= s.
    let t_align = min_scalar_two_adicity.div_ceil(3);
    let step = U1024::ONE.shl(t_align as usize);
    // Snap t_min up to the next multiple of step
    let t_base = {
        let rem = t_min.div_rem(&step).1;
        if rem.is_zero() {
            t_min
        } else {
            t_min.carrying_add(&step).0.borrowing_sub(&rem).0
        }
    };
    let t_steps = t_max.borrowing_sub(&t_base).0.div_rem(&step).0;

    for attempt in 0..max_attempts {
        let timer = std::time::Instant::now();

        // Pick a random multiple of `step` in [t_base, t_max]
        let t_val = t_base
            .carrying_add(&U1024::rand(&t_steps).widening_mul(&step).0)
            .0;
        let r = cyclotomic_phi18(&t_val);

        if bit_length(&r) != target_r_bits {
            continue;
        }
        if !passes_small_prime_sieve(&r) {
            continue;
        }
        if !is_prime(&r) {
            continue;
        }

        let r_two_adicity = two_adicity(&r.borrowing_sub(&U1024::ONE).0);
        println!(
            "[attempt {attempt}] Found prime r ({} bits, two-adicity={}), {:.2?}",
            bit_length(&r),
            r_two_adicity,
            timer.elapsed()
        );

        let neg_d = r.borrowing_sub(d).0;
        let sqrt_neg_d = match sqrt_mod(&neg_d, &r) {
            Some(v) => v,
            None => {
                println!("[attempt {attempt}] -D is not a QR mod r, skipping");
                continue;
            }
        };

        let r_field = RuntimeField::new(r);

        for i in 1..k {
            if gcd(i, k) != 1 {
                continue;
            }

            let t0 = r_field.add(&r_field.pow(&t_val, &U1024::from(i)), &U1024::ONE);
            let t0_minus_2 = t0.borrowing_sub(&U1024::from(2)).0;
            let y0 = r_field.mul(&t0_minus_2, &r_field.inv(&sqrt_neg_d));

            if let Some(params) = try_lift_to_prime(
                k,
                d,
                &r,
                &LiftParams {
                    t0: &t0,
                    y0: &y0,
                    target_p_bits,
                    min_base_two_adicity,
                    r_two_adicity,
                },
            ) {
                println!("[attempt {attempt}] SUCCESS, {:.2?}", timer.elapsed());
                return Some(params);
            }
        }

        println!(
            "[attempt {attempt}] lift_to_prime failed, {:.2?}",
            timer.elapsed()
        );
    }

    None
}

fn main() {
    println!("=== Cocks-Pinch Curve Parameter Generator ===\n");

    let k = 18u64;
    let d = U1024::from(3);
    let target_r_bits = 512;
    let target_p_bits = 1024;
    let max_attempts = 100_000u64;
    let min_scalar_two_adicity = 32u32;
    let min_base_two_adicity = 32u32;

    println!("k={k}, D={d}, target: r~{target_r_bits} bits, p~{target_p_bits} bits");
    println!("Scalar field NTT two-adicity >= {min_scalar_two_adicity}");
    println!("Base field NTT two-adicity   >= {min_base_two_adicity}");
    println!("Max attempts: {max_attempts}\n");

    let start = std::time::Instant::now();

    match cocks_pinch(
        k,
        &d,
        target_r_bits,
        target_p_bits,
        min_scalar_two_adicity,
        min_base_two_adicity,
        max_attempts,
    ) {
        Some(params) => {
            let r_adicity = two_adicity(&params.r.borrowing_sub(&U1024::ONE).0);
            let p_adicity = two_adicity(&params.p.borrowing_sub(&U1024::ONE).0);
            println!("Found pairing-friendly curve!");
            println!("  p = {} ({} bits)", params.p, bit_length(&params.p));
            println!("  p two-adicity: 2^{p_adicity} | (p-1)");
            println!("  r = {} ({} bits)", params.r, bit_length(&params.r));
            println!("  r two-adicity: 2^{r_adicity} | (r-1)");
            println!("  t = {}", params.t);
            println!("  y = {}", params.y);
            println!("  k = {}", params.k);
            println!("  D = {}", params.d);

            std::fs::write(CONFIG_PATH, params.to_toml()).expect("Failed to write config");
            println!("\nConfig written to {CONFIG_PATH}");
        }
        None => {
            println!("No curve found after {max_attempts} attempts.");
        }
    }

    println!("Total time: {:.2?}", start.elapsed());
}

struct RuntimeField {
    modulus: U1024,
    n_prime: U1024,
    r2: U1024,
}

impl RuntimeField {
    fn new(modulus: U1024) -> Self {
        let n_prime = Self::compute_n_prime(&modulus);
        let r2 = Self::compute_r2(&modulus);
        Self {
            modulus,
            n_prime,
            r2,
        }
    }

    fn reduce(&self, a: &U1024) -> U1024 {
        a.div_rem(&self.modulus).1
    }

    fn add(&self, a: &U1024, b: &U1024) -> U1024 {
        let (sum, carry) = a.carrying_add(b);
        let (sub_res, borrow) = sum.borrowing_sub(&self.modulus);
        if carry || !borrow { sub_res } else { sum }
    }

    fn mul(&self, a: &U1024, b: &U1024) -> U1024 {
        let a_mont = self.encode_mont(a);
        let b_mont = self.encode_mont(b);
        self.decode_mont(&self.mont_mul(&a_mont, &b_mont))
    }

    fn pow(&self, base: &U1024, exp: &U1024) -> U1024 {
        if self.modulus == U1024::ONE {
            return U1024::ZERO;
        }

        let mut result = self.encode_mont(&U1024::ONE);
        let mut b = self.encode_mont(&self.reduce(base));

        let top = (0..LIMBS)
            .rev()
            .find(|&i| exp.0[i] != 0)
            .map_or(0, |i| i + 1);

        for i in 0..top {
            let mut limb = exp.0[i];
            for _ in 0..64 {
                if limb & 1 == 1 {
                    result = self.mont_mul(&result, &b);
                }
                b = self.mont_mul(&b, &b);
                limb >>= 1;
            }
        }

        self.decode_mont(&result)
    }

    fn inv(&self, a: &U1024) -> U1024 {
        assert!(!a.is_zero());
        let exp = self.modulus.borrowing_sub(&U1024::from(2)).0;
        self.pow(a, &exp)
    }

    fn encode_mont(&self, a: &U1024) -> U1024 {
        self.mont_mul(a, &self.r2)
    }

    fn decode_mont(&self, a: &U1024) -> U1024 {
        Self::mont_reduce(a, &U1024::ZERO, &self.modulus, &self.n_prime)
    }

    fn mont_mul(&self, a: &U1024, b: &U1024) -> U1024 {
        let (lo, hi) = a.widening_mul(b);
        Self::mont_reduce(&lo, &hi, &self.modulus, &self.n_prime)
    }

    fn mont_reduce(lo: &U1024, hi: &U1024, n: &U1024, n_prime: &U1024) -> U1024 {
        let (m, _) = lo.widening_mul(n_prime);
        let (mn_lo, mn_hi) = m.widening_mul(n);
        let (_, carry_lo) = lo.carrying_add(&mn_lo);
        let (res_hi, carry_hi) = hi.carrying_add(&mn_hi);
        let (t, carry_final) = if carry_lo {
            res_hi.carrying_add(&U1024::ONE)
        } else {
            (res_hi, false)
        };
        if carry_hi || carry_final {
            return t.borrowing_sub(n).0;
        }
        let (sub, borrow) = t.borrowing_sub(n);
        if !borrow { sub } else { t }
    }

    fn compute_n_prime(n: &U1024) -> U1024 {
        let mut inv = U1024::ONE;
        let two = U1024::from(2);
        for _ in 0..10 {
            let ni = n.widening_mul(&inv).0;
            let diff = two.borrowing_sub(&ni).0;
            inv = inv.widening_mul(&diff).0;
        }
        U1024::ZERO.borrowing_sub(&inv).0
    }

    fn compute_r2(n: &U1024) -> U1024 {
        let neg_n = U1024::ZERO.borrowing_sub(n).0;
        let r_mod_n = neg_n.div_rem(n).1;
        Self::mod_mul_raw(&r_mod_n, &r_mod_n, n)
    }

    // Schoolbook modular multiply for bootstrapping (before Montgomery is available)
    fn mod_mul_raw(a: &U1024, b: &U1024, n: &U1024) -> U1024 {
        let (lo, hi) = a.widening_mul(b);
        if hi == U1024::ZERO {
            return lo.div_rem(n).1;
        }
        let r = U1024::ZERO.borrowing_sub(n).0.div_rem(n).1; // 2^1024 mod n
        let hi_reduced = Self::mod_mul_raw(&r, &hi, n);
        let lo_reduced = lo.div_rem(n).1;
        let (sum, carry) = lo_reduced.carrying_add(&hi_reduced);
        let (sub_res, borrow) = sum.borrowing_sub(n);
        if carry || !borrow { sub_res } else { sum }
    }
}
