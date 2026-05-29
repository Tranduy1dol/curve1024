use curve1024::{U1024, u1024::LIMBS};
use serde::Deserialize;

const CONFIG_PATH: &str = "config/curve1024.toml";

#[derive(Deserialize)]
struct ConfigFile {
    field: FieldConfig,
    cm: CmConfig,
}

#[derive(Deserialize)]
struct FieldConfig {
    modulus: String,
    order: String,
}

#[derive(Deserialize)]
struct CmConfig {
    trace: String,
    y: String,
    discriminant: String,
}

struct CurveConfig {
    p: U1024,
    r: U1024,
    t: U1024,
    y: U1024,
    d: U1024,
}

fn read_config() -> CurveConfig {
    let content = std::fs::read_to_string(CONFIG_PATH)
        .expect("config/curve1024.toml not found");
    let file: ConfigFile = toml::from_str(&content).expect("Failed to parse TOML");
    CurveConfig {
        p: U1024::from_hex(&file.field.modulus),
        r: U1024::from_hex(&file.field.order),
        t: U1024::from_hex(&file.cm.trace),
        y: U1024::from_hex(&file.cm.y),
        d: U1024::from_hex(&file.cm.discriminant),
    }
}

fn append_to_config(a: &U1024, b: &U1024, g_x: &U1024, g_y: &U1024) {
    let content = std::fs::read_to_string(CONFIG_PATH).expect("Failed to read config");
    let mut new_content = String::new();
    for line in content.lines() {
        if line == "[curve]" {
            break;
        }
        new_content.push_str(line);
        new_content.push('\n');
    }
    new_content.push_str(&format!(
        "\n[curve]\na = \"{a}\"\nb = \"{b}\"\n\n[generator]\nx = \"{g_x}\"\ny = \"{g_y}\"\n"
    ));
    std::fs::write(CONFIG_PATH, new_content).expect("Failed to write to config");
}

#[derive(Clone, PartialEq)]
struct Point {
    x: U1024,
    y: U1024,
    is_infinity: bool,
}

impl Point {
    fn infinity() -> Self {
        Self {
            x: U1024::ZERO,
            y: U1024::ZERO,
            is_infinity: true,
        }
    }

    fn new(x: U1024, y: U1024) -> Self {
        Self {
            x,
            y,
            is_infinity: false,
        }
    }
}

fn point_add(p: &Point, q: &Point, a: &U1024, f: &RuntimeField) -> Point {
    if p.is_infinity {
        return q.clone();
    }
    if q.is_infinity {
        return p.clone();
    }
    if p.x == q.x {
        return if p.y == q.y && p.y != U1024::ZERO {
            point_double(p, a, f)
        } else {
            Point::infinity()
        };
    }

    let lambda = f.mul(&f.sub(&q.y, &p.y), &f.inv(&f.sub(&q.x, &p.x)));
    let x3 = f.sub(&f.sub(&f.mul(&lambda, &lambda), &p.x), &q.x);
    let y3 = f.sub(&f.mul(&lambda, &f.sub(&p.x, &x3)), &p.y);
    Point::new(x3, y3)
}

fn point_double(p: &Point, a: &U1024, f: &RuntimeField) -> Point {
    if p.is_infinity || p.y == U1024::ZERO {
        return Point::infinity();
    }

    // λ = (3x² + a) / (2y)
    let num = f.add(&f.mul(&U1024::from(3), &f.mul(&p.x, &p.x)), a);
    let lambda = f.mul(&num, &f.inv(&f.add(&p.y, &p.y)));
    let x3 = f.sub(&f.mul(&lambda, &lambda), &f.add(&p.x, &p.x));
    let y3 = f.sub(&f.mul(&lambda, &f.sub(&p.x, &x3)), &p.y);
    Point::new(x3, y3)
}

fn scalar_mul(p: &Point, scalar: &U1024, a: &U1024, f: &RuntimeField) -> Point {
    let mut result = Point::infinity();
    let mut base = p.clone();
    for i in 0..1024 {
        if (scalar.0[i / 64] >> (i % 64)) & 1 == 1 {
            result = point_add(&result, &base, a, f);
        }
        base = point_double(&base, a, f);
    }
    result
}

// Tonelli-Shanks square root mod p. Returns None if a is not a QR.
fn sqrt_mod(a: &U1024, f: &RuntimeField) -> Option<U1024> {
    let p_minus_1 = f.modulus.borrowing_sub(&U1024::ONE).0;
    let half = p_minus_1.shr(1);

    if f.pow(a, &half) != U1024::ONE {
        return None;
    }

    let mut q = p_minus_1;
    let mut s = 0u32;
    while q.0[0] & 1 == 0 {
        q = q.shr(1);
        s += 1;
    }

    let mut z = U1024::from(2);
    while f.pow(&z, &half) != p_minus_1 {
        z = z.carrying_add(&U1024::ONE).0;
    }

    let mut m = s;
    let mut c = f.pow(&z, &q);
    let mut t = f.pow(a, &q);
    let mut root = f.pow(a, &q.carrying_add(&U1024::ONE).0.shr(1));

    loop {
        if t == U1024::ONE {
            return Some(root);
        }

        let mut i = 1u32;
        let mut temp = f.mul(&t, &t);
        while temp != U1024::ONE && i < m {
            temp = f.mul(&temp, &temp);
            i += 1;
        }
        if i == m {
            return None;
        }

        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = f.mul(&b, &b);
        }
        m = i;
        c = f.mul(&b, &b);
        t = f.mul(&t, &c);
        root = f.mul(&root, &b);
    }
}

// Returns j-invariant from the Hilbert class polynomial H_D(x) mod p.
fn j_invariant(d: &U1024, f: &RuntimeField) -> U1024 {
    if *d == U1024::from(3) {
        U1024::ZERO // H_{-3}(x) = x  →  j = 0
    } else if *d == U1024::from(4) {
        f.reduce(&U1024::from(1728)) // H_{-4}(x) = x - 1728  →  j = 1728
    } else {
        panic!("Hilbert class polynomial not implemented for |D|={d}");
    }
}

// Returns initial (a, b) from j-invariant using the standard CM formulas.
fn curve_coeffs_from_j(j: &U1024, f: &RuntimeField) -> (U1024, U1024) {
    if *j == U1024::ZERO {
        (U1024::ZERO, U1024::ONE)
    } else if *j == U1024::from(1728) {
        (U1024::ONE, U1024::ZERO)
    } else {
        let c = f.mul(j, &f.inv(&f.sub(&U1024::from(1728), j)));
        (f.mul(&U1024::from(3), &c), f.mul(&U1024::from(2), &c))
    }
}

fn find_point_on_curve(a: &U1024, b: &U1024, f: &RuntimeField) -> Point {
    let mut x = U1024::ONE;
    loop {
        let rhs = f.add(&f.add(&f.mul(&f.mul(&x, &x), &x), &f.mul(a, &x)), b);
        if let Some(y) = sqrt_mod(&rhs, f) {
            return Point::new(x, y);
        }
        x = x.carrying_add(&U1024::ONE).0;
    }
}

fn find_twist(a: &U1024, config: &CurveConfig, f: &RuntimeField) -> (U1024, U1024, U1024) {
    let p1 = config.p.carrying_add(&U1024::ONE).0;
    let three_y = config.y.widening_mul(&U1024::from(3)).0;
    let t3p = config.t.carrying_add(&three_y).0.shr(1);
    let t3m = if config.t > three_y {
        config.t.borrowing_sub(&three_y).0.shr(1)
    } else {
        three_y.borrowing_sub(&config.t).0.shr(1)
    };

    let candidates = [
        p1.borrowing_sub(&config.t).0,
        p1.carrying_add(&config.t).0,
        p1.borrowing_sub(&t3p).0,
        p1.carrying_add(&t3p).0,
        p1.borrowing_sub(&t3m).0,
        p1.carrying_add(&t3m).0,
    ];

    let order = candidates
        .iter()
        .find_map(|&ord| {
            let (cof, rem) = ord.div_rem(&config.r);
            (rem == U1024::ZERO && cof != U1024::ZERO).then_some(ord)
        })
        .expect("None of the 6 CM twist orders is divisible by r");

    for b_val in 1u64..=2000 {
        let b = f.reduce(&U1024::from(b_val));
        let pt = find_point_on_curve(a, &b, f);
        if scalar_mul(&pt, &order, a, f).is_infinity {
            println!("  Found b = {b_val}");
            let cofactor = order.div_rem(&config.r).0;
            return (*a, b, cofactor);
        }
    }

    panic!("Could not find a valid b for any twist");
}

fn find_generator(a: &U1024, b: &U1024, cofactor: &U1024, r: &U1024, f: &RuntimeField) -> Point {
    loop {
        let g = scalar_mul(&find_point_on_curve(a, b, f), cofactor, a, f);
        if !g.is_infinity && scalar_mul(&g, r, a, f).is_infinity {
            return g;
        }
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

fn main() {
    println!("=== Complex Multiplication: Curve & Generator ===\n");

    let config = read_config();
    println!("p:   {} bits", bit_length(&config.p));
    println!("r:   {} bits", bit_length(&config.r));
    println!("|D|: {}\n", config.d);

    let f = RuntimeField::new(config.p);

    let j = j_invariant(&config.d, &f);
    println!("j-invariant: {j}\n");

    let (a_init, _) = curve_coeffs_from_j(&j, &f);

    println!("Finding correct twist (b)...");
    let (a, b, cofactor) = find_twist(&a_init, &config, &f);
    println!("  y² = x³ + {a}·x + {b}");
    println!("  cofactor = {cofactor} ({} bits)\n", bit_length(&cofactor));

    let timer = std::time::Instant::now();
    println!("Finding generator...");
    let g = find_generator(&a, &b, &cofactor, &config.r, &f);
    println!("  Done in {:.2?}", timer.elapsed());
    println!("  G.x = {}", g.x);
    println!("  G.y = {}\n", g.y);

    let rg = scalar_mul(&g, &config.r, &a, &f);
    assert!(rg.is_infinity, "[r]G must be the point at infinity");
    println!("✓ [r]G = O");

    let lhs = f.mul(&g.y, &g.y);
    let rhs = f.add(
        &f.add(&f.mul(&f.mul(&g.x, &g.x), &g.x), &f.mul(&a, &g.x)),
        &b,
    );
    assert!(lhs == rhs, "G must lie on the curve");
    println!("✓ G is on y² = x³ + ax + b\n");

    append_to_config(&a, &b, &g.x, &g.y);
    println!("Written to {CONFIG_PATH}");
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

    fn sub(&self, a: &U1024, b: &U1024) -> U1024 {
        let (diff, borrow) = a.borrowing_sub(b);
        if borrow {
            diff.carrying_add(&self.modulus).0
        } else {
            diff
        }
    }

    fn mul(&self, a: &U1024, b: &U1024) -> U1024 {
        self.decode_mont(&self.mont_mul(&self.encode_mont(a), &self.encode_mont(b)))
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
        self.pow(a, &self.modulus.borrowing_sub(&U1024::from(2)).0)
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
            inv = inv.widening_mul(&two.borrowing_sub(&ni).0).0;
        }
        U1024::ZERO.borrowing_sub(&inv).0
    }

    fn compute_r2(n: &U1024) -> U1024 {
        let r_mod_n = U1024::ZERO.borrowing_sub(n).0.div_rem(n).1;
        Self::mod_mul_raw(&r_mod_n, &r_mod_n, n)
    }

    fn mod_mul_raw(a: &U1024, b: &U1024, n: &U1024) -> U1024 {
        let (lo, hi) = a.widening_mul(b);
        if hi == U1024::ZERO {
            return lo.div_rem(n).1;
        }
        let r = U1024::ZERO.borrowing_sub(n).0.div_rem(n).1;
        let hi_reduced = Self::mod_mul_raw(&r, &hi, n);
        let lo_reduced = lo.div_rem(n).1;
        let (sum, carry) = lo_reduced.carrying_add(&hi_reduced);
        let (sub_res, borrow) = sum.borrowing_sub(n);
        if carry || !borrow { sub_res } else { sum }
    }
}
