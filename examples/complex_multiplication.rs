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
        .expect("Run cock_pinch example first to generate config");

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

fn point_add(p: &Point, q: &Point, coeff_a: &U1024, f: &RuntimeField) -> Point {
    if p.is_infinity {
        return q.clone();
    }
    if q.is_infinity {
        return p.clone();
    }
    if p.x == q.x {
        if p.y == q.y && p.y != U1024::ZERO {
            return point_double(p, coeff_a, f);
        }
        return Point::infinity();
    }

    let num = f.sub(&q.y, &p.y);
    let den_inv = f.inv(&f.sub(&q.x, &p.x));
    let lambda = f.mul(&num, &den_inv);

    let x3 = f.sub(&f.sub(&f.mul(&lambda, &lambda), &p.x), &q.x);
    let y3 = f.sub(&f.mul(&lambda, &f.sub(&p.x, &x3)), &p.y);

    Point::new(x3, y3)
}

fn point_double(p: &Point, coeff_a: &U1024, f: &RuntimeField) -> Point {
    if p.is_infinity || p.y == U1024::ZERO {
        return Point::infinity();
    }

    // λ = (3x² + a) / (2y)
    let x_sq = f.mul(&p.x, &p.x);
    let num = f.add(&f.mul(&U1024::from(3), &x_sq), coeff_a);
    let den_inv = f.inv(&f.add(&p.y, &p.y));
    let lambda = f.mul(&num, &den_inv);

    let x3 = f.sub(&f.mul(&lambda, &lambda), &f.add(&p.x, &p.x));
    let y3 = f.sub(&f.mul(&lambda, &f.sub(&p.x, &x3)), &p.y);

    Point::new(x3, y3)
}

fn scalar_mul(p: &Point, scalar: &U1024, coeff_a: &U1024, f: &RuntimeField) -> Point {
    let mut result = Point::infinity();
    let mut base = p.clone();

    for i in 0..1024 {
        if (scalar.0[i / 64] >> (i % 64)) & 1 == 1 {
            result = point_add(&result, &base, coeff_a, f);
        }
        base = point_double(&base, coeff_a, f);
    }

    result
}

fn sqrt_mod(a: &U1024, f: &RuntimeField) -> Option<U1024> {
    let p = &f.modulus;
    let p_minus_1 = p.borrowing_sub(&U1024::ONE).0;
    let half_p_minus_1 = p_minus_1.shr(1);

    if f.pow(a, &half_p_minus_1) != U1024::ONE {
        return None;
    }

    let mut q = p_minus_1;
    let mut s = 0u32;
    while q.0[0] & 1 == 0 {
        q = q.shr(1);
        s += 1;
    }

    let mut z = U1024::from(2);
    while f.pow(&z, &half_p_minus_1) != p_minus_1 {
        z = z.carrying_add(&U1024::ONE).0;
    }

    let mut m = s;
    let mut c = f.pow(&z, &q);
    let mut t = f.pow(a, &q);
    let q_plus_1_half = q.carrying_add(&U1024::ONE).0.shr(1);
    let mut root = f.pow(a, &q_plus_1_half);

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

fn hilbert_class_poly_root(d: &U1024, f: &RuntimeField) -> U1024 {
    // For |D|=3 (D=-3): H_{-3}(x) = x, so j = 0
    // For |D|=4 (D=-4): H_{-4}(x) = x - 1728, so j = 1728
    // General D requires precomputed Hilbert class polynomials
    if *d == U1024::from(3) {
        U1024::ZERO
    } else if *d == U1024::from(4) {
        f.reduce(&U1024::from(1728))
    } else {
        panic!("Hilbert class polynomial not implemented for |D|={d}");
    }
}

fn compute_curve_coefficients(j: &U1024, f: &RuntimeField) -> (U1024, U1024) {
    if *j == U1024::ZERO {
        // j=0 → y² = x³ + b, start with b=1
        (U1024::ZERO, U1024::ONE)
    } else if *j == U1024::from(1728) {
        // j=1728 → y² = x³ + ax, start with a=1
        (U1024::ONE, U1024::ZERO)
    } else {
        // General: c = j / (1728 - j), a = 3c, b = 2c
        let c = f.mul(j, &f.inv(&f.sub(&U1024::from(1728), j)));
        let a = f.mul(&U1024::from(3), &c);
        let b = f.mul(&U1024::from(2), &c);
        (a, b)
    }
}

fn find_point_on_curve(a: &U1024, b: &U1024, f: &RuntimeField) -> Point {
    let mut x = U1024::ONE;
    loop {
        // rhs = x³ + ax + b
        let x_sq = f.mul(&x, &x);
        let x_cubed = f.mul(&x_sq, &x);
        let ax = f.mul(a, &x);
        let rhs = f.add(&f.add(&x_cubed, &ax), b);

        if let Some(y) = sqrt_mod(&rhs, f) {
            return Point::new(x, y);
        }

        x = x.carrying_add(&U1024::ONE).0;
    }
}

// For |D| = 3 (j=0), we can just try small values of b.
// We need #E = p + 1 - t
// For |D| = 3 (j=0), we can just try small values of b.
// We need the group order to be a multiple of r. There are 6 twist orders.
fn twist_test(a: &U1024, config: &CurveConfig, f: &RuntimeField) -> (U1024, U1024, U1024) {
    println!("  Trying small values of b since |D|=3 (j=0) has 6 twists...");

    let p_plus_1 = config.p.carrying_add(&U1024::ONE).0;

    // The 6 possible traces are ±t, ±(t+3y)/2, ±(t-3y)/2
    let three_y = config.y.widening_mul(&U1024::from(3)).0;

    let trace_3_plus = config.t.carrying_add(&three_y).0.shr(1);
    let trace_3_minus = if config.t > three_y {
        config.t.borrowing_sub(&three_y).0.shr(1)
    } else {
        three_y.borrowing_sub(&config.t).0.shr(1)
    };

    let candidates = [
        p_plus_1.borrowing_sub(&config.t).0,
        p_plus_1.carrying_add(&config.t).0,
        p_plus_1.borrowing_sub(&trace_3_plus).0,
        p_plus_1.carrying_add(&trace_3_plus).0,
        p_plus_1.borrowing_sub(&trace_3_minus).0,
        p_plus_1.carrying_add(&trace_3_minus).0,
    ];

    let mut correct_order = U1024::ZERO;

    for (i, &order) in candidates.iter().enumerate() {
        let (cofactor, rem) = order.div_rem(&config.r);
        if rem == U1024::ZERO && cofactor != U1024::ZERO {
            println!("  Found that twist order index {} is divisible by r!", i);
            correct_order = order;
            break;
        }
    }

    if correct_order == U1024::ZERO {
        panic!("None of the 6 formal CM twists have an order divisible by r!");
    }

    for b_val in 1u64..=2000 {
        let b = f.reduce(&U1024::from(b_val));

        let p_pt = find_point_on_curve(a, &b, f);
        let result = scalar_mul(&p_pt, &correct_order, a, f);

        if result.is_infinity {
            println!("  Found correct curve with b = {b_val}");
            let cofactor = correct_order.div_rem(&config.r).0;
            return (*a, b, cofactor);
        }
    }

    panic!("Could not find a valid b twist value");
}

fn find_generator(
    a: &U1024,
    b: &U1024,
    cofactor: &U1024,
    order_r: &U1024,
    f: &RuntimeField,
) -> Point {
    loop {
        let p = find_point_on_curve(a, b, f);
        let g = scalar_mul(&p, cofactor, a, f);

        if g.is_infinity {
            continue;
        }

        // Verify [r]G = O
        let check = scalar_mul(&g, order_r, a, f);
        if check.is_infinity {
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
    println!("=== CM Curve Equation & Generator Finder ===\n");

    let config = read_config();
    println!("Loaded config from {CONFIG_PATH}");
    println!("  p: {} bits", bit_length(&config.p));
    println!("  r: {} bits", bit_length(&config.r));
    println!("  |D| = {}\n", config.d);

    let field = RuntimeField::new(config.p);

    // Step 1: Hilbert class polynomial → j-invariant
    println!("Step 1: Computing j-invariant from Hilbert class polynomial...");
    let j = hilbert_class_poly_root(&config.d, &field);
    println!("  j = {j}\n");

    // Step 2: Curve coefficients from j-invariant
    println!("Step 2: Computing initial curve coefficients...");
    let (a_init, b_init) = compute_curve_coefficients(&j, &field);
    println!("  a = {a_init}");
    println!("  b = {b_init}\n");

    // Step 3: Twist test
    println!("Step 3: Twist test...");
    let (a, b, cofactor) = twist_test(&a_init, &config, &field);
    println!("  Final: y² = x³ + {a}·x + {b}\n");

    // Step 4: Cofactor and generator
    println!("Step 4: Finding generator point...");
    println!("  cofactor = {cofactor}");
    println!("  cofactor bits = {}", bit_length(&cofactor));

    let timer = std::time::Instant::now();
    let g = find_generator(&a, &b, &cofactor, &config.r, &field);
    println!("  Generator found in {:.2?}", timer.elapsed());
    println!("  G.x = {}", g.x);
    println!("  G.y = {}", g.y);

    // Verify
    println!("\nStep 5: Verification...");
    let rg = scalar_mul(&g, &config.r, &a, &field);
    assert!(rg.is_infinity, "[r]G must be point at infinity");
    println!("  ✓ [r]G = O (generator has order r)");

    // Verify point is on curve
    let lhs = field.mul(&g.y, &g.y);
    let x3 = field.mul(&field.mul(&g.x, &g.x), &g.x);
    let ax = field.mul(&a, &g.x);
    let rhs = field.add(&field.add(&x3, &ax), &b);
    assert!(lhs == rhs, "Generator must be on curve");
    println!("  ✓ G is on curve y² = x³ + ax + b");

    println!("\nWriting parameters to {}...", CONFIG_PATH);
    append_to_config(&a, &b, &g.x, &g.y);

    println!("=== Done ===");
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
