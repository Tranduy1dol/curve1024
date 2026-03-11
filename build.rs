use std::{env, fs, path::Path};

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use serde::Deserialize;

const CONFIG_PATH: &str = "config/curve1024.toml";

#[derive(Deserialize)]
struct CurveConfig {
    field: FieldConfig,
    #[allow(dead_code)]
    cm: Option<CmConfig>,
    #[allow(dead_code)]
    curve: Option<CurveEqConfig>,
    #[allow(dead_code)]
    generator: Option<GeneratorConfig>,
}

#[derive(Deserialize)]
struct FieldConfig {
    modulus: String,
    order: String,
}

#[derive(Deserialize)]
struct CmConfig {
    #[allow(dead_code)]
    trace: String,
    #[allow(dead_code)]
    y: String,
    #[allow(dead_code)]
    embedding_degree: u64,
    #[allow(dead_code)]
    discriminant: String,
}

#[derive(Deserialize)]
struct CurveEqConfig {
    #[allow(dead_code)]
    a: String,
    #[allow(dead_code)]
    b: String,
}

#[derive(Deserialize)]
struct GeneratorConfig {
    #[allow(dead_code)]
    x: String,
    #[allow(dead_code)]
    y: String,
}

impl Default for CurveConfig {
    fn default() -> Self {
        Self {
            field: FieldConfig {
                modulus: "0x0".to_string(),
                order: "0x0".to_string(),
            },
            cm: None,
            curve: None,
            generator: None,
        }
    }
}

struct U1024(pub [u64; 16]);

impl U1024 {
    fn from_hex(hex: &str) -> Self {
        let hex = hex.trim_start_matches("0x");
        let mut limbs = [0u64; 16];
        let mut limb_idx = 0;
        let mut char_idx = hex.len();

        while char_idx > 0 && limb_idx < 16 {
            let start = char_idx.saturating_sub(16);
            let chunk = &hex[start..char_idx];
            limbs[limb_idx] =
                u64::from_str_radix(chunk, 16).expect("Invalid hex in curve1024.toml");
            limb_idx += 1;
            char_idx = start;
        }
        Self(limbs)
    }
}

impl std::fmt::Display for U1024 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let inner: Vec<String> = self.0.iter().map(|l| format!("{l:#018x}")).collect();
        write!(f, "U1024([{}])", inner.join(", "))
    }
}

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest = Path::new(&out_dir).join("constants.rs");

    let config = if Path::new(CONFIG_PATH).exists() {
        let content = fs::read_to_string(CONFIG_PATH).expect("Failed to read curve1024.toml");
        let cfg: CurveConfig = toml::from_str(&content).expect("Failed to parse curve1024.toml");
        println!("cargo::rerun-if-changed=curve1024.toml");
        cfg
    } else {
        println!("cargo:warning=curve1024.toml not found, using zero placeholders");
        CurveConfig::default()
    };

    let field_cfg = config.field;
    let modulus = U1024::from_hex(&field_cfg.modulus);
    let order = U1024::from_hex(&field_cfg.order);

    let mod_biguint =
        BigUint::parse_bytes(field_cfg.modulus.trim_start_matches("0x").as_bytes(), 16)
            .unwrap_or(BigUint::ZERO);
    let r_biguint: BigUint = BigUint::from(1u32) << 1024usize;
    let r2_biguint = (&r_biguint * &r_biguint) % &mod_biguint;

    let p_bigint = BigInt::from_biguint(Sign::Plus, mod_biguint.clone());
    let r_bigint = BigInt::from_biguint(Sign::Plus, r_biguint.clone());

    let ext_gcd = p_bigint.extended_gcd(&r_bigint);
    let mut n_prime_bigint = (-ext_gcd.x) % &r_bigint;
    if n_prime_bigint < BigInt::from(0) {
        n_prime_bigint += &r_bigint;
    }
    let n_prime_biguint = n_prime_bigint.to_biguint().unwrap();

    let curve_a = config.curve.as_ref().map(|c| c.a.as_str()).unwrap_or("0x0");
    let curve_b = config.curve.as_ref().map(|c| c.b.as_str()).unwrap_or("0x0");
    let gen_x = config
        .generator
        .as_ref()
        .map(|g| g.x.as_str())
        .unwrap_or("0x0");
    let gen_y = config
        .generator
        .as_ref()
        .map(|g| g.y.as_str())
        .unwrap_or("0x0");

    let a_biguint = BigUint::parse_bytes(curve_a.trim_start_matches("0x").as_bytes(), 16)
        .unwrap_or(BigUint::ZERO);
    let a_mont = (&a_biguint * &r_biguint) % &mod_biguint;

    let b_biguint = BigUint::parse_bytes(curve_b.trim_start_matches("0x").as_bytes(), 16)
        .unwrap_or(BigUint::ZERO);
    let b_mont = (&b_biguint * &r_biguint) % &mod_biguint;

    let gx_biguint = BigUint::parse_bytes(gen_x.trim_start_matches("0x").as_bytes(), 16)
        .unwrap_or(BigUint::ZERO);
    let gx_mont = (&gx_biguint * &r_biguint) % &mod_biguint;

    let gy_biguint = BigUint::parse_bytes(gen_y.trim_start_matches("0x").as_bytes(), 16)
        .unwrap_or(BigUint::ZERO);
    let gy_mont = (&gy_biguint * &r_biguint) % &mod_biguint;

    let a = U1024::from_hex(&format!("{a_mont:x}"));
    let b = U1024::from_hex(&format!("{b_mont:x}"));
    let gx = U1024::from_hex(&format!("{gx_mont:x}"));
    let gy = U1024::from_hex(&format!("{gy_mont:x}"));
    let r2 = U1024::from_hex(&format!("{r2_biguint:x}"));
    let n_prime = U1024::from_hex(&format!("{n_prime_biguint:x}"));

    // Compute R2 and N_PRIME for the scalar field (curve order)
    let order_biguint =
        BigUint::parse_bytes(field_cfg.order.trim_start_matches("0x").as_bytes(), 16)
            .unwrap_or(BigUint::ZERO);
    let order_r2_biguint = (&r_biguint * &r_biguint) % &order_biguint;
    let order_bigint = BigInt::from_biguint(Sign::Plus, order_biguint.clone());
    let ext_gcd_order = order_bigint.extended_gcd(&r_bigint);
    let mut order_n_prime_bigint = (-ext_gcd_order.x) % &r_bigint;
    if order_n_prime_bigint < BigInt::from(0) {
        order_n_prime_bigint += &r_bigint;
    }
    let order_n_prime_biguint = order_n_prime_bigint.to_biguint().unwrap();

    let order_r2 = U1024::from_hex(&format!("{order_r2_biguint:x}"));
    let order_n_prime = U1024::from_hex(&format!("{order_n_prime_biguint:x}"));

    let code = format!(
        r#"// Auto-generated from curve1024.toml — do not edit manually

pub const CURVE_MODULUS: U1024 = {modulus};
pub const CURVE_ORDER: U1024 = {order};
pub const CURVE_A: U1024 = {a};
pub const CURVE_B: U1024 = {b};
pub const CURVE_GENERATOR_X: U1024 = {gx};
pub const CURVE_GENERATOR_Y: U1024 = {gy};

// Computed from CURVE_MODULUS
// R2 = R^2 mod p  where R = 2^1024
// N_PRIME = -p^(-1) mod 2^1024
pub const CURVE_R2: U1024 = {r2};
pub const CURVE_N_PRIME: U1024 = {n_prime};

// Computed from CURVE_ORDER (for scalar field Montgomery arithmetic)
pub const CURVE_ORDER_R2: U1024 = {order_r2};
pub const CURVE_ORDER_N_PRIME: U1024 = {order_n_prime};
"#,
        modulus = modulus,
        order = order,
        a = a,
        b = b,
        gx = gx,
        gy = gy,
        r2 = r2,
        n_prime = n_prime,
        order_r2 = order_r2,
        order_n_prime = order_n_prime,
    );

    fs::write(&dest, code).expect("Failed to write constants.rs");
}
