use serde::Deserialize;
use std::{path::Path, env, fs};

#[derive(Deserialize)]
struct CurveConfig {
    modulus: String,
    order: String,
    #[allow(dead_code)]
    trace: Option<String>,
    #[allow(dead_code)]
    cm_y: Option<String>,
    #[allow(dead_code)]
    embedding_degree: Option<u64>,
    #[allow(dead_code)]
    discriminant: Option<String>,
}

fn hex_to_limbs(hex: &str) -> [u64; 16] {
    let hex = hex.trim_start_matches("0x");
    let mut limbs = [0u64; 16];
    let mut limb_idx = 0;
    let mut char_idx = hex.len();

    while char_idx > 0 && limb_idx < 16 {
        let start = char_idx.saturating_sub(16);
        let chunk = &hex[start..char_idx];
        limbs[limb_idx] = u64::from_str_radix(chunk, 16).expect("Invalid hex in curve.toml");
        limb_idx += 1;
        char_idx = start;
    }
    limbs
}

fn limbs_to_rust(limbs: &[u64; 16]) -> String {
    let inner: Vec<String> = limbs.iter().map(|l| format!("{l:#018x}")).collect();
    format!("U1024([{}])", inner.join(", "))
}

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest = Path::new(&out_dir).join("constants.rs");

    let config = if Path::new("curve.toml").exists() {
        let content = fs::read_to_string("curve.toml").expect("Failed to read curve.toml");
        let cfg: CurveConfig = toml::from_str(&content).expect("Failed to parse curve.toml");
        println!("cargo::rerun-if-changed=curve.toml");
        cfg
    } else {
        println!("cargo:warning=curve.toml not found, using zero placeholders");
        CurveConfig {
            modulus: "0x0".to_string(),
            order: "0x0".to_string(),
            trace: None,
            cm_y: None,
            embedding_degree: None,
            discriminant: None,
        }
    };

    let modulus_limbs = hex_to_limbs(&config.modulus);
    let order_limbs = hex_to_limbs(&config.order);

    let code = format!(
        r#"// Auto-generated from curve.toml — do not edit manually

pub const CURVE_MODULUS: U1024 = {modulus};
pub const CURVE_ORDER: U1024 = {order};

// TODO: These need to be computed from CURVE_MODULUS
// R2 = R^2 mod p  where R = 2^1024
// N_PRIME = -p^(-1) mod 2^1024
pub const CURVE_R2: U1024 = U1024::ZERO;
pub const CURVE_N_PRIME: U1024 = U1024::ZERO;
"#,
        modulus = limbs_to_rust(&modulus_limbs),
        order = limbs_to_rust(&order_limbs),
    );

    fs::write(&dest, code).expect("Failed to write constants.rs");
}
