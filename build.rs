use std::{env, fs, path::Path};

use serde::Deserialize;

const CONFIG_PATH: &str = "config/curve1024.toml";

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

impl Default for CurveConfig {
    fn default() -> Self {
        Self {
            modulus: "0x0".to_string(),
            order: "0x0".to_string(),
            trace: None,
            cm_y: None,
            embedding_degree: None,
            discriminant: None,
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

    fn to_string(&self) -> String {
        let inner: Vec<String> = self.0.iter().map(|l| format!("{l:#018x}")).collect();
        format!("U1024([{}])", inner.join(", "))
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

    let modulus = U1024::from_hex(&config.modulus);
    let order = U1024::from_hex(&config.order);

    let code = format!(
        r#"// Auto-generated from curve1024.toml — do not edit manually

pub const CURVE_MODULUS: U1024 = {modulus};
pub const CURVE_ORDER: U1024 = {order};

// TODO: These need to be computed from CURVE_MODULUS
// R2 = R^2 mod p  where R = 2^1024
// N_PRIME = -p^(-1) mod 2^1024
pub const CURVE_R2: U1024 = U1024::ZERO;
pub const CURVE_N_PRIME: U1024 = U1024::ZERO;
"#,
        modulus = modulus.to_string(),
        order = order.to_string(),
    );

    fs::write(&dest, code).expect("Failed to write constants.rs");
}
