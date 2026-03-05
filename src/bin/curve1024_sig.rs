use std::fs;
use std::path::PathBuf;

use clap::{Parser, Subcommand};

use curve1024::{
    AffinePoint, EcdsaSignature, KeyPair, PrimeFieldConfig, PrimeFieldElement, SWCurveConfig,
    SchnorrSignature, U1024,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Curve1024Field;

impl PrimeFieldConfig for Curve1024Field {
    const MODULUS: U1024 = U1024::ZERO;
    const R2: U1024 = U1024::ZERO;
    const N_PRIME: U1024 = U1024::ZERO;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Curve1024Config;

impl SWCurveConfig for Curve1024Config {
    type BaseField = Curve1024Field;

    const COEFF_A: PrimeFieldElement<Curve1024Field> = PrimeFieldElement::ZERO;
    const COEFF_B: PrimeFieldElement<Curve1024Field> = PrimeFieldElement::ZERO;
    const ORDER: U1024 = U1024::ZERO;

    fn generator() -> AffinePoint<Self> {
        todo!("Define generator point G for Curve1024")
    }
}

#[derive(Parser)]
#[command(name = "curve1024-sig")]
#[command(about = "A GPG-like tool for elliptic curve key management and digital signatures")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Keygen {
        #[arg(short, long, default_value = "key.priv")]
        output: PathBuf,
    },

    Sign {
        #[arg(short, long)]
        file: PathBuf,

        #[arg(short, long, default_value = "key.priv")]
        key: PathBuf,

        #[arg(short, long)]
        output: Option<PathBuf>,

        #[arg(long, default_value = "schnorr")]
        scheme: String,
    },

    Verify {
        #[arg(short, long)]
        file: PathBuf,

        #[arg(short, long)]
        sig: PathBuf,

        #[arg(short, long, default_value = "key.pub")]
        key: PathBuf,

        #[arg(long, default_value = "schnorr")]
        scheme: String,
    },

    ExportPub {
        #[arg(short, long, default_value = "key.priv")]
        key: PathBuf,

        #[arg(short, long, default_value = "key.pub")]
        output: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Keygen { output } => {
            cmd_keygen(&output);
        }
        Commands::Sign {
            file,
            key,
            output,
            scheme,
        } => {
            let sig_path = output.unwrap_or_else(|| {
                let mut p = file.clone();
                p.set_extension("sig");
                p
            });
            cmd_sign(&file, &key, &sig_path, &scheme);
        }
        Commands::Verify {
            file,
            sig,
            key,
            scheme,
        } => {
            cmd_verify(&file, &sig, &key, &scheme);
        }
        Commands::ExportPub { key, output } => {
            cmd_export_pub(&key, &output);
        }
    }
}

fn cmd_keygen(output: &PathBuf) {
    println!("Generating keypair...");
    todo!("Generate keypair, save private key to {output:?} and public key to {output:?}.pub")
}

fn cmd_sign(file: &PathBuf, key: &PathBuf, sig_output: &PathBuf, scheme: &str) {
    println!("Signing {:?} with scheme '{}'...", file, scheme);
    todo!("Read file, load private key, sign, write signature to {sig_output:?}")
}

fn cmd_verify(file: &PathBuf, sig: &PathBuf, key: &PathBuf, scheme: &str) {
    println!("Verifying {:?} with scheme '{}'...", file, scheme);
    todo!("Read file, load signature and public key, verify, print result")
}

fn cmd_export_pub(key: &PathBuf, output: &PathBuf) {
    println!("Exporting public key...");
    todo!("Load private key from {key:?}, derive public key, save to {output:?}")
}
