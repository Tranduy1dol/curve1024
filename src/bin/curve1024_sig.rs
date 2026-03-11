use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

use clap::{Parser, Subcommand, ValueEnum};

use curve1024::{
    AffinePoint, Curve1024BaseField, Curve1024Config, EcdsaSignature, KeyPair, PrimeFieldElement,
    SchnorrSignature, U1024,
};

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

        #[arg(long, default_value = "schnorr")]
        scheme: Scheme,
    },

    Verify {
        #[arg(short, long)]
        file: PathBuf,

        #[arg(short, long, default_value = "key.pub")]
        key: PathBuf,
    },

    ExportPub {
        #[arg(short, long, default_value = "key.priv")]
        key: PathBuf,

        #[arg(short, long, default_value = "key.pub")]
        output: PathBuf,
    },
}

#[derive(Clone, ValueEnum)]
enum Scheme {
    Schnorr,
    Ecdsa,
}

impl Scheme {
    fn tag(&self) -> u8 {
        match self {
            Scheme::Schnorr => 0x00,
            Scheme::Ecdsa => 0x01,
        }
    }

    fn from_tag(tag: u8) -> Self {
        match tag {
            0x00 => Scheme::Schnorr,
            0x01 => Scheme::Ecdsa,
            _ => panic!("Unknown signature scheme tag: {tag}"),
        }
    }

    fn sig_size(&self) -> usize {
        match self {
            Scheme::Schnorr => 384, // R.x(128) + R.y(128) + s(128)
            Scheme::Ecdsa => 256,   // r(128) + s(128)
        }
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Keygen { output } => {
            cmd_keygen(&output);
        }
        Commands::Sign { file, key, scheme } => {
            cmd_sign(&file, &key, &scheme);
        }
        Commands::Verify { file, key } => {
            cmd_verify(&file, &key);
        }
        Commands::ExportPub { key, output } => {
            cmd_export_pub(&key, &output);
        }
    }
}

fn cmd_keygen(output: &Path) {
    let start = Instant::now();
    let keypair = KeyPair::<Curve1024Config>::generate();

    keypair
        .save(output.to_str().expect("Output dir not valid"))
        .expect("Failed to save keypair");

    let pub_path = output.with_extension("pub");
    let pub_x = keypair.public_key.x.to_bytes();
    let pub_y = keypair.public_key.y.to_bytes();
    fs::write(pub_path.clone(), [pub_x, pub_y].concat()).expect("Failed to save keypair");

    let elapsed = start.elapsed();
    println!("Keypair generated in {elapsed:.2?}");

    println!("Private key saved to {:?}", output.to_str());
    println!("Public key saved to {:?}", pub_path.to_str())
}

fn cmd_sign(file: &Path, key: &Path, scheme: &Scheme) {
    let start = Instant::now();

    let message = fs::read(file).expect("Failed to read file");
    let keypair = KeyPair::<Curve1024Config>::load(key.to_str().expect("Key dir not valid"))
        .expect("Failed to load keypair");

    let buf = match scheme {
        Scheme::Schnorr => {
            let sig = SchnorrSignature::<Curve1024Config>::sign(&keypair.private_key, &message);
            [
                sig.r_point.x.to_bytes(),
                sig.r_point.y.to_bytes(),
                sig.s.to_be_bytes(),
            ]
            .concat()
        }
        Scheme::Ecdsa => {
            let sig = EcdsaSignature::sign::<Curve1024Config>(&keypair.private_key, &message);
            [sig.r.to_be_bytes(), sig.s.to_be_bytes()].concat()
        }
    };

    // Format: [message] [sig_bytes] [scheme_tag: 1 byte] [magic: 4 bytes]
    let mut output = message;
    output.extend_from_slice(&buf);
    output.push(scheme.tag());
    output.extend_from_slice(b"C1K\n");
    fs::write(file, output).expect("Failed to write signed file");

    println!("Signed in {:.2?}", start.elapsed());
}

fn cmd_verify(file: &Path, key: &Path) {
    let start = Instant::now();
    let data = fs::read(file).expect("Failed to read file");
    let len = data.len();

    // Check magic footer
    if len < 5 || &data[len - 4..] != b"C1K\n" {
        panic!("File is not signed");
    }

    // Read scheme from footer
    let scheme = Scheme::from_tag(data[len - 5]);
    let sig_size = scheme.sig_size();
    let footer_size = 5 + sig_size;

    if len < footer_size {
        panic!("File too small to contain a valid signature");
    }

    // Split: [message] [sig_bytes] [tag] [magic]
    let message = &data[..len - footer_size];
    let sig_bytes = &data[len - footer_size..len - 5];

    // Load public key (256 bytes: x || y)
    let pub_data = fs::read(key).expect("Failed to read public key");
    assert!(pub_data.len() == 256, "Invalid public key file");
    let pub_x =
        PrimeFieldElement::<Curve1024BaseField>::from_bytes(pub_data[..128].try_into().unwrap());
    let pub_y =
        PrimeFieldElement::<Curve1024BaseField>::from_bytes(pub_data[128..].try_into().unwrap());
    let public_key = AffinePoint::<Curve1024Config>::new(pub_x, pub_y);

    // Verify
    let valid = match scheme {
        Scheme::Schnorr => {
            let rx = PrimeFieldElement::<Curve1024BaseField>::from_bytes(
                sig_bytes[..128].try_into().unwrap(),
            );
            let ry = PrimeFieldElement::<Curve1024BaseField>::from_bytes(
                sig_bytes[128..256].try_into().unwrap(),
            );
            let s = U1024::from_be_bytes(&sig_bytes[256..384]);
            let r_point = AffinePoint::<Curve1024Config>::new(rx, ry);
            let sig = SchnorrSignature { r_point, s };
            sig.verify(&public_key, message)
        }
        Scheme::Ecdsa => {
            let r = U1024::from_be_bytes(&sig_bytes[..128]);
            let s = U1024::from_be_bytes(&sig_bytes[128..256]);
            let sig = EcdsaSignature { r, s };
            sig.verify::<Curve1024Config>(&public_key, message)
        }
    };

    println!("Verified in {:.2?}", start.elapsed());
    if valid {
        println!("✓ Signature is VALID");
    } else {
        println!("✗ Signature is INVALID");
    }
}

fn cmd_export_pub(key: &Path, output: &Path) {
    let start = Instant::now();
    let keypair = KeyPair::<Curve1024Config>::load(key.to_str().expect("Key path not valid"))
        .expect("Failed to load private key");

    let pub_x = keypair.public_key.x.to_bytes();
    let pub_y = keypair.public_key.y.to_bytes();
    fs::write(output, [pub_x, pub_y].concat()).expect("Failed to write public key");

    println!("Public key exported in {:.2?}", start.elapsed());
    println!("Public key saved to {:?}", output.to_str());
}
