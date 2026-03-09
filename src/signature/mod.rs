pub mod ecdsa;
pub mod schnorr;

use sha2::{Digest, Sha256};

use crate::U1024;

pub fn hash_message(message: &[u8]) -> U1024 {
    let mut buf = [0u8; 128];
    for i in 0..4 {
        let mut hasher = Sha256::new();
        hasher.update([i as u8]);
        hasher.update(message);
        let hash_result = hasher.finalize();
        buf[i * 32..(i + 1) * 32].copy_from_slice(&hash_result);
    }
    U1024::from_be_bytes(&buf)
}
