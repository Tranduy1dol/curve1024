#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;

use std::fs;

use crate::{AffinePoint, SWCurveConfig, U1024};

pub struct KeyPair<C: SWCurveConfig> {
    pub private_key: U1024,
    pub public_key: AffinePoint<C>,
}

impl<C: SWCurveConfig> KeyPair<C> {
    pub fn generate() -> Self {
        let private_key = U1024::rand(&C::ORDER);
        let public_key = C::generator().mul(&private_key);
        Self {
            private_key,
            public_key,
        }
    }

    pub fn from_private_key(private_key: U1024) -> Self {
        let public_key = C::generator().mul(&private_key);
        Self {
            private_key,
            public_key,
        }
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let bytes = self.private_key.to_be_bytes();
        fs::write(path, bytes)?;

        #[cfg(unix)]
        fs::set_permissions(path, fs::Permissions::from_mode(0o600))?;

        Ok(())
    }

    pub fn load(path: &str) -> std::io::Result<Self> {
        let bytes = fs::read(path)?;
        let private_key = U1024::from_be_bytes(&bytes);
        Ok(Self::from_private_key(private_key))
    }
}
