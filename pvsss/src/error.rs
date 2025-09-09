#[derive(Debug)]
pub enum DecryptionError {
    InvalidLength,
    InvalidCipher
}