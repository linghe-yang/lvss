use crate::oneshot::{decrypt, encrypt, ring_lwe_dec, ring_lwe_enc, setup, verify};
use crate::r_ring::R;
use crate::util::concat_vectors;
use nalgebra::{DMatrix, DVector};

mod gadget;
mod oneshot;
pub mod r_ring;
pub mod rp;
pub mod rq;
pub mod util;
mod zetas;
mod primefind;
mod barrett;

pub use gadget::build_b_matrix;
pub use gadget::calculate_u;
pub use gadget::split_m_bar;
pub use gadget::gadget_decompose;
pub use gadget::gadget_reconstruct;
pub use util::*;
pub use oneshot::SecretKey;
pub use oneshot::PublicKey;

#[derive(Debug,Clone)]
pub struct Store {
    r: DVector<R>,
    e: DVector<R>,
    e_prime: DVector<R>,

    m: DVector<R>,
}

/// Verifiable Encryption (VE) implementation based on the one-shot scheme from the paper
/// "One-Shot Verifiable Encryption from Lattices" by Lyubashevsky and Neven (2017).
/// This provides a relaxed verifiable encryption for relations of the form Bm = u mod p,
/// where m is a short-norm vector.
/// we further adapts it for custom relations v = R * x + y mod p via gadget decomposition.
pub struct VE;

impl VE {

    /// Generates a key pair for the verifiable encryption scheme.
    pub fn gen_keypair() -> (PublicKey, SecretKey) {
        setup()
    }

    /// Encrypts vectors x and y under the public key using Ring-LWE encryption after gadget decomposition.
    /// This adapts the encryption to handle large-norm vectors by decomposing them into small-norm components
    /// (as per the gadget technique mentioned in the paper for handling relations like Equation (1)),
    /// concatenates them into m, and encrypts m. Returns the ciphertext (v, w) and a Store for later proof generation.
    pub fn encrypt(
        pk: &PublicKey,
        x: &DVector<R>,
        y: &DVector<R>,
    ) -> (DVector<R>, DVector<R>, Store) {
        let s_x = gadget_decompose(x);
        let s_y = gadget_decompose(y);
        let m = concat_vectors(&[s_x, s_y]);
        let (v, w, r, e, e_prime) = ring_lwe_enc(pk, &m);
        let store = Store { r, e, e_prime, m };
        (v, w, store)
    }

    /// Generates a zero-knowledge proof of knowledge for the encrypted witness.
    /// Uses the stored randomness and errors from encryption to produce a Fiat-Shamir with Aborts proof (c, z)
    /// showing that the ciphertext encrypts a witness m satisfying Bm = u mod p, where B and u are provided.
    pub fn prove(
        pk: &PublicKey,
        store: &Store,
        b: &DMatrix<R>,
        u: &DVector<R>,
        v: &DVector<R>,
        w: &DVector<R>,
    ) -> (R, DVector<R>) {
        encrypt(pk, b, u, &store.m, &store.r, &store.e, &store.e_prime, v, w)
    }

    /// Verifies the proof for the verifiable encryption.
    /// Checks if the proof (c, z) is valid for the ciphertext (v, w) and relation Bm = u mod p under the public key.
    pub fn verify(
        pk: &PublicKey,
        b: &DMatrix<R>,
        u: &DVector<R>,
        t: (&DVector<R>, &DVector<R>, &R, &DVector<R>),
    ) -> bool {
        verify(pk, b, u, t)
    }

    /// Decrypts the ciphertext using the secret key and the proof.
    /// This recovers a witness m' (possibly different from the original m but satisfying
    /// the relaxed relation Bm' = c'u mod p as per Section 1.3). Leveraging Lemma 2.2,
    /// when p ≡ 5 mod 8, small-norm elements are invertible in Rp, enabling deterministic
    /// recovery of the original m via m = m' * (c / c') mod p. Returns Some(m) on success
    /// or None on failure.
    pub fn try_decrypt(
        pk: &PublicKey,
        sk: &SecretKey,
        t: (&DVector<R>, &DVector<R>, &R, &DVector<R>),
    ) -> Option<DVector<R>> {
        decrypt(pk, sk, t)
    }

    pub fn decrypt(
        sk: &SecretKey,
        v: &DVector<R>,
        w: &DVector<R>,
    ) -> DVector<R> {
        ring_lwe_dec(sk,v,w)
    }
}

#[cfg(test)]
mod tests {
    use crate::r_ring::R;
    use crate::rp::{LOG_P, P_PRIME};
    use crate::{VE, build_b_matrix, calculate_u, split_m_bar, gadget_reconstruct};
    use nalgebra::DMatrix;
    use crate::util::random_dvector;

    #[test]
    fn test() {
        let (pk, sk) = VE::gen_keypair();
        let xl = 2;
        let yl = 1;
        let half_p = P_PRIME / 2;
        let x = random_dvector(xl);
        let y = random_dvector(yl);
        let (v, w, st) = VE::encrypt(&pk, &x, &y);

        let r = DMatrix::from_fn(yl, xl, |_, _| R::random_uniform(-half_p, half_p));
        let u = calculate_u(&r, &x, &y);
        let b = build_b_matrix(&r);

        let (c, z) = VE::prove(&pk, &st, &b, &u, &v, &w);

        let t = (&v, &w, &c, &z);

        let vs = VE::verify(&pk, &b, &u, t);
        assert!(vs, "Invalid CipherText");
        let m_bar = VE::try_decrypt(&pk, &sk, t);
        match m_bar {
            Some(m_bar) => {
                assert_eq!(m_bar, st.m, "Inconsistent decryption");
                assert_eq!(m_bar.len(), xl * LOG_P + yl * LOG_P);

                let (sx,sy) = split_m_bar(m_bar, xl, yl).unwrap();
                let x_recover = gadget_reconstruct(&sx);
                let y_recover = gadget_reconstruct(&sy);
                for (i,ele) in x_recover.iter().enumerate(){
                    assert!(ele.equal_mod_p(&x[i]))
                }
                for (i,ele) in y_recover.iter().enumerate(){
                    assert!(ele.equal_mod_p(&y[i]))
                }

                println!("Decryption Success!");

            }
            None => {
                panic!("Decrypt failed");
            }
        }
    }
}
