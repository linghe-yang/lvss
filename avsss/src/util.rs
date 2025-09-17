use hrcrypto::hash::Hash;
use merkle_light::hash::Algorithm;
use merkle_light::proof::Proof;
use nalgebra::{DMatrix, DVector};
use ve::r_ring::R;
use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::SeedableRng;

pub(crate) fn cipher_to_bytes(tuple: (&DVector<R>, &DVector<R>)) -> Vec<u8> {
    let (vec1, vec2) = tuple;
    let mut bytes = Vec::new();
    for elem in vec1.iter() {
        bytes.extend_from_slice(&elem.to_bytes());
    }
    for elem in vec2.iter() {
        bytes.extend_from_slice(&elem.to_bytes());
    }
    bytes
}

pub fn generate_r_matrix(seed: Hash, xl: usize, yl: usize, _sigma: f64) -> DMatrix<R> {
    let mut rng = ChaCha12Rng::from_seed(seed);
    DMatrix::from_fn(yl, xl, |_, _| R::random_constant_ternary(&mut rng))
}

pub(crate) fn verify_merkle<T: Ord + Clone + AsRef<[u8]>, A: Algorithm<T>>(root:&T, leaf: T, proof: &Proof<T>, alg: &mut A ) -> bool {
    proof.item() == A::leaf(alg, leaf) && &proof.root() == root && proof.validate::<A>()
}