use nalgebra::{DMatrix, DVector};
use rand::{Rng, RngCore, SeedableRng};
use rand::prelude::StdRng;
use rand::seq::SliceRandom;
use rand_distr::num_traits::Zero;
use sha2::{Digest, Sha256};
use crate::oneshot::{PublicKey, MAX_C_L1};
use crate::r_ring::{N, R};


pub fn vector_scalar_mul_mod_q(vec: &DVector<R>, scalar: R) -> DVector<R> {
    let len = vec.len();
    let mut result = DVector::zeros(len);
    for i in 0..len {
        result[i] = vec[i].scalar_mul_mod_q(scalar.coeffs[0]); // Assuming scalar is constant polynomial
    }
    result
}

pub fn vector_scalar_mul_mod_p(vec: &DVector<R>, scalar: R) -> DVector<R> {
    let len = vec.len();
    let mut result = DVector::zeros(len);
    for i in 0..len {
        result[i] = vec[i].scalar_mul_mod_p(scalar.coeffs[0]); // Assuming scalar is constant polynomial
    }
    result
}

pub fn matrix_vector_mul_mod_q(mat: &DMatrix<R>, vec: &DVector<R>) -> DVector<R> {
    let rows = mat.nrows();
    let mut result = DVector::zeros(rows);
    for i in 0..rows {
        let mut sum = R::zero();
        for j in 0..mat.ncols() {
            sum += mat[(i, j)].mul_mod_q(&vec[j]);
        }
        sum.mod_q();
        result[i] = sum;
    }
    result
}

pub fn matrix_vector_mul_mod_p(mat: &DMatrix<R>, vec: &DVector<R>) -> DVector<R> {
    assert_eq!(mat.ncols(), vec.len());
    let rows = mat.nrows();
    let mut result = DVector::zeros(rows);
    for i in 0..rows {
        let mut sum = R::zero();
        for j in 0..mat.ncols() {
            sum = sum.add_mod_p(&mat[(i, j)].mul_mod_p(&vec[j]));
        }
        // sum.mod_p();
        result[i] = sum;
    }
    result
}

pub fn matrix_vector_mul(mat: &DMatrix<R>, vec: &DVector<R>) -> DVector<R> {
    assert_eq!(mat.ncols(), vec.len());
    let rows = mat.nrows();
    let mut result = DVector::zeros(rows);
    for i in 0..rows {
        let mut sum = R::zero();
        for j in 0..mat.ncols() {
            sum += mat[(i, j)] * vec[j];
        }
        result[i] = sum;
    }
    result
}

pub fn matrix_matrix_mul_mod_p(mat1: &DMatrix<R>, mat2: &DMatrix<R>) -> DMatrix<R> {
    let (rows1, cols1) = mat1.shape();
    let (cols2, _) = mat2.shape();
    assert_eq!(cols1, cols2, "Matrix dimensions must be compatible for multiplication");

    let mut result = DMatrix::zeros(rows1, cols2);
    for i in 0..rows1 {
        for j in 0..cols2 {
            let mut sum = R::zero();
            for k in 0..cols1 {
                sum = sum.add_mod_p(&mat1[(i, k)].mul_mod_p(&mat2[(k, j)]));
            }
            result[(i, j)] = sum;
        }
    }
    result
}

pub fn matrix_matrix_mul_mod_q(mat1: &DMatrix<R>, mat2: &DMatrix<R>) -> DMatrix<R> {
    let (rows1, cols1) = mat1.shape();
    let (cols2, _) = mat2.shape();
    assert_eq!(cols1, cols2, "Matrix dimensions must be compatible for multiplication");

    let mut result = DMatrix::zeros(rows1, cols2);
    for i in 0..rows1 {
        for j in 0..cols2 {
            let mut sum = R::zero();
            for k in 0..cols1 {
                sum = sum.add_mod_q(&mat1[(i, k)].mul_mod_q(&mat2[(k, j)]));
            }
            result[(i, j)] = sum;
        }
    }
    result
}

// Helper: Concat multiple DVector<R>
pub fn concat_vectors(vecs: &[DVector<R>]) -> DVector<R> {
    let mut combined = Vec::new();
    for v in vecs {
        combined.extend_from_slice(v.as_slice());
    }
    DVector::from_vec(combined)
}

// Helper: L2 norm squared over all coefficients in vector of R
pub fn vector_norm_l2_squared(vec: &DVector<R>) -> f64 {
    vec.iter().map(|r| r.norm_l2().powi(2)).sum()
}

pub fn euclidean_norm(vec: &DVector<R>) -> f64 {
    vector_norm_l2_squared(vec).sqrt()
}

// // Helper: Dot product over all coefficients (flattened)
pub fn vector_dot_product(vec1: &DVector<R>, vec2: &DVector<R>) -> f64 {
    let mut dot = 0.0;
    for (r1, r2) in vec1.iter().zip(vec2.iter()) {
        for (c1, c2) in r1.coeffs.iter().zip(r2.coeffs.iter()) {
            dot += (*c1 as f64) * (*c2 as f64);
        }
    }
    dot
}

pub fn vector_norm_inf(vec: &DVector<R>) -> i64 {
    vec.iter().map(|r| r.norm_inf()).max().unwrap_or(0)
}
pub fn vector_add(vec1: &DVector<R>, vec2: &DVector<R>) -> DVector<R> {
    assert_eq!(vec1.len(), vec2.len(), "Vector lengths must match");
    DVector::from_fn(vec1.len(), |i, _| vec1[i] + vec2[i])
}

pub fn vector_add_mod_q(vec1: &DVector<R>, vec2: &DVector<R>) -> DVector<R> {
    assert_eq!(vec1.len(), vec2.len(), "Vector lengths must match");
    DVector::from_fn(vec1.len(), |i, _| vec1[i].add_mod_q(&vec2[i]))
}
pub fn vector_add_mod_p(vec1: &DVector<R>, vec2: &DVector<R>) -> DVector<R> {
    assert_eq!(vec1.len(), vec2.len(), "Vector lengths must match");
    DVector::from_fn(vec1.len(), |i, _| vec1[i].add_mod_p(&vec2[i]))
}

pub fn random_dvector(len: usize) -> DVector<R> {
    DVector::from_iterator(
        len,
        (0..len).map(|_| R::random_in_p()),
    )
}

pub fn random_gaussian_dvector(len: usize, sigma: f64) -> DVector<R> {
    let rng = &mut rand::thread_rng();
    DVector::from_iterator(
        len,
        (0..len).map(|_| R::random_gaussian(rng, sigma)),
    )
}
pub fn generate_low_l1_r(rng: &mut impl RngCore) -> R {
    let mut coeffs = [0i64; N];
    let max_non_zero = MAX_C_L1.min(N);

    // Randomly choose number of non-zero coefficients
    let num_non_zero = rng.gen_range(0..=max_non_zero);

    // Select num_non_zero random indices without replacement
    let mut indices: Vec<usize> = (0..N).collect();
    indices.shuffle(rng);

    // Set coefficients to ±1
    for &idx in indices.iter().take(num_non_zero) {
        coeffs[idx] = if rng.gen_bool(0.5) { 1 } else { -1 };
    }

    R { coeffs }
}

pub fn h(pk: &PublicKey, b: &DMatrix<R>, target: &DVector<R>, commitment: &DVector<R>) -> R {

    let mut hasher = Sha256::new();

    hasher.update(pk.a.to_bytes());
    hasher.update(pk.t.to_bytes());
    hasher.update(pk.p.to_le_bytes());
    hasher.update(pk.q.to_le_bytes());


    hasher.update((b.nrows() as u32).to_le_bytes());
    hasher.update((b.ncols() as u32).to_le_bytes());
    for i in 0..b.nrows() {
        for j in 0..b.ncols() {
            hasher.update(b[(i, j)].to_bytes());
        }
    }

    hasher.update((target.len() as u32).to_le_bytes());
    for i in 0..target.len() {
        hasher.update(target[i].to_bytes());
    }


    // hasher.update((commitment.len() as u32).to_le_bytes());
    for i in 0..commitment.len() {
        hasher.update(commitment[i].to_bytes());
    }

    let hash = hasher.finalize();

    // Deterministically seed RNG from hash
    let seed: [u8; 32] = hash.into();
    let mut rng = StdRng::from_seed(seed);

    // Generate c using the seeded RNG
    generate_low_l1_r(&mut rng)
}



pub fn poly_extended_euclid(a: Vec<i64>, b: Vec<i64>, p: i64) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    if b.iter().all(|&x| x == 0) {
        (trim_trailing_zeros(a), vec![1], vec![0])
    } else {
        let (q, r) = poly_div(&a, &b, p);
        let (g, s, t) = poly_extended_euclid(b, r, p);
        let new_t = poly_sub(&s, &poly_mul(&t, &q, p), p);
        (g, t, new_t)
    }
}
fn mod_inverse(a: i64, p: i64) -> Option<i64> {
    let (g, s, _) = extended_euclid_int(a, p);
    if g != 1 {
        None
    } else {
        Some((s % p + p) % p)
    }
}

fn extended_euclid_int(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_euclid_int(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

// Helper for polynomial operations (Vec<i64> low to high degree)
fn trim_trailing_zeros(mut v: Vec<i64>) -> Vec<i64> {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    if v.is_empty() {
        v.push(0);
    }
    v
}

// fn poly_add(a: &Vec<i64>, b: &Vec<i64>, p: i64) -> Vec<i64> {
//     let max_len = a.len().max(b.len());
//     let mut result = vec![0; max_len];
//     for i in 0..max_len {
//         let aa = if i < a.len() { a[i] } else { 0 };
//         let bb = if i < b.len() { b[i] } else { 0 };
//         result[i] = ((aa + bb) % p + p) % p;
//     }
//     trim_trailing_zeros(result)
// }

fn poly_sub(a: &Vec<i64>, b: &Vec<i64>, p: i64) -> Vec<i64> {
    let max_len = a.len().max(b.len());
    let mut result = vec![0; max_len];
    for i in 0..max_len {
        let aa = if i < a.len() { a[i] } else { 0 };
        let bb = if i < b.len() { b[i] } else { 0 };
        result[i] = ((aa - bb) % p + p) % p;
    }
    trim_trailing_zeros(result)
}

fn poly_mul(a: &Vec<i64>, b: &Vec<i64>, p: i64) -> Vec<i64> {
    let len = a.len() + b.len() - 1;
    let mut result = vec![0; len];
    for i in 0..a.len() {
        for j in 0..b.len() {
            result[i + j] = ((result[i + j] + a[i] * b[j]) % p + p) % p;
        }
    }
    trim_trailing_zeros(result)
}

fn poly_div(dividend: &Vec<i64>, divisor: &Vec<i64>, p: i64) -> (Vec<i64>, Vec<i64>) {
    let mut dividend = dividend.clone();
    let divisor = divisor.clone();
    if divisor.is_empty() || divisor.last().unwrap_or(&0) == &0 {
        panic!("Divisor cannot be zero");
    }
    if dividend.len() < divisor.len() {
        return (vec![0], dividend);
    }
    let mut quotient = vec![0; dividend.len() - divisor.len() + 1];
    let lead_divisor_inv = mod_inverse(*divisor.last().unwrap(), p).unwrap();
    for i in (0..quotient.len()).rev() {
        let lead_dividend = if i + divisor.len() - 1 < dividend.len() {
            dividend[i + divisor.len() - 1]
        } else {
            0
        };
        let q_term = ((lead_dividend * lead_divisor_inv) % p + p) % p;
        quotient[i] = q_term;
        for j in 0..divisor.len() {
            if i + j < dividend.len() {
                dividend[i + j] = ((dividend[i + j] - q_term * divisor[j]) % p + p) % p;
            }
        }
    }
    (trim_trailing_zeros(quotient), trim_trailing_zeros(dividend))
}
