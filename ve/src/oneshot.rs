use crate::r_ring::R;
use crate::rp::P_PRIME;
use crate::rq::Q_PRIME;
use crate::util::*;
use nalgebra::{DMatrix, DVector};
use rand::{Rng, thread_rng};
use rand_distr::num_traits::Zero;
use std::f64::consts::E;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PublicKey {
    pub a: R,
    pub t: R,
    pub p: i64,
    pub q: i64,
}
impl PublicKey {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.a.to_bytes());
        bytes.extend_from_slice(&self.t.to_bytes());
        bytes.extend_from_slice(&self.p.to_le_bytes());
        bytes.extend_from_slice(&self.q.to_le_bytes());
        bytes
    }
}
impl Default for PublicKey {
    fn default() -> Self {
        PublicKey {
            a: R::default(),
            t: R::default(),
            p: 0,
            q: 0,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecretKey {
    pub s1: R,
}

impl Default for SecretKey {
    fn default() -> Self {
        SecretKey {
            s1: R::default(),
        }
    }
}

pub const MAX_C_L1: usize = 36;

/// Setup function to generate keypair based on Ring-LWE parameters.
/// In production, use larger q (e.g., 12289 or similar 5 mod 8 prime) and adjust.
pub fn setup() -> (PublicKey, SecretKey) {
    let p: i64 = P_PRIME; // Message modulus >2
    let q: i64 = Q_PRIME; // Ring modulus, prime congruent to 5 mod 8

    // s1, s2 ← S1: uniform in {-1,0,1}
    let s1 = R::random_uniform(-1, 1);
    let s2 = R::random_uniform(-1, 1);

    // a ← Rq: uniform in [-(q-1)/2, (q-1)/2]
    let half_q = (q - 1) / 2;
    let a = R::random_uniform(-half_q, half_q);

    // t ← a * s1 + s2 mod q
    let t = a.mul_mod_q(&s1).add_mod_q(&s2);

    let pk = PublicKey { a, t, p, q };
    let sk = SecretKey { s1 };

    (pk, sk)
}

pub fn ring_lwe_enc(
    pk: &PublicKey,
    m: &DVector<R>,
) -> (DVector<R>, DVector<R>, DVector<R>, DVector<R>, DVector<R>) {
    let k = m.len();
    // Step 1: r, e, e' ← S^k_1
    let r = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    let e = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    let e_prime = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    // Step 2: v = p (a r + e) mod q, w = p (t r + e') + m mod q (vector-wise)
    let mut v = DVector::zeros(k);
    let mut w = DVector::zeros(k);
    for i in 0..k {
        v[i] = pk.a.mul_mod_q(&r[i]) + e[i];
        v[i] = v[i].scalar_mul_mod_q(pk.p);
        v[i].caddq();

        w[i] = pk.t.mul_mod_q(&r[i]) + e_prime[i];
        w[i] = w[i].scalar_mul_mod_q(pk.p).add_mod_q(&m[i]);
        w[i].caddq();
    }
    (v, w, r, e, e_prime)
}

/// Encrypt function: Returns t = (v, w, c, z)
pub fn encrypt(
    pk: &PublicKey,
    b: &DMatrix<R>,
    u: &DVector<R>,
    m: &DVector<R>,
    r: &DVector<R>,
    e: &DVector<R>,
    e_prime: &DVector<R>,
    v: &DVector<R>,
    w: &DVector<R>,
) -> (R, DVector<R>) {
    let k = m.len();
    let ell = b.nrows();
    let sigma = 50688f64;

    // // Step 1: r, e, e' ← S^k_1
    // let r = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    // let e = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    // let e_prime = DVector::from_fn(k, |_, _| R::random_uniform(-1, 1));
    //
    // // Step 2: v = p (a r + e) mod q, w = p (t r + e') + m mod q (vector-wise)
    // let mut v = DVector::zeros(k);
    // let mut w = DVector::zeros(k);
    // for i in 0..k {
    //     v[i] = pk.a.mul_mod_q(r[i]) + e[i];
    //     v[i] = v[i].scalar_mul_mod_q(pk.p);
    //     v[i].caddq();
    //
    //     w[i] = pk.t.mul_mod_q(r[i]) + e_prime[i];
    //     w[i] = w[i].scalar_mul_mod_q(pk.p).add_mod_q(m[i]);
    //     w[i].caddq();
    // }

    let mut rng = rand::thread_rng();
    loop {
        // Step 3: y ← DR^{4k},0,σ
        let y_r = DVector::from_fn(k, |_, _| R::random_gaussian(&mut rng, sigma));
        let y_e = DVector::from_fn(k, |_, _| R::random_gaussian(&mut rng, sigma));
        let y_e_prime = DVector::from_fn(k, |_, _| R::random_gaussian(&mut rng, sigma));
        let y_m = DVector::from_fn(k, |_, _| R::random_gaussian(&mut rng, sigma));
        let y = concat_vectors(&[y_r.clone(), y_e.clone(), y_e_prime.clone(), y_m.clone()]);

        // Compute commitment = b' y mod (q,q,p)
        let mut pa_yr = DVector::zeros(k);
        let mut p_ye = DVector::zeros(k);
        for i in 0..k {
            pa_yr[i] = pk.a.mul_mod_q(&y_r[i]).scalar_mul_mod_q(pk.p);
            p_ye[i] = y_e[i].scalar_mul_mod_q(pk.p);
        }
        let mut first = vector_add_mod_q(&pa_yr, &p_ye);

        let mut pt_yr = DVector::zeros(k);
        let mut p_ye_prime = DVector::zeros(k);
        for i in 0..k {
            pt_yr[i] = pk.t.mul_mod_q(&y_r[i]).scalar_mul_mod_q(pk.p);
            p_ye_prime[i] = y_e_prime[i].scalar_mul_mod_q(pk.p);
        }
        let mut second = vector_add_mod_q(&vector_add_mod_q(&pt_yr, &p_ye_prime), &y_m);
        let mut third = matrix_vector_mul_mod_p(b, &y_m);

        for i in 0..k {
            first[i].caddq();
            second[i].caddq();
        }
        for i in 0..ell {
            third[i].caddp();
        }
        let commitment = concat_vectors(&[first, second, third]);

        // Target = [v; w; u]
        let target = concat_vectors(&[v.clone(), w.clone(), u.clone()]);

        // Step 4: c ← H(...)
        let c = h(pk, b, &target, &commitment);

        // Step 5: s ← [r; e; e'; m] c
        let concat_secret = concat_vectors(&[r.clone(), e.clone(), e_prime.clone(), m.clone()]);
        // let s = vector_scalar_mul_mod(&concat_secret, c, pk.q);
        let mut s = DVector::zeros(concat_secret.len());
        for i in 0..concat_secret.len() {
            s[i] = concat_secret[i] * c;
        }

        // Step 6: z ← s + y
        let z = vector_add(&s, &y);

        // Step 7: Rejection sampling
        let dot_z_s = vector_dot_product(&z, &s);
        let norm_s_sq = vector_norm_l2_squared(&s);
        let exp_arg = -dot_z_s / sigma.powi(2) + norm_s_sq / (2.0 * sigma.powi(2));
        let prob = E.powf(exp_arg) / 3.0;
        if rng.gen_range(0.0..1.0) > prob {
            continue;
        }

        // Step 8: If ||z||_inf > 6 sigma
        if vector_norm_inf(&z) as f64 > 6.0 * sigma {
            continue;
        }
        // Step 9: Return
        return (c, z);
    }
}

pub fn verify(
    pk: &PublicKey,
    b: &DMatrix<R>,
    u: &DVector<R>,
    t: (&DVector<R>, &DVector<R>, &R, &DVector<R>),
) -> bool {
    let k = t.0.len();
    let ell = b.nrows();
    let sigma = 50688f64;

    let (v, w, c, z) = t;

    // Step 1: if ||z||_inf > 6 sigma, return false
    if (vector_norm_inf(&z) as f64) > 6.0 * sigma {
        return false;
    }

    // Split z into z_r, z_e, z_e_prime, z_m
    let z_r = DVector::from_vec(z.as_slice()[0..k].to_vec());
    let z_e = DVector::from_vec(z.as_slice()[k..2 * k].to_vec());
    let z_e_prime = DVector::from_vec(z.as_slice()[2 * k..3 * k].to_vec());
    let z_m = DVector::from_vec(z.as_slice()[3 * k..4 * k].to_vec());

    // Compute b' z - c * target mod (q, q, p)
    // first: p * (a * z_r + z_e) - c * v mod q
    let mut first = DVector::zeros(k);
    for i in 0..k {
        let pa_zr_i = pk.a.mul_mod_q(&z_r[i]).scalar_mul_mod_q(pk.p);
        let p_ze_i = z_e[i].scalar_mul_mod_q(pk.p);
        let sum = pa_zr_i + p_ze_i;
        let c_v_i = v[i].mul_mod_q(c);
        first[i] = sum.sub_mod_q(&c_v_i);
    }
    // second: p * (t * z_r + z_e_prime) + z_m - c * w mod q
    let mut second = DVector::zeros(k);
    for i in 0..k {
        let pt_zr_i = pk.t.mul_mod_q(&z_r[i]).scalar_mul_mod_q(pk.p);
        let p_ze_prime_i = z_e_prime[i].scalar_mul_mod_q(pk.p);
        let temp = pt_zr_i + p_ze_prime_i + z_m[i];
        let c_w_i = w[i].mul_mod_q(c);
        second[i] = temp.sub_mod_q(&c_w_i);
    }
    // third: b * z_m - c * u mod p
    let mut third = DVector::zeros(ell);
    for i in 0..ell {
        let mut sum = R::zero();
        for j in 0..k {
            sum = sum.add_mod_p(&b[(i, j)].mul_mod_p(&z_m[j]));
            // sum += b[(i, j)] * z_m[j];
        }
        let c_u_i = u[i].mul_mod_p(c);
        third[i] = sum.sub_mod_p(&c_u_i);
    }
    for i in 0..k {
        first[i].caddq();
        second[i].caddq();
    }
    for i in 0..ell {
        third[i].caddp();
    }
    // commitment = [first; second; third]
    let commitment = concat_vectors(&[first, second, third]);

    // target = [v; w; u]
    let target = concat_vectors(&[v.clone(), w.clone(), u.clone()]);

    let c_prime = h(pk, b, &target, &commitment);
    // Step 2: if c != H(pk, b, target, commitment), return false
    if c != &c_prime {
        return false;

    }

    // Step 3: return true
    true
}

pub fn ring_lwe_dec(
    sk: &SecretKey,
    v: &DVector<R>,
    w: &DVector<R>,
) -> DVector<R> {
    let k = v.len();
    let mut m = DVector::zeros(k);
    for i in 0..k {
        let v_s1 = v[i].mul_mod_q(&sk.s1);
        m[i] = w[i].sub_mod_q(&v_s1);
        m[i].mod_p();
    }
    m
}

pub fn decrypt(
    pk: &PublicKey,
    sk: &SecretKey,
    // B: &DMatrix<R>,
    // u: &DVector<R>,
    t: (&DVector<R>, &DVector<R>, &R, &DVector<R>),
) -> Option<DVector<R>> {
    let k = t.0.len();
    let sigma = 50688f64;
    let c_max = MAX_C_L1 * 2; // max ||c - c'||_1 for c, c' in C

    // Step 1: Verify ciphertext
    // if !verify(pk, B, u, t.clone()) {
    //     return None;
    // }

    let (v, w, c, _z) = t;

    // Step 2: guess the honest encryptor case.
    let c_bar = R::one();
    // Step 3: Compute m = (w - v s1) c_bar mod q
    let mut m = DVector::zeros(k);
    for i in 0..k {
        let v_s1 = v[i].mul_mod_q(&sk.s1);
        let w_minus_vs1 = w[i] - v_s1;
        m[i] = w_minus_vs1.mul_mod_q(&c_bar);
    }

    let mut valid = true;
    // Step 4: Check norms
    for i in 0..k {
        if m[i].norm_inf() as f64 > (pk.q as f64) / (2.0 * c_max as f64) {
            valid = false;
        }
        m[i].mod_p();
        if m[i].norm_inf() as f64 >= 12.0 * sigma {
            valid = false;
        }
    }
    if valid {
        // Step 5: Return (m mod p, c_bar)
        // let m_bar = DVector::from_fn(k, |i, _| {
        //     let mut temp = m[i];
        //     temp.mod_p();
        //     temp
        // });
        return Some(m);
    } else {
        panic!("decryption error");
    }
    // otherwise, malicious encryptor provides a mˉ satisfying m/c=mˉ/cˉmod p, we have to guess and recover m
    //TODO: the malicious encryptor case has not been tested
    // let mut rng = thread_rng();
    // loop {
    //     // Guess c' ∈ C
    //     let c_prime = generate_low_l1_r(&mut rng);
    //
    //     // Compute c_bar = c - c'
    //     let c_bar = *c - c_prime;
    //
    //     // Compute m = (w - v s1) c_bar mod q
    //     let mut m = DVector::zeros(k);
    //     valid = true;
    //     for i in 0..k {
    //         let v_s1 = v[i] * sk.s1;
    //         let w_minus_vs1 = w[i] - v_s1;
    //         m[i] = w_minus_vs1.mul_mod_q(&c_bar);
    //         let mut temp = m[i];
    //         temp.mod_p();
    //         if m[i].norm_inf() as f64 > (pk.q as f64) / (2.0 * c_max as f64)
    //             || temp.norm_inf() as f64 >= 12.0 * sigma
    //         {
    //             valid = false;
    //             break;
    //         }
    //     }
    //     if valid {
    //         let m_bar = DVector::from_fn(k, |i, _| {
    //             let mut temp = m[i];
    //             temp.mod_p();
    //             temp
    //         });
    //         // Compute original m = m_bar * (c * inv(c_bar)) mod p : p must be 5 mod 8 according to Lemma 2.2 (making c_bar invertible)
    //         let inv_c_bar = c_bar.inverse_mod(pk.p).unwrap();
    //         let c_inv_c_bar = c.mul_mod_p(&inv_c_bar);
    //         let original_m = DVector::from_fn(k, |i, _| m_bar[i].mul_mod_p(&c_inv_c_bar));
    //         return Some(original_m);
    //     }
    // }
}
