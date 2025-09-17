use crate::barrett::symmetric_barrett_reduce;
use crate::r_ring::{N, R};
use crate::zetas::ZETAS_256;

// pub const Q_PRIME: i64 = (1 << 23) - (1 << 13) + 1;
// pub const Q_INV: i64 = 1732267787797143553; // q^(-1) mod 2^32
//
// const MULT_Q: i128 = 2201172575746; // Precomputed for q_prime = 3
// const HALF_Q: i64 = 4190208; // Precomputed for q_prime = 3

pub const Q_PRIME: i64 = 34359724033;
pub const Q_INV: i64 = -5903494919661537279; // q^(-1) mod 2^64

const MULT_Q: i128 = 536871136;
const Q_HALF: i64 = 17179862016;


pub fn reduce32_q(a: i64) -> i64 {
    // let mut t = (a + (1 << 22)) >> 23;
    // t = a - t.wrapping_mul(Q_PRIME);
    // t
    symmetric_barrett_reduce(a as i128, Q_PRIME, MULT_Q, Q_HALF)
}
pub fn reduce64_q(a: i128) -> i64 {
    // let q = Q_PRIME as i128;
    // let r = a % q;
    // let r_i64 = r as i64;
    // if r_i64 > Q_HALF {
    //     r_i64 - Q_PRIME
    // } else if r_i64 < -Q_HALF {
    //     r_i64 + Q_PRIME
    // } else {
    //     r_i64
    // }
    symmetric_barrett_reduce(a, Q_PRIME, MULT_Q, Q_HALF)
}

fn caddq(a: i64) -> i64 {
    // In C right-shift of negative signed integers is implementation-defined, so C reference implementation contains bug.
    // In Rust if a < 0 right-shift is defined to fill with 1s, 0s otherwise, so we're bug free here.
    a + ((a >> 63) & Q_PRIME)
}

/// For integer a with -2^{31} * Q <= a <= 2^31 * Q,
/// compute r \equiv 2^{-32} * a (mod Q) such that -Q < r < Q.
///
/// Returns r.
pub fn montgomery_reduce(a: i128) -> i64 {
    let mut t = (a as i64).wrapping_mul(Q_INV) as i128;
    t = (a as i128 - t.wrapping_mul(Q_PRIME as i128)) >> 64;
    t as i64
}

impl R {
    /// Reduces each coefficient of the polynomial modulo p to [-p/2, p/2].
    pub fn mod_q(&mut self) {
        for i in 0..N {
            self.coeffs[i] = reduce32_q(self.coeffs[i]);
        }
    }

    pub fn caddq(&mut self) {
        for i in 0..N {
            self.coeffs[i] = caddq(self.coeffs[i]);
        }
    }

    /// Multiplies the polynomial by a scalar, reducing each coefficient to [-p/2, p/2].
    pub fn scalar_mul_mod_q(&self, scalar: i64) -> Self {
        let mut result = R { coeffs: [0; N] };
        let scalar_64 = scalar as i128;
        for i in 0..N {
            // Use i128 to avoid overflow in multiplication
            let product = (self.coeffs[i] as i128) * scalar_64;
            // Directly reduce to [-p/2, p/2] using reduce32_p
            result.coeffs[i] = reduce64_q(product);
        }
        result
    }

    /// Multiplies two polynomials in Z[x]/(x^N + 1) modulo p, reducing coefficients to [-p/2, p/2]. N must be 256
    pub fn mul_mod_q(&self, other: &Self) -> Self {
        poly_mul_q_256(self, &other)
    }

    pub fn mul_mod_q_conv(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        let mut temp = [0i128; N];

        for i in 0..N {
            for j in 0..N {
                let k = (i + j) % N;
                let sign = if i + j >= N { -1 } else { 1 };
                let product = (self.coeffs[i] as i128) * (other.coeffs[j] as i128) * (sign as i128);
                temp[k] = temp[k].wrapping_add(product);
            }
        }

        for i in 0..N {
            result.coeffs[i] = reduce64_q(temp[i]);
        }

        result
    }

    /// Adds two polynomials modulo p, reducing coefficients to [-p/2, p/2].
    pub fn add_mod_q(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        for i in 0..N {
            result.coeffs[i] = reduce32_q(self.coeffs[i].wrapping_add(other.coeffs[i]));
        }
        result
    }

    /// Subtracts other polynomial from self modulo p, reducing coefficients to [-p/2, p/2].
    pub fn sub_mod_q(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        for i in 0..N {
            // Subtract and reduce to [-p/2, p/2]
            result.coeffs[i] = reduce32_q(self.coeffs[i].wrapping_sub(other.coeffs[i]));
        }
        result
    }

    pub fn equal_mod_q(&self, other: &Self) -> bool {
        let mut a = self.clone();
        let mut b = other.clone();
        a.caddq();
        b.caddq();
        a == b
    }
}

fn pointwise_montgomery(c: &mut R, a: &R, b: &R) {
    for i in 0..N {
        c.coeffs[i] = montgomery_reduce(a.coeffs[i] as i128 * b.coeffs[i] as i128);
    }
}

fn ntt_256(a: &mut [i64]) {
    let mut k: usize = 0;
    let mut len: usize = 128;

    while len > 0 {
        let mut start: usize = 0;
        while start < N {
            k += 1;
            let zeta = ZETAS_256[k] as i128;
            let mut j = start;
            while j < (start + len) {
                let t = montgomery_reduce(zeta.wrapping_mul(a[j + len] as i128));
                a[j + len] = a[j] - t;
                a[j] += t;
                j += 1;
            }
            start = j + len;
        }
        len >>= 1;
    }
}

fn invntt_tomont_256(a: &mut [i64]) {
    let mut k: usize = 256;
    let mut len: usize = 1;
    const F: i128 = 7320386123; // mont^2/256
    // const F: i128 = 41978; // mont^2/256

    while len < N {
        let mut start: usize = 0;
        while start < 256 {
            k -= 1;
            let zeta = -ZETAS_256[k] as i128;
            let mut j = start;
            while j < (start + len) {
                let t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = montgomery_reduce(zeta.wrapping_mul(a[j + len] as i128));
                j += 1
            }
            start = j + len;
        }
        len <<= 1;
    }
    for j in 0..N {
        a[j] = montgomery_reduce(F.wrapping_mul(a[j] as i128));
    }
}

fn reduce(a: &mut R) {
    // Bad C style
    // for i in 0..N {
    //     a.coeffs[i] = reduce::reduce32(a.coeffs[i]);
    // }
    // Nice Rust style
    for coeff in a.coeffs.iter_mut() {
        *coeff = reduce32_q(*coeff);
    }
}

fn poly_mul_q_256(a: &R, b: &R) -> R {
    let mut a_ntt = a.clone();
    let mut b_ntt = b.clone();
    ntt_256(&mut a_ntt.coeffs);
    ntt_256(&mut b_ntt.coeffs);
    let mut c_ntt = R::default();
    pointwise_montgomery(&mut c_ntt, &a_ntt, &b_ntt);
    invntt_tomont_256(&mut c_ntt.coeffs);
    reduce(&mut c_ntt);
    c_ntt
}
