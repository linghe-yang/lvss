use crate::r_ring::{N, R};

pub const P_PRIME: i32 = 32789;

pub const LOG_P: usize = 16;
// pub const P_PRIME: i32 = 13;

pub fn reduce32_p(a: i32) -> i32 {
    let mut t = (a + (1 << 14)) >> 15;
    t = a - t.wrapping_mul(P_PRIME);
    t
}

/// Reduces an i64 input modulo p to [-p/2, p/2].
pub const P_HALF: i32 = P_PRIME >> 1;
pub fn reduce64_p(a: i64) -> i32 {
    let q = P_PRIME as i64;
    let r = a % q;
    let r_i32 = r as i32;
    if r_i32 > P_HALF {
        r_i32 - P_PRIME
    } else if r_i32 < -P_HALF {
        r_i32 + P_PRIME
    } else {
        r_i32
    }
}

fn caddp(a: i32) -> i32 {
    a + ((a >> 31) & P_PRIME)
}

impl R {
    /// Reduces each coefficient of the polynomial modulo p to [-p/2, p/2].
    pub fn mod_p(&mut self) {
        for i in 0..N {
            self.coeffs[i] = reduce32_p(self.coeffs[i]);
        }
    }

    pub fn caddp(&mut self) {
        for i in 0..N {
            self.coeffs[i] = caddp(self.coeffs[i]);
        }
    }

    /// Multiplies the polynomial by a scalar, reducing each coefficient to [-p/2, p/2].
    pub fn scalar_mul_mod_p(&self, scalar: i32) -> Self {
        let mut result = R { coeffs: [0; N] };
        let scalar_64 = scalar as i64;
        for i in 0..N {
            let product = (self.coeffs[i] as i64) * scalar_64;
            result.coeffs[i] = reduce64_p(product);
        }
        result
    }

    /// Multiplies two polynomials in Z[x]/(x^N + 1) modulo p, reducing coefficients to [-p/2, p/2].
    pub fn mul_mod_p(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        let mut temp = [0i64; N];

        for i in 0..N {
            for j in 0..N {
                let k = (i + j) % N;
                let sign = if i + j >= N { -1 } else { 1 };
                let product = (self.coeffs[i] as i64) * (other.coeffs[j] as i64) * (sign as i64);
                temp[k] = temp[k].wrapping_add(product);
            }
        }

        for i in 0..N {
            result.coeffs[i] = reduce64_p(temp[i]);
        }

        result
    }

    /// Adds two polynomials modulo p, reducing coefficients to [-p/2, p/2].
    pub fn add_mod_p(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        for i in 0..N {
            // Add and reduce to [-p/2, p/2]
            result.coeffs[i] = reduce32_p(self.coeffs[i].wrapping_add(other.coeffs[i]));
        }
        result
    }

    /// Subtracts other polynomial from self modulo p, reducing coefficients to [-p/2, p/2].
    pub fn sub_mod_p(&self, other: &Self) -> Self {
        let mut result = R { coeffs: [0; N] };
        for i in 0..N {
            // Subtract and reduce to [-p/2, p/2]
            result.coeffs[i] = reduce32_p(self.coeffs[i].wrapping_sub(other.coeffs[i]));
        }
        result
    }

    pub fn equal_mod_p(&self, other: &Self) -> bool {
        let mut a = self.clone();
        let mut b = other.clone();
        a.caddp();
        b.caddp();
        a == b
    }

}