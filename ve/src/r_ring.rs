use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Neg};
use std::default::Default;
use std::fmt::{self, Display, Formatter};
use rand::Rng;
use rand_distr::{Normal, Distribution};
use rand_distr::num_traits::Zero;
use crate::rp::P_PRIME;
use crate::rq::Q_PRIME;
use crate::util::poly_extended_euclid;
use serde::{Serialize, Serializer, Deserialize, Deserializer};
use serde::de::{Error, Visitor};

pub const N: usize = 256;
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct R {
    pub coeffs: [i64; N],
}

impl Zero for R {
    fn zero() -> Self {
        Self { coeffs: [0; N] }
    }

    fn is_zero(&self) -> bool {
        self.coeffs == [0; N]
    }
}

impl R {

    /// Returns the multiplicative identity (one) of the ring.
    pub fn one() -> Self {
        let mut coeffs = [0; N];
        coeffs[0] = 1;
        Self { coeffs }
    }

    /// Generates a random polynomial with coefficients uniformly sampled from [a, b].
    pub fn random_uniform(a: i64, b: i64) -> Self {
        let mut rng = rand::thread_rng();
        let mut coeffs = [0; N];
        for coeff in coeffs.iter_mut() {
            *coeff = rng.gen_range(a..=b);
        }
        Self { coeffs }
    }

    pub fn random_in_p() -> Self {
        Self::random_uniform(-P_PRIME/2, P_PRIME/2)
    }

    pub fn random_in_q() -> Self {
        Self::random_uniform(-Q_PRIME/2, Q_PRIME/2)
    }

    /// Generates a random polynomial with coefficients sampled from a discrete Gaussian distribution
    /// with mean 0 and standard deviation sigma:D_{R,0,σ} Samples are rounded to nearest integer.
    pub fn random_gaussian<RR: Rng>(mut rng: &mut RR, sigma: f64) -> Self {
        let normal = Normal::new(0.0, sigma).unwrap();
        let mut coeffs = [0; N];
        for coeff in coeffs.iter_mut() {
            *coeff = normal.sample(&mut rng).round() as i64;
        }
        Self { coeffs }
    }

    /// Computes the 1-norm (sum of absolute values of coefficients).
    pub fn norm_l1(&self) -> i64 {
        self.coeffs.iter().map(|&x| x.abs()).sum()
    }

    /// Computes the 2-norm (square root of sum of squares of coefficients).
    pub fn norm_l2(&self) -> f64 {
        (self.coeffs.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()).sqrt()
    }

    /// Computes the infinity norm (maximum absolute value of coefficients).
    pub fn norm_inf(&self) -> i64 {
        self.coeffs.iter().map(|&x| x.abs()).max().unwrap_or(0)
    }
    pub fn inverse_mod(&self, p: i64) -> Option<R> {
        let g_coeffs = self.coeffs.iter().map(|&c| c % p).collect::<Vec<i64>>();
        let mut f_coeffs = vec![1];
        f_coeffs.extend(vec![0; N - 1]);
        f_coeffs.push(1);

        let (d, s, _) = poly_extended_euclid(g_coeffs, f_coeffs, p);
        if d != vec![1] {
            None
        } else {
            let mut inv_coeffs = [0; N];
            for (i, &c) in s.iter().enumerate() {
                if i >= N {
                    return None; // Degree too high
                }
                inv_coeffs[i] = c;
            }
            Some(R { coeffs: inv_coeffs })
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(N * 8);
        for &coeff in self.coeffs.iter() {
            bytes.extend_from_slice(&coeff.to_le_bytes());
        }
        bytes
    }
}

impl Default for R {
    fn default() -> Self {
        Self::zero()
    }
}

impl From<i64> for R {
    /// Converts a constant integer into a constant polynomial in the ring.
    fn from(val: i64) -> Self {
        let mut coeffs = [0; N];
        coeffs[0] = val;
        Self { coeffs }
    }
}

impl Add for R {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut res = [0; N];
        for i in 0..N {
            res[i] = self.coeffs[i] + other.coeffs[i];
        }
        Self { coeffs: res }
    }
}

impl AddAssign for R {
    fn add_assign(&mut self, other: Self) {
        for i in 0..N {
            self.coeffs[i] += other.coeffs[i];
        }
    }
}

impl Sub for R {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut res = [0; N];
        for i in 0..N {
            res[i] = self.coeffs[i] - other.coeffs[i];
        }
        Self { coeffs: res }
    }
}

impl SubAssign for R {
    fn sub_assign(&mut self, other: Self) {
        for i in 0..N {
            self.coeffs[i] -= other.coeffs[i];
        }
    }
}

impl Neg for R {
    type Output = Self;

    fn neg(self) -> Self {
        let mut res = [0; N];
        for i in 0..N {
            res[i] = -self.coeffs[i];
        }
        Self { coeffs: res }
    }
}

impl Mul for R {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut res = [0; N];
        for i in 0..N {
            for j in 0..N {
                let idx = (i + j) % N;
                let sign = if (i + j) / N == 0 { 1 } else { -1 };
                res[idx] += sign * self.coeffs[i] * other.coeffs[j];
            }
        }
        Self { coeffs: res }
    }
}

impl MulAssign for R {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Mul<i64> for R {
    type Output = Self;

    fn mul(self, scalar: i64) -> Self {
        let mut coeffs = [0; N];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            coeffs[i] = coeff.wrapping_mul(scalar);
        }
        Self { coeffs }
    }
}

impl Mul<R> for i64 {
    type Output = R;

    fn mul(self, rq: R) -> R {
        rq * self
    }
}

impl Display for R {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.coeffs)
    }
}


impl Serialize for R {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(
                self.coeffs.as_ptr() as *const u8,
                N * size_of::<i64>(),
            )
        };
        let encoded = base64::encode(bytes);
        serializer.serialize_str(&encoded)
    }
}

struct RVisitor;

impl<'de> Visitor<'de> for RVisitor {
    type Value = R;

    fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
        formatter.write_str("a base64 string representing 256 i64 values")
    }

    fn visit_str<E>(self, v: &str) -> Result<R, E>
    where
        E: Error,
    {
        let decoded = base64::decode(v).map_err(E::custom)?;
        if decoded.len() != N * size_of::<i64>() {
            return Err(E::custom("invalid length"));
        }
        let mut coeffs = [0i64; N];
        unsafe {
            std::ptr::copy_nonoverlapping(
                decoded.as_ptr(),
                coeffs.as_mut_ptr() as *mut u8,
                decoded.len(),
            );
        }

        Ok(R { coeffs })
    }
}

impl<'de> Deserialize<'de> for R {
    fn deserialize<D>(deserializer: D) -> Result<R, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(RVisitor)
    }
}

#[test]
fn test() {
    let r = R::random_uniform(-1000,1000);

    let encoded = serde_json::to_string(&r).unwrap();
    println!("Encoded: {}", encoded);

    let decoded: R = serde_json::from_str(&encoded).unwrap();
    println!("Equal: {}", decoded == r);
}