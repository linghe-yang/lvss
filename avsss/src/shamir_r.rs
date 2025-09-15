use std::collections::HashSet;
use ve::r_ring::{N, R};
use ve::rp::P_PRIME;

// Helper function to compute (a % m) in positive range
fn mod_positive(a: i64, m: i64) -> i64 {
    ((a % m) + m) % m
}

// Extended Euclidean algorithm
fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

// Modular inverse using extended gcd
fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(a, m);
    if g != 1 {
        None
    } else {
        Some(mod_positive(x, m))
    }
}

// Lagrange interpolation at a target point for scalar values in mod p
fn lagrange_interpolate_at(target: i32, points: &[(i32, i32)], p: i32) -> Option<i32> {
    let n = points.len();
    let pp = p as i64;
    let target_mod = mod_positive(target as i64, pp);

    let mut xs_mod: Vec<i64> = Vec::with_capacity(n);
    let mut ys_mod: Vec<i64> = Vec::with_capacity(n);

    for &(x, y) in points {
        xs_mod.push(mod_positive(x as i64, pp));
        ys_mod.push(mod_positive(y as i64, pp));
    }

    // Check for distinct x_mod
    let mut seen = HashSet::new();
    for &x in &xs_mod {
        if !seen.insert(x) {
            return None; // Duplicates mod p
        }
    }

    let mut result: i64 = 0;

    for i in 0..n {
        let x_i = xs_mod[i];
        let y_i = ys_mod[i];

        let mut num: i64 = 1;
        let mut den: i64 = 1;

        for j in 0..n {
            if i == j {
                continue;
            }
            let x_j = xs_mod[j];

            let diff_num = mod_positive(target_mod - x_j, pp);
            num = mod_positive(num * diff_num, pp);

            let diff_den = mod_positive(x_i - x_j, pp);
            den = mod_positive(den * diff_den, pp);
        }

        let inv_den = match mod_inverse(den, pp) {
            Some(inv) => inv,
            None => return None,
        };

        let l_i = mod_positive(num * inv_den, pp);
        let contrib = mod_positive(y_i * l_i, pp);
        result = mod_positive(result + contrib, pp);
    }

    Some(result as i32)
}

// Interpolate R at target by interpolating each coefficient independently
fn interpolate_r_at(target: i32, points: &[(i32, R)], p: i32) -> Option<R> {
    let n = points.len();
    let mut coeffs = [0i32; N];

    for j in 0..N {
        let scalar_points: Vec<(i32, i32)> = points.iter().map(|&(x, ref r)| (x, r.coeffs[j])).collect();
        coeffs[j] = match lagrange_interpolate_at(target, &scalar_points, p) {
            Some(val) => val,
            None => return None,
        };
    }

    Some(R { coeffs })
}

// Modular reduction for i32
fn mod_i32(a: i32, p: i32) -> i32 {
    let pp = p as i64;
    let aa = a as i64;
    mod_positive(aa, pp) as i32
}

// Compare two R modulo p
fn eq_mod_r(a: &R, b: &R, p: i32) -> bool {
    for j in 0..N {
        if mod_i32(a.coeffs[j], p) != mod_i32(b.coeffs[j], p) {
            return false;
        }
    }
    true
}

pub fn shamir_reconstruct_r(shares: &[(i32, R)], t: usize) -> Option<R> {
    let p = P_PRIME;
    let n_shares = shares.len();
    let threshold = t + 1;

    if n_shares < threshold {
        return None;
    }

    // Check distinct x
    let mut xs: Vec<i32> = shares.iter().map(|&(x, _)| x).collect();
    xs.sort();
    for i in 1..n_shares {
        if xs[i] == xs[i - 1] {
            return None;
        }
    }

    // Sort shares by x for determinism
    let mut sorted_shares: Vec<(i32, R)> = shares.to_vec();
    sorted_shares.sort_by_key(|&(x, _)| x);

    // Reconstruct using first threshold points
    let recon_points = &sorted_shares[0..threshold];
    let secret = match interpolate_r_at(0, recon_points, p) {
        Some(s) => s,
        None => return None,
    };

    // If exactly threshold, return secret
    if n_shares == threshold {
        return Some(secret);
    }

    // Verify remaining shares
    for extra in &sorted_shares[threshold..] {
        let (x, ref y) = *extra;
        let computed_y = match interpolate_r_at(x, recon_points, p) {
            Some(cy) => cy,
            None => return None,
        };
        if !eq_mod_r(&computed_y, y, p) {
            return None;
        }
    }

    Some(secret)
}