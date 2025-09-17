use std::collections::HashSet;
use ve::r_ring::{N, R};
use ve::rp::P_PRIME;
use crate::components::ID;

// Helper function to compute (a % m) in positive range
fn mod_positive(a: i128, m: i128) -> i128 {
    ((a % m) + m) % m
}

// Extended Euclidean algorithm
fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

// Modular inverse using extended gcd
fn mod_inverse(a: i128, m: i128) -> Option<i128> {
    let (g, x, _) = extended_gcd(a, m);
    if g != 1 {
        None
    } else {
        Some(mod_positive(x, m))
    }
}

// Solve linear system A x = b mod p using Gaussian elimination
fn solve_system_mod_p(mut a: Vec<Vec<i128>>, mut b: Vec<i128>, p: i128) -> Option<Vec<i128>> {
    let rows = a.len();
    let cols = if rows > 0 { a[0].len() } else { 0 };

    // Augment A with b
    for i in 0..rows {
        a[i].push(mod_positive(b[i], p));
    }

    // Gaussian elimination with partial pivoting
    let mut h = 0; // pivot row
    let mut k = 0; // pivot column
    while h < rows && k < cols {
        // Find pivot
        let mut pivot = h;
        for i in h + 1..rows {
            if a[i][k].abs() > a[pivot][k].abs() {
                pivot = i;
            }
        }
        if a[pivot][k] == 0 {
            k += 1;
            continue;
        }
        // Swap rows
        a.swap(h, pivot);
        // Get inverse of pivot
        let piv = a[h][k];
        let inv = match mod_inverse(piv, p) {
            Some(i) => i,
            None => return None,
        };
        // Normalize pivot row
        for j in k..cols + 1 {
            a[h][j] = mod_positive(a[h][j] * inv, p);
        }
        // Eliminate other rows
        for i in 0..rows {
            if i == h {
                continue;
            }
            let factor = a[i][k];
            for j in k..cols + 1 {
                a[i][j] = mod_positive(a[i][j] - factor * a[h][j], p);
            }
        }
        h += 1;
        k += 1;
    }

    // Check consistency and extract solution
    let mut x = vec![0i128; cols];
    let mut rank = 0;
    for i in 0..cols {
        if i >= rows {
            break;
        }
        if a[i][i] == 1 {
            x[i] = a[i][cols];
            rank += 1;
        } else {
            // Check if inconsistent
            let mut all_zero = true;
            for j in i..cols {
                if a[i][j] != 0 {
                    all_zero = false;
                    break;
                }
            }
            if all_zero && a[i][cols] != 0 {
                return None;
            }
        }
    }

    if rank < cols {
        return None; // Not unique solution
    }

    for val in &mut x {
        *val = mod_positive(*val, p);
    }

    Some(x)
}

// Polynomial division mod p: num / den = quot + rem / den, deg(rem) < deg(den)
fn poly_div_mod_p(num: &[i128], den: &[i128], p: i128) -> Option<(Vec<i128>, Vec<i128>)> {
    if den.is_empty() || *den.last().unwrap() == 0 {
        return None;
    }
    let mut rem: Vec<i128> = num.iter().map(|&v| mod_positive(v, p)).collect();
    let d_deg = den.len() - 1;
    let ld_den = mod_positive(*den.last().unwrap(), p);
    let inv_ld = match mod_inverse(ld_den, p) {
        Some(inv) => inv,
        None => return None,
    };
    let q_deg = if rem.len() <= d_deg { 0 } else { rem.len() - d_deg - 1 };
    let mut quot = vec![0i128; q_deg + 1];

    for ii in (0..=q_deg).rev() {
        // Trim leading zeros from rem
        while !rem.is_empty() && *rem.last().unwrap() == 0 {
            rem.pop();
        }
        if rem.len() <= d_deg + ii {
            continue;
        }
        let r_deg = rem.len() - 1;
        let ld_rem = *rem.last().unwrap();
        let coeff = mod_positive(ld_rem * inv_ld, p);
        quot[ii] = coeff;
        for jj in 0..den.len() {
            let rem_idx = r_deg - (d_deg - jj);
            let sub = mod_positive(coeff * mod_positive(den[jj], p), p);
            rem[rem_idx] = mod_positive(rem[rem_idx] - sub, p);
        }
    }
    // Trim leading zeros from rem
    while !rem.is_empty() && *rem.last().unwrap() == 0 {
        rem.pop();
    }
    // Trim leading zeros from quot
    while quot.len() > 1 && *quot.last().unwrap() == 0 {
        quot.pop();
    }
    Some((quot, rem))
}

pub fn shamir_reconstruct_rs(shares: &[(ID, R)], t: usize) -> Option<R> {
    let p = P_PRIME;
    let n = shares.len();
    let k = t + 1;
    if n < 2 * t + 1 {
        return None;
    }

    // Sort shares by x
    let mut sorted_shares: Vec<(ID, R)> = shares.to_vec();
    sorted_shares.sort_by_key(|&(x, _)| x);

    // Check distinct x
    let mut xs: Vec<ID> = sorted_shares.iter().map(|&(x, _)| x).collect();
    let mut seen = HashSet::new();
    for &x in &xs {
        if !seen.insert(x) {
            return None;
        }
    }

    let pp = p as i128;
    let xs_mod: Vec<i128> = xs.iter().map(|&x| mod_positive(x as i128, pp)).collect();

    // Check distinct mod p
    let mut seen_mod = HashSet::new();
    for &x in &xs_mod {
        if !seen_mod.insert(x) {
            return None;
        }
    }

    let d = t;
    let e = (n - k) / 2;
    let v = d + 2 * e + 1; // number of variables

    let mut secret_coeffs = [0i64; N];

    for j in 0..N {
        let ys: Vec<i128> = sorted_shares.iter().map(|&(_, r)| mod_positive(r.coeffs[j] as i128, pp)).collect();

        // Build the system A x = b
        let mut a: Vec<Vec<i128>> = vec![vec![0i128; v]; n];
        let mut bb: Vec<i128> = vec![0i128; n];

        for i in 0..n {
            let x = xs_mod[i];
            let y = ys[i];

            // Coefficients for q_0 to q_{d+e}
            let mut xp = 1i128;
            for jj in 0..=(d + e) {
                a[i][jj] = xp;
                xp = mod_positive(xp * x, pp);
            }

            // Coefficients for e_0 to e_{e-1}
            let mut xp_e = 1i128;
            for jj in 0..e {
                let col = d + e + 1 + jj;
                a[i][col] = mod_positive(-y * xp_e, pp);
                xp_e = mod_positive(xp_e * x, pp);
            }

            // b_i = y * x^e
            bb[i] = mod_positive(y * xp_e, pp);
        }

        // Solve the system
        let sol = match solve_system_mod_p(a, bb, pp) {
            Some(s) => s,
            None => return None,
        };

        // Extract Q and E
        let q_vec: Vec<i128> = sol[0..=d + e].to_vec();
        let mut e_vec: Vec<i128> = if e > 0 { sol[d + e + 1..].to_vec() } else { vec![] };
        e_vec.push(1); // leading coefficient 1 for E

        // Polynomial division Q / E
        let (quot, rem) = match poly_div_mod_p(&q_vec, &e_vec, pp) {
            Some((q, r)) => (q, r),
            None => return None,
        };

        // Check if division is exact and deg(quot) <= d
        if !rem.iter().all(|&r| r == 0) || quot.len() > d + 1 {
            return None;
        }

        // Secret coefficient is quot[0] mod p
        secret_coeffs[j] = mod_positive(quot[0], pp) as i64;
    }

    Some(R { coeffs: secret_coeffs })
}