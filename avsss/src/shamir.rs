use nalgebra::DVector;
use ve::r_ring::R;
use ve::rp::P_PRIME;
use ve::util::random_dvector;

pub fn shamir_share(s: &DVector<R>, n: usize, t: usize) -> Vec<(i64, DVector<R>)> {
    let dim = s.len();
    let p = P_PRIME;
    let mut coeffs: Vec<DVector<R>> = vec![s.clone()];
    for _ in 1..=t {
        coeffs.push(random_dvector(dim));
    }

    let mut shares = Vec::with_capacity(n);
    for i in 1..=n {
        let x = i as i64;
        let mut eval = coeffs[0].clone();
        let mut x_pow = 1i64;
        for k in 1..=t {
            x_pow = mod_mul(x_pow, x, p);
            let scalar = x_pow;
            let term = coeffs[k].map(|r| r.scalar_mul_mod_p(scalar));
            eval = eval.zip_map(&term, |a, b| a.add_mod_p(&b));
        }
        shares.push((x, eval));
    }
    shares
}
pub fn shamir_reconstruct(shares: &[(i64, DVector<R>)], t: usize) -> Option<DVector<R>> {
    let p = P_PRIME;
    let m = shares.len();
    if m < t + 1 {
        return None;
    }

    // Use first t+1 shares to reconstruct the secret
    let base_shares: Vec<_> = shares.iter().take(t + 1).cloned().collect();
    let reconstructed = reconstruct_base(&base_shares, p);
    // If exactly t+1, no need to verify
    if m == t + 1 {
        return Some(reconstructed);
    }

    // For additional shares, verify they lie on the same polynomial
    for share in shares.iter().skip(t + 1) {
        let x = share.0;
        let y = &share.1;
        let expected = lagrange_eval(&base_shares, x, p);

        // Check if expected == y element-wise
        if !expected.iter().zip(y.iter()).all(|(e, yy)| e.equal_mod_p(yy)) {
            return None;
        }
    }

    Some(reconstructed)
}


fn reconstruct_base(shares: &[(i64, DVector<R>)], p: i64) -> DVector<R> {
    let m = shares.len(); // Exactly t+1
    let dim = shares[0].1.len();
    let zero_r = R::default();
    let mut secret = DVector::from_fn(dim, |_, _| zero_r);

    for j in 0..m {
        let xj = shares[j].0;
        let yj = &shares[j].1;

        let mut l_num = 1i64;
        let mut l_den = 1i64;
        for k in 0..m {
            if k == j {
                continue;
            }
            let xk = shares[k].0;
            let num_factor = mod_diff(0, xk, p); // For f(0)
            l_num = mod_mul(l_num, num_factor, p);
            let den_factor = mod_diff(xj, xk, p);
            l_den = mod_mul(l_den, den_factor, p);
        }

        let inv_den = mod_inverse(l_den, p);
        let l0 = mod_mul(l_num, inv_den, p);

        let scaled = yj.map(|r| r.scalar_mul_mod_p(l0));
        secret = secret.zip_map(&scaled, |a, b| a.add_mod_p(&b));
    }

    secret
}


// Helper function to evaluate the polynomial at a new point x_new using exactly t+1 shares
fn lagrange_eval(shares: &[(i64, DVector<R>)], x_new: i64, p: i64) -> DVector<R> {
    let m = shares.len(); // Exactly t+1
    let dim = shares[0].1.len();
    let zero_r = R::default();
    let mut result = DVector::from_fn(dim, |_, _| zero_r);

    for i in 0..m {
        let xi = shares[i].0;
        let yi = &shares[i].1;

        let mut num = 1i64;
        let mut den = 1i64;
        for j in 0..m {
            if i == j {
                continue;
            }
            let xj = shares[j].0;
            num = mod_mul(num, mod_diff(x_new, xj, p), p);
            den = mod_mul(den, mod_diff(xi, xj, p), p);
        }

        let inv_den = mod_inverse(den, p);
        let li = mod_mul(num, inv_den, p);

        let term = yi.map(|r| r.scalar_mul_mod_p(li));
        result = result.zip_map(&term, |a, b| a.add_mod_p(&b));
    }

    result
}

fn mod_mul(a: i64, b: i64, p: i64) -> i64 {
    let product = (a as i128 * b as i128) % (p as i128);
    ((product + p as i128) % p as i128) as i64
}

fn mod_diff(a: i64, b: i64, p: i64) -> i64 {
    let diff = (a as i128 - b as i128) % (p as i128);
    ((diff + p as i128) % p as i128) as i64
}

fn mod_inverse(a: i64, p: i64) -> i64 {
    let mut t = 0i128;
    let mut nt = 1i128;
    let mut r = p as i128;
    let mut nr = (a as i128).abs() % (p as i128);
    let sign = if a < 0 { -1i128 } else { 1i128 };

    while nr != 0 {
        let q = r / nr;
        let tmp = nt;
        nt = t - q * nt;
        t = tmp;
        let tmp_nr = nr;
        nr = r - q * nr;
        r = tmp_nr;
    }

    let mut res = (t * sign) % (p as i128);
    res = (res + p as i128) % (p as i128);
    res as i64
}



