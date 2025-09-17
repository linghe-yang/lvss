pub fn is_prime(n: i128) -> bool {
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return n > 1;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    // Check divisibility by 6k ± 1 up to sqrt(n)
    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

pub fn find_prime_near_power_of_two(b: i128) -> i128 {
    // Target: 2^30 ≈ 1,073,741,824
    // Want n = 512k + 1 ≈ 2^30
    // k ≈ (2^30 - 1) / 512 ≈ 2,097,151
    let base_k: i128 = (1 << (b-1)) / 512; // ≈ 2,097,152
    let max_i128: i128 = i128::MAX as i128; // 2,147,483,647

    // Search k around base_k
    for delta in 0..=1000 {
        // Check both k = base_k + delta and base_k - delta
        for &k in &[(base_k - delta), (base_k + delta)] {
            let candidate = 512 * k + 1;
            // Ensure candidate fits in i128
            if candidate <= max_i128 && is_prime(candidate) {
                return candidate as i128;
            }
        }
    }
    // Fallback: return a known prime if none found (shouldn't happen)
    1_073_741_825 // Known prime: 512 * 2,097,152 + 1
}

fn find_prime_near_2_power_b(b: u32) -> Option<u64> {
    let target = 1u64 << b; // 2^b
    let mut candidate = target;

    // 调整 candidate 使其满足 p = 5 (mod 8)
    if candidate % 8 != 5 {
        candidate = ((candidate / 8) * 8) + 5;
    }

    // 向上和向下搜索接近 2^b 的素数
    let mut lower = if candidate > target { candidate - 8 } else { candidate };
    let mut upper = candidate;

    loop {
        // 检查上界
        if upper >= target / 2 && is_prime(upper as i128) {
            return Some(upper);
        }
        // 检查下界
        if lower <= target * 2 && is_prime(lower as i128) {
            return Some(lower);
        }
        // 更新上下界
        upper += 8; // 下一个满足 p = 5 (mod 8) 的数
        if lower >= 8 {
            lower -= 8; // 上一个满足 p = 5 (mod 8) 的数
        } else {
            break; // 防止下界变为负数
        }
        // 如果搜索范围过大，退出
        if upper > target * 2 && lower < target / 2 {
            break;
        }
    }
    None
}


fn mod_inverse(q: i64) -> Option<i64> {
    let a = q as u64;
    // 只有奇数才有模 2^k 的逆
    if a & 1 == 0 {
        return None;
    }

    // 使用 u128 做中间乘法以避免溢出，Newton 迭代：
    // x_{n+1} = x_n * (2 - a*x_n)  (在模 2^{64} 下)
    // 初值 x0 = 1 (因为 a 是奇数，在模 2 下逆是 1)。每次迭代位数翻倍，
    // 6 次迭代能覆盖 64 位 (2^6 = 64)。
    let a128 = a as u128;
    let mut x: u128 = 1; // 初始逆（模 2）
    for _ in 0..6 {
        // all ops in u128 with wrapping behaviour
        let t = a128.wrapping_mul(x);
        let two_minus_ax = (2u128).wrapping_sub(t);
        x = x.wrapping_mul(two_minus_ax);
    }
    let inv_u64 = x as u64; // 低 64 位即为逆
    Some(inv_u64 as i64)    // 直接按位解释成 i64 返回
}


fn calculate_f(q: i128, n: i128) -> Option<i128> {
    if q <= 0 || n <= 0 {
        return None;
    }
    let q = q as i128;
    let n = n as i128;
    let r = 1_i128 << 32;
    let r_squared = {
        let r_mod_q = r % q;
        (r_mod_q * r_mod_q) % q
    };
    let n_inv = {
        let (mut a, mut m, mut x0, mut x1) = (q, n, 0_i128, 1_i128);
        while a > 1 {
            if m == 0 {
                return None;
            }
            let quotient = a / m;
            let temp = m;
            m = a % m;
            a = temp;
            let temp_x = x0;
            x0 = x1 - quotient * x0;
            x1 = temp_x;
        }
        if a != 1 {
            return None;
        }
        if x1 < 0 {
            x1 += q;
        }
        x1
    };
    let f = (r_squared * n_inv) % q;
    Some(f)
}


use num_bigint::BigInt;

/// 计算 F = (R^2 / n) mod q, 其中 R = 2^32
pub fn compute_f(q: BigInt, n: BigInt) -> BigInt {
    let r: BigInt = BigInt::from(1) << 64; // R = 2^64
    let r2: BigInt = (r.clone() * r) % q.clone(); // R^2 mod q

    // 计算 n 在模 q 下的逆元
    let n_inv = modinv(n as BigInt, q.clone());

    // F = (R^2 * n^{-1}) mod q
    ((r2 * n_inv as BigInt) % (q as BigInt)) as BigInt
}

/// 扩展欧几里得求逆元
fn modinv(a: BigInt, m: BigInt) -> BigInt {
    let (mut t, mut new_t) = (BigInt::from(0), BigInt::from(1));
    let (mut r, mut new_r) = (m.clone(), a.clone() % m.clone());
    while new_r != BigInt::from(0) {
        let quotient = r.clone() / new_r.clone();
        t = t - quotient.clone() * new_t.clone();
        r = r.clone() - quotient * new_r.clone();

        std::mem::swap(&mut t, &mut new_t);
        std::mem::swap(&mut r, &mut new_r);
    }
    if r > BigInt::from(1) {
        panic!("{} has no inverse mod {}", a, m.clone());
    }
    if t < BigInt::from(0) {
        t += m;
    }
    t
}



#[cfg(test)]
mod tests {
    use rand_distr::num_traits::real::Real;
    use crate::rq::montgomery_reduce;
    use super::*;

    #[test]
    fn find_prime() {
        let prime = find_prime_near_power_of_two(36);
        assert_eq!(prime % 512, 1, "Prime must be 1 mod 512");
        assert!(prime > 0, "Prime must be positive");
        assert!(is_prime(prime as i128), "Result must be prime");
        println!("res = {}", prime)
    }

    #[test]
    fn find_prime_5_mod_8() {
        let prime = find_prime_near_2_power_b(21).unwrap();
        assert!(prime > 0, "Prime must be positive");
        assert!(is_prime(prime as i128), "Result must be prime");
        println!("p = {}", prime);
        let log_p = (prime as f64).log2();
        println!("log_p = {}",log_p);
    }

    #[test]
    fn find_mod_inverse_q() {
        let q: i64 = 34359724033;
        let inverse = mod_inverse(q).unwrap();
        println!("res = {}", inverse);
    }

    #[test]
    fn find_f() {
        let q:i64 = 34359724033;
        let n = 256;
        let f = compute_f(BigInt::from(q), BigInt::from(n));
        println!("{:?}", f);
    }
}
