use ve::r_ring::{N, R};
use ve::rp::P_PRIME;

// 模 p 下的加法
fn mod_add(a: i32, b: i32, p: i32) -> i32 {
    let mut sum = (a + b) % p;
    if sum < 0 {
        sum += p;
    }
    sum
}

// 模 p 下的乘法
fn mod_mul(a: i32, b: i32, p: i32) -> i32 {
    let mut prod = ((a as i64 * b as i64) % p as i64) as i32;
    if prod < 0 {
        prod += p;
    }
    prod
}

// 模 p 下的模逆（使用扩展欧几里得算法）
fn mod_inverse(a: i32, p: i32) -> Option<i32> {
    let (mut r0, mut r1) = (a, p);
    let (mut t0, mut t1) = (0, 1);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    if r0 != 1 {
        return None; // 模逆不存在
    }
    if t0 < 0 {
        t0 += p;
    }
    Some(t0)
}

// Lagrange 插值在模 p 下计算 f(0)
fn lagrange_interpolate(points: &[(i32, i32)], p: i32) -> Option<i32> {
    let n = points.len();
    let mut result = 0;

    for i in 0..n {
        let (xi, yi) = points[i];
        let mut li = 1; // Lagrange 基多项式 l_i(0)
        for j in 0..n {
            if i != j {
                let xj = points[j].0;
                let numerator = mod_neg(xj, p); // -x_j
                let denominator = mod_sub(xi, xj, p); // x_i - x_j
                let inv_denominator = mod_inverse(denominator, p)?;
                li = mod_mul(li, mod_mul(numerator, inv_denominator, p), p);
            }
        }
        result = mod_add(result, mod_mul(yi, li, p), p);
    }
    Some(result)
}

// 模 p 下的减法
fn mod_sub(a: i32, b: i32, p: i32) -> i32 {
    mod_add(a, -b, p)
}

// 模 p 下的负数
fn mod_neg(a: i32, p: i32) -> i32 {
    if a == 0 {
        0
    } else {
        p - a
    }
}

// Shamir 重构函数
pub fn shamir_reconstruct_r(shares: &[(i32, R)], t: usize) -> Option<R> {
    let p = P_PRIME;
    // 检查输入份额数量
    if shares.len() < t + 1 {
        return None;
    }

    // 初始化结果多项式
    let mut result = R { coeffs: [0; N] };

    // 对每个系数进行 Lagrange 插值
    for coeff_idx in 0..N {
        // 提取第 coeff_idx 个系数的点 (x_i, y_i)
        let points: Vec<(i32, i32)> = shares
            .iter()
            .take(t + 1)
            .map(|&(x, r)| (x, r.coeffs[coeff_idx]))
            .collect();

        // 计算 f(0) 即秘密的第 coeff_idx 个系数
        result.coeffs[coeff_idx] = lagrange_interpolate(&points, p)?;
    }

    // 如果份额数量大于 t+1，验证剩余份额
    if shares.len() > t + 1 {
        for &(x, ref r) in shares.iter().skip(t + 1) {
            // 对每个剩余份额，重新计算 f(x) 并比较
            let mut computed_r = R { coeffs: [0; N] };
            for coeff_idx in 0..N {
                let points: Vec<(i32, i32)> = shares
                    .iter()
                    .take(t + 1)
                    .map(|&(x_i, r_i)| (x_i, r_i.coeffs[coeff_idx]))
                    .collect();
                let mut fx = 0;
                for i in 0..points.len() {
                    let (xi, yi) = points[i];
                    let mut li = 1;
                    for j in 0..points.len() {
                        if i != j {
                            let xj = points[j].0;
                            let numerator = mod_sub(x, xj, p);
                            let denominator = mod_sub(xi, xj, p);
                            let inv_denominator = mod_inverse(denominator, p)?;
                            li = mod_mul(li, mod_mul(numerator, inv_denominator, p), p);
                        }
                    }
                    fx = mod_add(fx, mod_mul(yi, li, p), p);
                }
                computed_r.coeffs[coeff_idx] = fx;
            }
            if computed_r != *r {
                return None; // 验证失败
            }
        }
    }

    Some(result)
}