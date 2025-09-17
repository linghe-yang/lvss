fn precompute_barrett_constants(q_prime: i64) {
    assert!(q_prime > 0 && q_prime % 2 != 0); // Assume positive odd prime
    let mult: i128 = ((1i128 << 64) + (q_prime as i128) - 1) / (q_prime as i128);
    let half_q: i64 = q_prime / 2;
    println!("const MULT: i128 = {}; // for q_prime = {}", mult, q_prime);
    println!("const HALF_Q: i64 = {}; // for q_prime = {}", half_q, q_prime);
}

pub(crate) fn symmetric_barrett_reduce(a: i128, q_prime: i64, mult: i128, half_q: i64) -> i64 {
    assert!(q_prime > 0 && q_prime % 2 != 0); // Assume positive odd prime
    let t: i128 = a * mult + (1i128 << 63);
    let quotient: i128 = t >> 64;
    let mut remainder: i128 = a - quotient * (q_prime as i128);

    while remainder > half_q as i128 {
        remainder -= q_prime as i128;
    }
    while remainder < -(half_q as i128) {
        remainder += q_prime as i128;
    }

    remainder as i64
}

#[test]
fn test() {
    let q_prime: i64 = 34359724033;
    // Run precomputation to get constants
    precompute_barrett_constants(q_prime);
}