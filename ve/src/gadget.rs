use crate::r_ring::{N, R};
use crate::rp::{LOG_P, P_PRIME};
use crate::util::*;
use nalgebra::{DMatrix, DVector};
use rand_distr::num_traits::Zero;

pub fn gadget_decompose(f: &DVector<R>) -> DVector<R> {
    let len = f.len();
    let mut s: DVector<R> = DVector::zeros(len * LOG_P);

    for i in 0..len {
        let mut poly = f[i];
        poly.caddp();

        for j in 0..N {
            let mut c = poly.coeffs[j];

            for k in 0..LOG_P {
                let u_k = c % 2;
                s[i * LOG_P + k].coeffs[j] = u_k;
                c /= 2;
            }
            if c != 0 {
                panic!("gadget error");
            }
        }
    }

    s
}

pub fn gadget_reconstruct(s: &DVector<R>) -> DVector<R> {
    let len_s = s.len();
    let len_f = len_s / LOG_P;
    let mut f:DVector<R> = DVector::zeros(len_f);

    for i in 0..len_f {
        let mut sum = R::zero();
        for k in 0..LOG_P {
            let s_ki = s[i * LOG_P + k];
            let term = s_ki * (1 << k);
            sum += term;
        }
        f[i] = sum;
        f[i].mod_p();
    }
    f
}




pub fn calculate_u(r: &DMatrix<R>, x: &DVector<R>, y: &DVector<R>) -> DVector<R> {
    let rx = matrix_vector_mul_mod_p(&r, &x);
    let v = vector_add_mod_p(&rx, &y);
    v
}


/// Builds the matrix [R * G_x | G_y] where R is yl * xl matrix, G_x and G_y are constructed based on R's dimensions.
/// Elements are of type R (polynomial ring).
pub fn build_b_matrix(r: &DMatrix<R>) -> DMatrix<R> {
    let (yl, xl) = r.shape();
    let log_p = LOG_P; // Assuming LOG_P from previous context

    // Construct G_x: xl * (xl * log_p) matrix with constant elements
    let mut g_x = DMatrix::zeros(xl, xl * log_p);
    for i in 0..xl {
        for j in 0..log_p {
            g_x[(i, i * log_p + j)] = R::from(1 << j); // g^T = (1, 2, 4, ..., 2^log_p)
        }
    }

    // Construct G_y: yl * (yl * log_p) matrix with constant elements
    let mut g_y = DMatrix::zeros(yl, yl * log_p);
    for i in 0..yl {
        for j in 0..log_p {
            g_y[(i, i * log_p + j)] = R::from(1 << j); // g^T = (1, 2, 4, ..., 2^log_p)
        }
    }

    // Manually compute R * G_x
    let mut rg_x = DMatrix::zeros(yl, xl * log_p);
    for i in 0..yl {
        for j in 0..(xl * log_p) {
            let mut sum = R::zero();
            for k in 0..xl {
                sum = sum.add_mod_p(&r[(i, k)].mul_mod_p(&g_x[(k, j)]));
            }
            rg_x[(i, j)] = sum;
        }
    }

    // Horizontally concatenate [R * G_x | G_y]
    let mut result = DMatrix::zeros(yl, xl * log_p + yl * log_p);
    result.index_mut((0..yl, 0..(xl * log_p))).copy_from(&rg_x);
    result.index_mut((0..yl, xl * log_p..(xl * log_p + yl * log_p))).copy_from(&g_y);

    result
}

pub fn split_m_bar(vector: DVector<R>, xl: usize, yl: usize) -> (DVector<R>, DVector<R>) {
    let xl_len = xl * LOG_P;
    let yl_len = yl * LOG_P;

    assert_eq!(vector.len(), xl_len + yl_len, "m_bar's length must equal xl * log_p + yl * log_p");

    let xl_vector = DVector::from_vec(vector.as_slice()[..xl_len].to_vec());
    let yl_vector = DVector::from_vec(vector.as_slice()[xl_len..].to_vec());

    (xl_vector, yl_vector)
}

#[test]
fn verify_bmu_relation() {
    // Generate random x, y with coefficients in [-p/2, p/2]
    let xl = 2;
    let yl = 2;

    let half_p = P_PRIME / 2;
    let r = DMatrix::from_fn(yl,xl,|_,_|R::random_uniform(-half_p,half_p));

    let x_vec = random_dvector(xl);
    let y_vec = random_dvector(yl);

    let rx = matrix_vector_mul_mod_p(&r, &x_vec);
    let v = vector_add_mod_p(&rx, &y_vec);
    println!("v = {:?}", v);


    // Compute s_x, s_y using gadget_decompose
    let s_x = gadget_decompose(&x_vec);
    let s_y = gadget_decompose(&y_vec);

    // Concatenate s_x and s_y into m
    let m = concat_vectors(&[s_x, s_y]);

    // Compute B = [R * G_x | G_y]
    let b = build_b_matrix(&r);

    // Compute u = Bm mod p
    let u = matrix_vector_mul_mod_p(&b, &m);
    println!("u = {:?}", u);

    // Verify u = v
    let res = u.iter().zip(v.iter()).all(|(u_i, v_i)| u_i == v_i);
    println!("res = {}", res);
}