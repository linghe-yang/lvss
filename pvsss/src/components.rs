use crate::shamir::{shamir_reconstruct, shamir_share};
use crate::util::{cipher_to_bytes, generate_r_matrix, verify_merkle};
use hrcrypto::hash::{do_hash, Hash};
use hrtypes::appxcon::HashingAlg;
use merkle_light::merkle::MerkleTree;
use merkle_light::proof::Proof;
use nalgebra::DVector;
use std::time::Instant;
use ve::r_ring::R;
use ve::util::{random_gaussian_dvector, vector_norm_inf};
use ve::{build_b_matrix, calculate_u, gadget_reconstruct, split_m_bar, PublicKey, SecretKey, VE};

pub const X_LEN: usize = 1;
pub const Y_LEN: usize = 1;
pub const X_SIGMA: f64 = 1.0;
pub const Y_SIGMA: f64 = 1.0;
pub const R_SIGMA: f64 = 1.0;

#[derive(Debug, Clone)]
pub struct PrivateShare {
    pub id: i32,
    pub v: DVector<R>,
    pub w: DVector<R>,
    pub c: R,
    pub z: DVector<R>,
    pub merkle_proof: Proof<Hash>,
}

#[derive(Debug, Clone)]
pub struct PublicShare {
    pub merkle_root: Hash,
    pub u_vec: Vec<(i32, DVector<R>)>,
}

pub fn share(
    x: DVector<R>,
    n: usize,
    t: usize,
    pks: &Vec<(i32, PublicKey)>,
) -> (Vec<PrivateShare>, PublicShare) {
    let x_shares = shamir_share(&x, n, t);
    let y = random_gaussian_dvector(Y_LEN, Y_SIGMA);
    let y_shares = shamir_share(&y, n, t);

    let mut ciphers = Vec::new();

    for (id, xi) in x_shares.iter() {
        let pk = &pks.iter().find(|(i, pk)| i == id).unwrap().1;
        let yi = &y_shares[*id as usize - 1].1;
        let (v, w, st) = VE::encrypt(pk, xi, yi);
        ciphers.push((id, v, w, st));
    }
    let mut leafs = Vec::new();
    for ele in ciphers.iter() {
        leafs.push(do_hash(&cipher_to_bytes((&ele.1, &ele.2))));
    }
    let mt = MerkleTree::<Hash, HashingAlg>::from_iter(leafs);
    let h = mt.root();

    let r = generate_r_matrix(h, X_LEN, Y_LEN, R_SIGMA);

    let mut u_vec = Vec::new();

    let b = build_b_matrix(&r);

    for (id, xi) in x_shares.iter() {
        let yi = &y_shares[*id as usize - 1].1;
        let ui = calculate_u(&r, xi, yi);
        u_vec.push((*id, ui));
    }
    let mut shares = Vec::new();
    for (id, _) in x_shares.iter() {
        let ui = &u_vec.iter().find(|(i, _)| i == id).unwrap().1;
        let pk = &pks.iter().find(|(i, _)| i == id).unwrap().1;
        let cipher = &ciphers.iter().find(|(i, _, _, _)| *i == id).unwrap();
        let (c, z) = VE::prove(pk, &cipher.3, &b, &ui, &cipher.1, &cipher.2);
        let mf = mt.gen_proof(*id as usize - 1);
        let share = PrivateShare {
            id: *id,
            v: cipher.1.clone(),
            w: cipher.2.clone(),
            c,
            z,
            merkle_proof: mf,
        };
        shares.push(share);
    }
    let ps = PublicShare {
        merkle_root: h,
        u_vec,
    };

    (shares, ps)
}

pub fn verify(pk: &PublicKey, share: &PrivateShare, u: &DVector<R>, h: Hash) -> bool {
    let leaf = do_hash(&cipher_to_bytes((&share.v, &share.w)));
    if !verify_merkle(&h, leaf, &share.merkle_proof, &mut HashingAlg::new()) {
        println!("Failed to verify merkle proof");
        return false;
    }
    let r = generate_r_matrix(h, X_LEN, Y_LEN, R_SIGMA);
    let b = build_b_matrix(&r);
    VE::verify(pk, &b, u, (&share.v, &share.w, &share.c, &share.z))
}

pub fn decrypt(pk: &PublicKey, sk: &SecretKey, share: &PrivateShare) -> Option<(DVector<R>,DVector<R>)> {
    let res = VE::decrypt(pk,sk,(&share.v, &share.w, &share.c, &share.z));
    match res {
        Some(m_bar) => {
            let (sx,sy) = split_m_bar(m_bar, X_LEN, Y_LEN);
            let x_recover = gadget_reconstruct(&sx);
            let y_recover = gadget_reconstruct(&sy);
            Some((x_recover, y_recover))
        }
        None => {
            None
        }
    }
}
pub fn test_share() {
    let n = 10;
    let t = 3;
    let x = random_gaussian_dvector(X_LEN, X_SIGMA);
    println!("x = {:?}", x);
    let mut pks = Vec::new();
    let mut sks = Vec::new();
    for i in 1..=n {
        let (pk, sk) = VE::gen_keypair();
        pks.push((i as i32, pk));
        sks.push((i as i32, sk));
    }

    let now = Instant::now();

    let shares = share(x.clone(), n, t, &pks);
    let duration = now.elapsed();
    println!("share took {} ms", duration.as_millis());

    let ps = &shares.0[0];
    let pps = &shares.1;
    let pk = &pks[0].1;
    let u = shamir_reconstruct(&pps.u_vec, t).unwrap();

    let res = verify(pk, ps, &pps.u_vec[0].1, pps.merkle_root);
    let u_norm = vector_norm_inf(&u);
    println!("u's norm: {:?}", u_norm);
    println!("verify result: {:?}", res);

    let mut x_decs = Vec::new();
    for id in 0..n{
        let dec = decrypt(&pks[id].1, &sks[id].1, &shares.0[id]).unwrap();
        x_decs.push((id as i32 + 1, dec.0));
    }
    let x_recon = shamir_reconstruct(&x_decs, t).unwrap();

    println!("x_recon: {:?}", x_recon);
    assert_eq!(x,x_recon);
}
