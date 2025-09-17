use crate::error::DecryptionError;
use crate::shamir::{shamir_reconstruct, shamir_share};
use crate::util::{cipher_to_bytes, generate_r_matrix, verify_merkle};
use hrcrypto::hash::{do_hash, Hash};
use hrtypes::appxcon::HashingAlg;
use merkle_light::merkle::MerkleTree;
use merkle_light::proof::Proof;
use nalgebra::{DMatrix, DVector};
use std::time::Instant;
use ve::r_ring::R;
use ve::util::random_gaussian_dvector;
use ve::{
    build_b_matrix, calculate_u, gadget_reconstruct, split_m_bar, PublicKey, SecretKey, Store, VE,
};

pub const X_LEN: usize = 2;
pub const Y_LEN: usize = 1;
pub const X_SIGMA: f64 = 1.0;
pub const Y_SIGMA: f64 = 1.0;
pub const R_SIGMA: f64 = 1.0;
pub const B_PRIME: f64 = 100000.0;
pub const BETA: i64 = 128;

pub type ID = i32;

#[derive(Debug, Clone)]
pub struct PrivateShare {
    pub id: ID,
    pub v: DVector<R>,
    pub w: DVector<R>,
    pub merkle_proof: Proof<Hash>,
}

#[derive(Debug, Clone, Default)]
pub struct PublicShare {
    pub merkle_root: Hash,
    pub u_vec: Vec<(ID, DVector<R>)>,
}

#[derive(Debug, Clone)]
pub struct SuppleShare {
    pub id: ID,
    pub v: DVector<R>,
    pub w: DVector<R>,
    pub c: R,
    pub z: DVector<R>,
    pub merkle_proof: Proof<Hash>,
}

#[derive(Debug, Clone)]
pub struct SharingStore {
    pub ciphers: Vec<(ID, DVector<R>, DVector<R>, Store)>,
    pub u_vec: Vec<(ID, DVector<R>)>,
    pub r: DMatrix<R>,
    pub merkle_proofs: Vec<(ID,Proof<Hash>)>,
}

impl Default for SharingStore {
    fn default() -> Self {
        SharingStore {
            ciphers: Vec::default(),
            u_vec: Vec::default(),
            r: DMatrix::from_fn(1,1,|_, _| R::default()),
            merkle_proofs: Vec::default(),
        }
    }
}
impl SharingStore {
    pub fn filter_by_indices(&self, indices: &[ID]) -> Self {
        // 检查 indices 是否为空
        if indices.is_empty() {
            panic!("Input indices array cannot be empty");
        }

        // 收集 ciphers 中匹配的元素
        let filtered_ciphers: Vec<(ID, DVector<R>, DVector<R>, Store)> = indices
            .iter()
            .filter_map(|&idx| {
                self.ciphers
                    .iter()
                    .find(|(i, _, _, _)| *i == idx)
                    .cloned()
            })
            .collect();

        // 收集 u_vec 中匹配的元素
        let filtered_u_vec: Vec<(ID, DVector<R>)> = indices
            .iter()
            .filter_map(|&idx| self.u_vec.iter().find(|(i, _)| *i == idx).cloned())
            .collect();

        let filtered_mfs: Vec<(ID,Proof<Hash>)> = indices.iter().filter_map(|&idx| {
            self.merkle_proofs.iter().find(|(i,_)| *i == idx).cloned()
        }).collect();

        // 检查是否所有 indices 都在 ciphers 和 u_vec 中找到匹配
        if filtered_ciphers.len() != indices.len() || filtered_u_vec.len() != indices.len() || filtered_mfs.len() != indices.len() {
            panic!("One or more indices not found");
        }

        // 返回新的 SharingStore
        SharingStore {
            ciphers: filtered_ciphers,
            u_vec: filtered_u_vec,
            r: self.r.clone(),
            merkle_proofs: filtered_mfs,
        }
    }
}

pub fn share(
    x: DVector<R>,
    n: usize,
    t: usize,
    pks: &Vec<(ID, PublicKey)>,
) -> (Vec<PrivateShare>, PublicShare, SharingStore) {
    let x_shares = shamir_share(&x, n, t);
    let y = random_gaussian_dvector(Y_LEN, Y_SIGMA);
    let y_shares = shamir_share(&y, n, t);

    let mut ciphers = Vec::new();

    for (id, xi) in x_shares.iter() {
        let pk = &pks.iter().find(|(i, pk)| i == id).unwrap().1;
        let yi = &y_shares[*id as usize - 1].1;
        let (v, w, st) = VE::encrypt(pk, xi, yi);
        ciphers.push((*id, v, w, st));
    }
    let mut leafs = Vec::new();
    for ele in ciphers.iter() {
        leafs.push(do_hash(&cipher_to_bytes((&ele.1, &ele.2))));
    }
    let mt = MerkleTree::<Hash, HashingAlg>::from_iter(leafs);
    let h = mt.root();

    let r = generate_r_matrix(h, X_LEN, Y_LEN, R_SIGMA);

    let mut u_vec = Vec::new();

    for (id, xi) in x_shares.iter() {
        let yi = &y_shares[*id as usize - 1].1;
        let ui = calculate_u(&r, xi, yi);
        u_vec.push((*id, ui));
    }
    let mut shares = Vec::new();
    let mut mfs = Vec::new();
    for (id, _) in x_shares.iter() {
        let cipher = &ciphers.iter().find(|(i, _, _, _)| i == id).unwrap();
        let mf = mt.gen_proof(*id as usize - 1);
        mfs.push((*id,mf.clone()));
        let share = PrivateShare {
            id: *id,
            v: cipher.1.clone(),
            w: cipher.2.clone(),
            merkle_proof: mf,
        };
        shares.push(share);
    }
    let ps = PublicShare {
        merkle_root: h,
        u_vec: u_vec.clone(),
    };

    let st = SharingStore {
        ciphers,
        u_vec,
        merkle_proofs: mfs,
        r,
    };

    (shares, ps, st)
}

pub fn supple_share(vss_st: SharingStore, pks: &Vec<(ID, PublicKey)>) -> Vec<SuppleShare> {
    let b = build_b_matrix(&vss_st.r);
    let mut shares = Vec::new();
    for (id, v, w, st) in vss_st.ciphers.iter() {
        let ui = &vss_st.u_vec.iter().find(|(i, _)| i == id).unwrap().1;
        let pk = &pks.iter().find(|(i, _)| i == id).unwrap().1;
        let (c, z) = VE::prove(pk, st, &b, &ui, v, w);
        let mf = &vss_st.merkle_proofs.iter().find(|(i, _)| i == id).unwrap().1;
        let share = SuppleShare {
            id: *id,
            v: v.clone(),
            w: w.clone(),
            c,
            z,
            merkle_proof: mf.clone(),
        };
        shares.push(share);
    }
    shares
}

pub fn is_merkle_valid(h: Hash, v: &DVector<R>, w: &DVector<R>, proof: &Proof<Hash>) -> bool {
    let leaf = do_hash(&cipher_to_bytes((v, w)));
    if !verify_merkle(&h, leaf, &proof, &mut HashingAlg::new()) {
        return false;
    }
    true
}

pub fn verify(pk: &PublicKey, share: &SuppleShare, u: &DVector<R>, h: Hash) -> bool {
    if !is_merkle_valid(h, &share.v, &share.w, &share.merkle_proof){
        return false;
    }
    let r = generate_r_matrix(h, X_LEN, Y_LEN, R_SIGMA);
    let b = build_b_matrix(&r);
    VE::verify(pk, &b, u, (&share.v, &share.w, &share.c, &share.z))
}

pub fn decrypt(
    sk: &SecretKey,
    share: &PrivateShare,
) -> Result<(DVector<R>, DVector<R>), DecryptionError> {
    let m = VE::decrypt(sk, &share.v, &share.w);
    let (sx, sy) = split_m_bar(m, X_LEN, Y_LEN).map_err(|e| DecryptionError::InvalidLength)?;
    let x_recover = gadget_reconstruct(&sx);
    let y_recover = gadget_reconstruct(&sy);
    Ok((x_recover, y_recover))
}

pub fn try_decrypt(pk: &PublicKey, sk: &SecretKey, share: &SuppleShare) -> Result<(DVector<R>, DVector<R>), DecryptionError> {
    let res = VE::try_decrypt(pk, sk, (&share.v, &share.w, &share.c, &share.z));
    match res {
        Some(m_bar) => {
            let (sx, sy) =
                split_m_bar(m_bar, X_LEN, Y_LEN).map_err(|e| DecryptionError::InvalidLength)?;
            let x_recover = gadget_reconstruct(&sx);
            let y_recover = gadget_reconstruct(&sy);
            Ok((x_recover, y_recover))
        }
        None => Err(DecryptionError::InvalidCipher),
    }
}
pub fn test_share_avsss() {
    let n = 10;
    let t = 3;
    let x = random_gaussian_dvector(X_LEN, X_SIGMA);
    let mut pks = Vec::new();
    let mut sks = Vec::new();
    for i in 1..=n {
        let (pk, sk) = VE::gen_keypair();
        pks.push((i as ID, pk));
        sks.push((i as ID, sk));
    }

    let now = Instant::now();

    let (pri_shares, pub_share, st) = share(x.clone(), n, t, &pks);

    let duration = now.elapsed();
    println!("share took {} ms", duration.as_millis());

    let u = shamir_reconstruct(&pub_share.u_vec, t).unwrap();

    let vss_st = st.filter_by_indices(&[1]);

    let sups = supple_share(vss_st, &pks);
    let vs = verify(&pks[0].1, &sups[0], &pub_share.u_vec[0].1, pub_share.merkle_root);
    println!("verify res = {}", vs);

    let mut x_decs = Vec::new();
    for id in 0..n {
        let dec = decrypt(&sks[id].1, &pri_shares[id]).unwrap();
        x_decs.push((id as ID + 1, dec.0));
    }
    let x_recon = shamir_reconstruct(&x_decs, t).unwrap();

    assert_eq!(x, x_recon);
}


