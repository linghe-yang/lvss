use crate::components::{PrivateShare, PublicShare, SuppleShare};
use hrcrypto::hash::Hash;
use merkle_light::proof::Proof;
use nalgebra::DVector;
use serde::de::{self, SeqAccess, Visitor};
use serde::ser::{SerializeSeq};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use ve::r_ring::R;

impl Serialize for PrivateShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(5))?;
        seq.serialize_element(&self.id)?;
        seq.serialize_element(&self.v.as_slice())?;
        seq.serialize_element(&self.w.as_slice())?;
        seq.serialize_element(&self.merkle_proof.lemma())?;
        seq.serialize_element(&self.merkle_proof.path())?;
        seq.end()
    }
}

struct PrivateShareVisitor;

impl<'de> Visitor<'de> for PrivateShareVisitor {
    type Value = PrivateShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct PrivateShare as sequence")
    }

    fn visit_seq<A>(self, mut seq: A) -> Result<PrivateShare, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let id: i32 = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &self))?;
        let v_vec: Vec<R> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &self))?;
        let w_vec: Vec<R> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(2, &self))?;
        let lemma: Vec<Hash> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(3, &self))?;
        let path: Vec<bool> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(4, &self))?;

        Ok(PrivateShare {
            id,
            v: DVector::from_vec(v_vec),
            w: DVector::from_vec(w_vec),
            merkle_proof: Proof::new(lemma, path),
        })
    }
}

impl<'de> Deserialize<'de> for PrivateShare {
    fn deserialize<D>(deserializer: D) -> Result<PrivateShare, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_seq(PrivateShareVisitor)
    }
}

impl Serialize for SuppleShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(7))?;
        seq.serialize_element(&self.id)?;
        seq.serialize_element(&self.v.as_slice())?;
        seq.serialize_element(&self.w.as_slice())?;
        seq.serialize_element(&self.c)?;
        seq.serialize_element(&self.z.as_slice())?;
        seq.serialize_element(&self.merkle_proof.lemma())?;
        seq.serialize_element(&self.merkle_proof.path())?;
        seq.end()
    }
}

struct SuppleShareVisitor;

impl<'de> Visitor<'de> for SuppleShareVisitor {
    type Value = SuppleShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct SuppleShare as sequence")
    }

    fn visit_seq<A>(self, mut seq: A) -> Result<SuppleShare, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let id: i32 = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &self))?;
        let v_vec: Vec<R> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &self))?;
        let w_vec: Vec<R> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(2, &self))?;
        let c: R = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(3, &self))?;
        let z_vec: Vec<R> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(4, &self))?;
        let lemma: Vec<Hash> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(5, &self))?;
        let path: Vec<bool> = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(6, &self))?;

        Ok(SuppleShare {
            id,
            v: DVector::from_vec(v_vec),
            w: DVector::from_vec(w_vec),
            c,
            z: DVector::from_vec(z_vec),
            merkle_proof: Proof::new(lemma, path),
        })
    }
}

impl<'de> Deserialize<'de> for SuppleShare {
    fn deserialize<D>(deserializer: D) -> Result<SuppleShare, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_seq(SuppleShareVisitor)
    }
}

impl Serialize for PublicShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(2))?;
        seq.serialize_element(&self.merkle_root)?;
        let u_vec_serializable: Vec<(i32, Vec<R>)> = self
            .u_vec
            .iter()
            .map(|(id, dv)| (*id, dv.as_slice().to_vec()))
            .collect();
        seq.serialize_element(&u_vec_serializable)?;

        seq.end()
    }
}

struct PublicShareVisitor;

impl<'de> Visitor<'de> for PublicShareVisitor {
    type Value = PublicShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct PublicShare as sequence")
    }

    fn visit_seq<A>(self, mut seq: A) -> Result<PublicShare, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let merkle_root: Hash = seq
            .next_element()?
            .ok_or_else(|| de::Error::custom("Missing merkle_root in sequence"))?;
        let u_vec_raw: Vec<(i32, Vec<R>)> = seq
            .next_element()?
            .ok_or_else(|| de::Error::custom("Missing u_vec in sequence"))?;
        if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
            return Err(de::Error::custom("Too many elements in sequence"));
        }
        let u_vec: Vec<(i32, DVector<R>)> = u_vec_raw
            .into_iter()
            .map(|(id, v)| (id, DVector::from_vec(v)))
            .collect();

        Ok(PublicShare { merkle_root, u_vec })
    }
}

impl<'de> Deserialize<'de> for PublicShare {
    fn deserialize<D>(deserializer: D) -> Result<PublicShare, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_seq(PublicShareVisitor)
    }
}