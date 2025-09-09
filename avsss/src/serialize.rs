use crate::components::{PrivateShare, PublicShare, SuppleShare};
use hrcrypto::hash::Hash;
use merkle_light::proof::Proof;
use nalgebra::DVector;
use serde::de::{Error, MapAccess, Visitor};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use ve::r_ring::R;

impl Serialize for PrivateShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut st = serializer.serialize_struct("PrivateShare", 4)?;
        st.serialize_field("id", &self.id)?;

        // DVector -> Vec
        st.serialize_field("v", &self.v.as_slice())?;
        st.serialize_field("w", &self.w.as_slice())?;

        // Proof 拆分
        st.serialize_field("lemma", &self.merkle_proof.lemma())?;
        st.serialize_field("path", &self.merkle_proof.path())?;

        st.end()
    }
}

struct PrivateShareVisitor;

impl<'de> Visitor<'de> for PrivateShareVisitor {
    type Value = PrivateShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct PrivateShare")
    }

    fn visit_map<A>(self, mut map: A) -> Result<PrivateShare, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut id = None;
        let mut v: Option<Vec<R>> = None;
        let mut w: Option<Vec<R>> = None;
        let mut lemma: Option<Vec<Hash>> = None;
        let mut path: Option<Vec<bool>> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "id" => id = Some(map.next_value()?),
                "v" => v = Some(map.next_value()?),
                "w" => w = Some(map.next_value()?),
                "lemma" => lemma = Some(map.next_value()?),
                "path" => path = Some(map.next_value()?),
                _ => {
                    return Err(A::Error::unknown_field(
                        &key,
                        &["id", "v", "w", "lemma", "path"],
                    ));
                }
            }
        }

        let id = id.ok_or_else(|| A::Error::missing_field("id"))?;
        let v = v.ok_or_else(|| A::Error::missing_field("v"))?;
        let w = w.ok_or_else(|| A::Error::missing_field("w"))?;
        let lemma = lemma.ok_or_else(|| A::Error::missing_field("lemma"))?;
        let path = path.ok_or_else(|| A::Error::missing_field("path"))?;

        Ok(PrivateShare {
            id,
            v: DVector::from_vec(v),
            w: DVector::from_vec(w),
            merkle_proof: Proof::new(lemma, path),
        })
    }
}

impl<'de> Deserialize<'de> for PrivateShare {
    fn deserialize<D>(deserializer: D) -> Result<PrivateShare, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct(
            "PrivateShare",
            &["id", "v", "w", "lemma", "path"],
            PrivateShareVisitor,
        )
    }
}

impl Serialize for SuppleShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut st = serializer.serialize_struct("SuppleShare", 6)?;
        st.serialize_field("id", &self.id)?;

        // DVector -> Vec
        st.serialize_field("v", &self.v.as_slice())?;
        st.serialize_field("w", &self.w.as_slice())?;
        st.serialize_field("c", &self.c)?;
        st.serialize_field("z", &self.z.as_slice())?;

        // Proof 拆分
        st.serialize_field("lemma", &self.merkle_proof.lemma())?;
        st.serialize_field("path", &self.merkle_proof.path())?;

        st.end()
    }
}

struct SuppleShareVisitor;

impl<'de> Visitor<'de> for SuppleShareVisitor {
    type Value = SuppleShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct SuppleShare")
    }

    fn visit_map<A>(self, mut map: A) -> Result<SuppleShare, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut id = None;
        let mut v: Option<Vec<R>> = None;
        let mut w: Option<Vec<R>> = None;
        let mut c = None;
        let mut z: Option<Vec<R>> = None;
        let mut lemma: Option<Vec<Hash>> = None;
        let mut path: Option<Vec<bool>> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "id" => id = Some(map.next_value()?),
                "v" => v = Some(map.next_value()?),
                "w" => w = Some(map.next_value()?),
                "c" => c = Some(map.next_value()?),
                "z" => z = Some(map.next_value()?),
                "lemma" => lemma = Some(map.next_value()?),
                "path" => path = Some(map.next_value()?),
                _ => {
                    return Err(A::Error::unknown_field(
                        &key,
                        &["id", "v", "w", "c", "z", "lemma", "path"],
                    ));
                }
            }
        }

        let id = id.ok_or_else(|| A::Error::missing_field("id"))?;
        let v = v.ok_or_else(|| A::Error::missing_field("v"))?;
        let w = w.ok_or_else(|| A::Error::missing_field("w"))?;
        let c = c.ok_or_else(|| A::Error::missing_field("c"))?;
        let z = z.ok_or_else(|| A::Error::missing_field("z"))?;
        let lemma = lemma.ok_or_else(|| A::Error::missing_field("lemma"))?;
        let path = path.ok_or_else(|| A::Error::missing_field("path"))?;

        Ok(SuppleShare {
            id,
            v: DVector::from_vec(v),
            w: DVector::from_vec(w),
            c,
            z: DVector::from_vec(z),
            merkle_proof: Proof::new(lemma, path),
        })
    }
}

impl<'de> Deserialize<'de> for SuppleShare {
    fn deserialize<D>(deserializer: D) -> Result<SuppleShare, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct(
            "SuppleShare",
            &["id", "v", "w", "c", "z", "lemma", "path"],
            SuppleShareVisitor,
        )
    }
}

impl Serialize for PublicShare {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut st = serializer.serialize_struct("PublicShare", 2)?;
        st.serialize_field("merkle_root", &self.merkle_root)?;

        // 把 Vec<(i32, DVector<R>)> 转成 Vec<(i32, Vec<R>)>
        let v_serializable: Vec<(i32, Vec<R>)> = self
            .u_vec
            .iter()
            .map(|(id, dv)| (*id, dv.as_slice().to_vec()))
            .collect();

        st.serialize_field("u_vec", &v_serializable)?;
        st.end()
    }
}

struct PublicShareVisitor;

impl<'de> Visitor<'de> for PublicShareVisitor {
    type Value = PublicShare;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct PublicShare")
    }

    fn visit_map<A>(self, mut map: A) -> Result<PublicShare, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut merkle_root: Option<Hash> = None;
        let mut u_vec: Option<Vec<(i32, Vec<R>)>> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "merkle_root" => merkle_root = Some(map.next_value()?),
                "u_vec" => u_vec = Some(map.next_value()?),
                _ => return Err(A::Error::unknown_field(&key, &["merkle_root", "u_vec"])),
            }
        }

        let merkle_root = merkle_root.ok_or_else(|| A::Error::missing_field("merkle_root"))?;
        let u_vec_raw = u_vec.ok_or_else(|| A::Error::missing_field("u_vec"))?;

        // Vec<(i32, Vec<R>)> → Vec<(i32, DVector<R>)>
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
        deserializer.deserialize_struct(
            "PublicShare",
            &["merkle_root", "u_vec"],
            PublicShareVisitor,
        )
    }
}
