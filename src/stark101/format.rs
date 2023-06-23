use crate::field::Field;
use crate::mulmodg::MulGroupMod;
use crate::poly::Poly;
use itertools::Itertools;
use num_bigint::BigUint;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::ops::Index;
use std::{convert::TryInto, marker::PhantomData};

use super::FF;

pub fn sha256hex(data: String) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data.as_bytes());
    hex::encode(hasher.finalize())
}
pub fn as_neg_str(v: FF) -> String {
    let m1 = FF::zero() - FF::one();
    let limit = (m1 / 2.into()).unwrap() + FF::one();
    if v > limit {
        let mv = m1 - v + FF::one();
        format!("-{}", mv)
    } else {
        format!("{}", v)
    }
}