use crate::field::Field;
use crate::poly::Poly;
use itertools::Itertools;
use sha2::{Digest, Sha256};
use std::ops::Index;
use std::{convert::TryInto, marker::PhantomData};

fn merkle_tree<const N: usize, F, H, D, E>(data: D, encode: E, hasher: H) -> [u8; N]
where
    F: Field,
    E: Fn(F) -> Vec<u8>,
    D: IntoIterator<Item = F>,
    H: Fn(&[&[u8]]) -> [u8; N],
{
    let mut v: Vec<[u8; N]> = data
        .into_iter()
        .map(|n| hasher(&[&encode(n)]))
        .collect();

    assert!(v.len().is_power_of_two());

    let mut vv = Vec::with_capacity(v.len() - 1);

    while v.len() > 1 {
        vv.clear();
        for (l, r) in v.iter().tuples() {
            vv.push(hasher(&[l, r]));
        }
        std::mem::swap(&mut v, &mut vv);
    }

    v[0]
}
