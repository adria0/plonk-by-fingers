use crate::{plonk::PlonkTypes, poly::Field};


use crate::field::{U64Field};

use super::pairing::PBHPairing;
use super::{g1, g2, gt};

pub type F101 = U64Field<101>;
pub const fn f101(x: u64) -> F101 {
    U64Field::<101>(x % 101)
}

pub type F17 = U64Field<17>;
pub const fn f17(x: u64) -> F17 {
    U64Field::<17>(x % 17)
}
#[derive(Debug, PartialEq)]
pub struct PlonkByHandTypes {}
impl PlonkTypes for PlonkByHandTypes {
    type G1 = g1::G1P;
    type G2 = g2::G2P;
    type GT = gt::GTP;
    type E = PBHPairing;
    type GF = F101;
    type HF = F17;
    const K1: Self::HF = f17(2);
    const K2: Self::HF = f17(3);
    const OMEGA: Self::HF = f17(4);
    fn gf(sg: F17) -> F101 {
        F101::from(sg.as_u64())
    }
}