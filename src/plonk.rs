#![allow(dead_code, unused_imports)]

use std::fmt::Display;

use super::field::Field;
use super::poly::Poly;

pub type F17 = Field<17>;
pub type P17 = Poly<17>;
pub type F101 = Field<101>;

pub fn f101(n: u64) -> F101 {
    F101::from(n)
}

pub fn f17(n: u64) -> F17 {
    F17::from(n)
}
pub fn p17(n: &[i64]) -> P17 {
    P17::from(n)
}
