#![allow(clippy::many_single_char_names)]

use std::{fmt::Display, ops::{Neg, Mul, Add}};
use crate::field::{Field};

pub trait G1:
    Copy
    + Display
    + PartialEq
    + std::hash::Hash
    + PartialOrd
    + Ord
    + Neg<Output = Self>
    + Add<Output = Self>
    + Mul<Self::S, Output = Self>
{
    type F: Field;
    type S: Field;

    fn generator() -> Self;
    fn generator_subgroup_size() -> Self::F;
    fn identity() -> Self;
    fn new(x: Self::F, y: Self::F) -> Self;
    fn x(&self) -> &Self::F;
    fn y(&self) -> &Self::F;
    fn in_curve(&self) -> bool;
    fn is_identity(&self) -> bool;
}

pub trait G2:
    Display + Copy + PartialEq + Neg<Output = Self> + Add + Mul<Self::S, Output = Self>
{
    type F: Field;
    type S: Field;

    fn x(&self) -> &Self::F;
    fn y(&self) -> &Self::F;

    fn new(x: Self::F, y: Self::F) -> Self;
    fn generator() -> Self;
    fn embeeding_degree() -> u64;
}

pub trait GT: Display + Copy + PartialEq + Mul {
    type S;
    fn pow(&self, n: Self::S) -> Self;
}

pub trait Pairing {
    type G1: G1;
    type G2: G2;
    type GT: GT;

    fn pairing(p: Self::G1, q: Self::G2) -> Self::GT;
}