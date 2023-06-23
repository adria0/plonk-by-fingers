use super::poly::Poly;
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign}, borrow::Cow,
};

mod u64field;

pub use u64field::U64Field;

pub trait Field:
    Sized
    + Debug
    + Copy
    + Display
    + PartialEq
    + From<i32>
    + From<i64>
    + From<u64>
    + Add<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + AddAssign
    + for<'a> AddAssign<&'a Self>
    + Div<Output = Option<Self>>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> MulAssign<&'a Self>
    + Sub<Output = Self>
    + for<'a> SubAssign<&'a Self>
    + Neg<Output = Self>
{
    type Order: Display;

    fn le_bytes(&self) -> Cow<[u8]>;
    fn order() -> Self::Order;
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn as_u64(&self) -> u64;
    fn in_field(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn pow(&self, exp: u64) -> Self;
    fn as_poly(&self) -> Poly<Self> {
        Poly::new(vec![*self])
    }
    fn carrying_mul(&self, rhs: &Self, carry: &mut Self) -> Self;
}


