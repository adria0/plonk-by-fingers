use super::Field;
use crate::poly::Poly;

use std::{
    borrow::Cow,
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    let (mut s, mut old_s) = (0, 1);
    let (mut t, mut old_t) = (1, 0);
    let (mut r, mut old_r) = (b, a);

    while r != 0 {
        let quotient = old_r / r;
        old_r -= quotient * r;
        std::mem::swap(&mut old_r, &mut r);
        old_s -= quotient * s;
        std::mem::swap(&mut old_s, &mut s);
        old_t -= quotient * t;
        std::mem::swap(&mut old_t, &mut t);
    }
    (old_r, old_s, old_t)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct U64Field<const M: u64>(pub(crate) u64);

#[allow(non_snake_case)]
impl<const M: u64> Field for U64Field<M> {
    type Order = u64;

    fn order() -> Self::Order {
        M
    }
    fn zero() -> Self {
        Self::from(0u64)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    fn one() -> Self {
        Self::from(1u64)
    }
    fn as_u64(&self) -> u64 {
        self.0
    }
    fn in_field(&self) -> bool {
        self.0 < M
    }
    fn inv(&self) -> Option<Self> {
        let (gcd, c, _) = extended_gcd(self.0 as i64, M as i64);
        if gcd == 1 {
            if c < 0 {
                Some(Self((M as i64 + c) as u64))
            } else {
                Some(Self(c as u64))
            }
        } else {
            None
        }
    }
    fn pow(&self, mut exp: u64) -> Self {
        let mut result = Self::one();
        let mut base = *self;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            exp >>= 1;
            base = base * base;
        }
        result
    }
    fn carrying_mul(&self, rhs: &Self, carry: &mut Self) -> Self {
        let r = self.0 as u128 * rhs.0 as u128 + carry.0 as u128;

        *carry = Self((r / M as u128) as u64);

        Self((r % M as u128) as u64)
    }
    fn le_bytes(&self) -> Cow<[u8]> {
        Cow::Owned(self.0.to_le_bytes().to_vec())
    }
}

impl<const M: u64> From<i64> for U64Field<M> {
    fn from(n: i64) -> Self {
        if n < 0 {
            -Self::from(-n as u64)
        } else {
            Self::from(n as u64)
        }
    }
}

impl<const M: u64> From<i32> for U64Field<M> {
    fn from(n: i32) -> Self {
        if n < 0 {
            -Self::from(-n as u64)
        } else {
            Self::from(n as u64)
        }
    }
}

impl<const M: u64> From<u64> for U64Field<M> {
    fn from(n: u64) -> Self {
        Self(n % M)
    }
}

impl<const M: u64> Display for U64Field<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const M: u64> Add for U64Field<M> {
    type Output = U64Field<M>;
    fn add(self, rhs: Self) -> Self::Output {
        U64Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<&U64Field<M>> for U64Field<M> {
    type Output = U64Field<M>;
    fn add(self, rhs: &U64Field<M>) -> Self::Output {
        U64Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<U64Field<M>> for &U64Field<M> {
    type Output = U64Field<M>;
    fn add(self, rhs: U64Field<M>) -> Self::Output {
        U64Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<&U64Field<M>> for &U64Field<M> {
    type Output = U64Field<M>;
    fn add(self, rhs: &U64Field<M>) -> Self::Output {
        U64Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> AddAssign<&U64Field<M>> for U64Field<M> {
    fn add_assign(&mut self, rhs: &U64Field<M>) {
        self.0 = (self.0 + rhs.0) % M;
    }
}

impl<const M: u64> AddAssign<U64Field<M>> for U64Field<M> {
    fn add_assign(&mut self, rhs: U64Field<M>) {
        self.0 = (self.0 + rhs.0) % M;
    }
}

impl<const M: u64> Sub for U64Field<M> {
    type Output = U64Field<M>;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<const M: u64> SubAssign<&U64Field<M>> for U64Field<M> {
    fn sub_assign(&mut self, rhs: &U64Field<M>) {
        *self += -rhs;
    }
}

impl<const M: u64> Neg for U64Field<M> {
    type Output = U64Field<M>;
    fn neg(self) -> Self::Output {
        U64Field((M - self.0) % M)
    }
}

impl<const M: u64> Neg for &U64Field<M> {
    type Output = U64Field<M>;
    fn neg(self) -> Self::Output {
        U64Field((M - self.0) % M)
    }
}

impl<const M: u64> Mul for U64Field<M> {
    type Output = U64Field<M>;
    fn mul(self, rhs: Self) -> Self::Output {
        U64Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<&U64Field<M>> for &U64Field<M> {
    type Output = U64Field<M>;
    fn mul(self, rhs: &U64Field<M>) -> Self::Output {
        U64Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<&U64Field<M>> for U64Field<M> {
    type Output = U64Field<M>;
    fn mul(self, rhs: &U64Field<M>) -> Self::Output {
        U64Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<U64Field<M>> for &U64Field<M> {
    type Output = U64Field<M>;
    fn mul(self, rhs: U64Field<M>) -> Self::Output {
        U64Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<Poly<U64Field<M>>> for &U64Field<M> {
    type Output = Poly<U64Field<M>>;
    fn mul(self, rhs: Poly<U64Field<M>>) -> Self::Output {
        rhs * self
    }
}

impl<const M: u64> Mul<Poly<U64Field<M>>> for U64Field<M> {
    type Output = Poly<U64Field<M>>;
    fn mul(self, rhs: Poly<U64Field<M>>) -> Self::Output {
        rhs * self
    }
}

impl<const M: u64> MulAssign<&U64Field<M>> for U64Field<M> {
    fn mul_assign(&mut self, rhs: &U64Field<M>) {
        self.0 = (self.0 * rhs.0) % M;
    }
}

impl<const M: u64> Div for U64Field<M> {
    type Output = Option<U64Field<M>>;

    fn div(self, rhs: Self) -> Self::Output {
        rhs.inv().map(|v| v * self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f101(n: u64) -> U64Field<101> {
        U64Field::from(n)
    }

    #[test]
    fn test_base() {
        assert_eq!(f101(200), f101(100) + f101(100));
        assert_eq!(f101(100), f101(0) - f101(1));
        assert_eq!(f101(100), f101(0) - f101(1));
        assert_eq!(None, f101(1) / f101(0));
        assert_eq!(f101(4), f101(12) * (f101(4) / f101(12)).unwrap());
    }
    #[test]
    fn test_vectors() {
        assert_eq!(f101(100), -f101(1));
        assert_eq!(f101(50), -(f101(1) / f101(2)).unwrap());
        assert_eq!(f101(20), -(f101(1) / f101(5)).unwrap());
        assert_eq!(f101(1), f101(100).pow(0));
        assert_eq!(f101(100) * f101(100), f101(100).pow(2));
        assert_eq!(f101(100) * f101(100) * f101(100), f101(100).pow(3));
    }
}
