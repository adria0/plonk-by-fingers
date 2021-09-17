#![allow(dead_code)]

use super::poly::Poly;
use std::{
    fmt::Display,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    let (mut s, mut old_s) = (0, 1);
    let (mut t, mut old_t) = (1, 0);
    let (mut r, mut old_r) = (b, a);

    while r != 0 {
        let quotient = &old_r / &r;
        old_r -= &quotient * &r;
        std::mem::swap(&mut old_r, &mut r);
        old_s -= &quotient * &s;
        std::mem::swap(&mut old_s, &mut s);
        old_t -= quotient * &t;
        std::mem::swap(&mut old_t, &mut t);
    }
    (old_r, old_s, old_t)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Field<const M: u64>(u64);

#[allow(non_snake_case)]
pub fn Fi64<const M: u64>(n: i64) -> Field<M> {
    if n < 0 {
        -Field::from(-n as u64)
    } else {
        Field::from(n as u64)
    }
}

impl<const M: u64> Field<M> {
    pub fn from(n: u64) -> Self {
        Self(n % M)
    }
    pub fn zero() -> Self {
        Self::from(0)
    }
    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }
    pub fn one() -> Self {
        Self::from(1)
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn in_field(&self) -> bool {
        self.0 < M
    }
    pub fn rebase<const M2: u64>(&self) -> Field<M2> {
        Field::<M2>::from(self.0)
    }
    pub fn inv(&self) -> Option<Self> {
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
    pub fn pow(&self, mut exp: u64) -> Self {
        let mut result = Self::one();
        let mut base = *self;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            exp = exp >> 1;
            base = base * base;
        }
        result
    }
    pub fn as_poly(&self) -> Poly<M> {
        Poly::new(vec![*self])
    }
}

impl<const M: u64> Display for Field<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const M: u64> Add for Field<M> {
    type Output = Field<M>;
    fn add(self, rhs: Self) -> Self::Output {
        Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<&Field<M>> for Field<M> {
    type Output = Field<M>;
    fn add(self, rhs: &Field<M>) -> Self::Output {
        Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<Field<M>> for &Field<M> {
    type Output = Field<M>;
    fn add(self, rhs: Field<M>) -> Self::Output {
        Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> Add<&Field<M>> for &Field<M> {
    type Output = Field<M>;
    fn add(self, rhs: &Field<M>) -> Self::Output {
        Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u64> AddAssign<&Field<M>> for Field<M> {
    fn add_assign(&mut self, rhs: &Field<M>) {
        self.0 = (self.0 + rhs.0) % M;
    }
}

impl<const M: u64> AddAssign<Field<M>> for Field<M> {
    fn add_assign(&mut self, rhs: Field<M>) {
        self.0 = (self.0 + rhs.0) % M;
    }
}

impl<const M: u64> Sub for Field<M> {
    type Output = Field<M>;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<const M: u64> SubAssign<&Field<M>> for Field<M> {
    fn sub_assign(&mut self, rhs: &Field<M>) {
        *self += -rhs;
    }
}

impl<const M: u64> Neg for Field<M> {
    type Output = Field<M>;
    fn neg(self) -> Self::Output {
        Field((M - self.0) % M)
    }
}

impl<const M: u64> Neg for &Field<M> {
    type Output = Field<M>;
    fn neg(self) -> Self::Output {
        Field((M - self.0) % M)
    }
}

impl<const M: u64> Mul for Field<M> {
    type Output = Field<M>;
    fn mul(self, rhs: Self) -> Self::Output {
        Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<&Field<M>> for &Field<M> {
    type Output = Field<M>;
    fn mul(self, rhs: &Field<M>) -> Self::Output {
        Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<&Field<M>> for Field<M> {
    type Output = Field<M>;
    fn mul(self, rhs: &Field<M>) -> Self::Output {
        Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<Field<M>> for &Field<M> {
    type Output = Field<M>;
    fn mul(self, rhs: Field<M>) -> Self::Output {
        Field((self.0 * rhs.0) % M)
    }
}

impl<const M: u64> Mul<Poly<M>> for &Field<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Poly<M>) -> Self::Output {
        rhs * self
    }
}

impl<const M: u64> Mul<Poly<M>> for Field<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Poly<M>) -> Self::Output {
        rhs * self
    }
}

impl<const M: u64> MulAssign<&Field<M>> for Field<M> {
    fn mul_assign(&mut self, rhs: &Field<M>) {
        self.0 = (self.0 * rhs.0) % M;
    }
}

impl<const M: u64> Div for Field<M> {
    type Output = Option<Field<M>>;

    fn div(self, rhs: Self) -> Self::Output {
        rhs.inv().map(|v| v * self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f101(n: u64) -> Field<101> {
        Field::from(n)
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
