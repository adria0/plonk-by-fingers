//! This module provides an implementation of polinomials over bls12_381::Field<M>

pub use super::field::Field;
use std::{
    cmp::max,
    fmt::{Display, Formatter, Result},
    iter,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, SubAssign},
};

/// A polinomial with bl12_381::Field<M> coeffs
#[derive(Clone, Debug, PartialEq)]
pub struct Poly<const M: u64>(Vec<Field<M>>);

impl<const M: u64> Poly<M> {
    /// Creates a new Poly from its `coeffs`icients, first element the coefficient for x^0
    /// for safetly, input value is normalized (trailing zeroes are removed)
    pub fn new(coeffs: Vec<Field<M>>) -> Self {
        let mut poly = Poly(coeffs);
        poly.normalize();
        poly
    }

    /// Creates a new polinomial where the `coeffs` fits in u64 values
    pub fn from(coeffs: &[i64]) -> Self {
        Poly::new(
            coeffs
                .iter()
                .map(|n| {
                    if *n >= 0 {
                        Field::<M>::from(*n as u64)
                    } else {
                        -Field::<M>::from(-*n as u64)
                    }
                })
                .collect::<Vec<Field<M>>>(),
        )
    }
    pub fn coeffs(&self) -> &[Field<M>] {
        &self.0
    }
    /// Returns p(x)=0
    pub fn zero() -> Self {
        Poly(vec![Field::<M>::zero()])
    }

    /// Returns p(x)=1
    pub fn one() -> Self {
        Poly(vec![Field::<M>::one()])
    }

    /// Creates a polinomial that contains a set of `p` points, by using lagrange
    /// see <https://en.wikipedia.org/wiki/Lagrange_polynomial>
    pub fn lagrange(p: &[(Field<M>, Field<M>)]) -> Self {
        let k = p.len();
        let mut l = Poly::zero();
        for j in 0..k {
            let mut l_j = Poly::one();
            for i in 0..k {
                if i != j {
                    let c = (p[j].0 - p[i].0).inv();
                    assert!(
                        bool::from(c.is_some()),
                        "lagrange polinomial x points must be unique"
                    );
                    let c = c.unwrap();
                    l_j = &l_j * &Poly::new(vec![-(c * p[i].0), c]);
                }
            }
            l += &(&l_j * &p[j].1);
        }
        l
    }
    /// Creates a polinomial that has roots at the selected points (x-p_1)(x-p_2)...(x-p_n)
    pub fn z(points: &[Field<M>]) -> Poly<M> {
        points.iter().fold(Poly::one(), |acc, x| {
            &acc * &Poly::new(vec![-x, Field::<M>::one()])
        })
    }

    /// Evals the polinomial at the desired point
    pub fn eval(&self, x: &Field<M>) -> Field<M> {
        let mut x_pow = Field::<M>::one();
        let mut y = self.0[0];
        for (i, _) in self.0.iter().enumerate().skip(1) {
            x_pow *= x;
            y += &(x_pow * self.0[i]);
        }
        y
    }

    /// Evals the polinomial suplying the `x_pows` x^0, x^1, x^2
    pub fn eval_with_pows(&self, x_pow: &[Field<M>]) -> Field<M> {
        let mut y = self.0[0];
        for (i, _) in self.0.iter().enumerate() {
            y += &(x_pow[i] * self.0[i]);
        }
        y
    }

    /// Returns the degree of the polinominal, degree(x+1) = 1
    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    /// Normalizes the coefficients, removing ending zeroes
    pub fn normalize(&mut self) {
        if self.0.len() > 1 && self.0[self.0.len() - 1].is_zero() {
            let first_non_zero = self.0.iter().rev().position(|p| *p != Field::<M>::zero());
            if let Some(first_non_zero) = first_non_zero {
                self.0
                    .resize(self.0.len() - first_non_zero, Field::<M>::zero());
            } else {
                self.0.resize(1, Field::<M>::zero());
            }
        }
    }

    /// Returns if p(x)=0
    pub fn is_zero(&self) -> bool {
        self.0.len() == 1 && self.0[0].is_zero()
    }

    /// Sets the `i`-th coefficient to the selected `p` value
    pub fn set(&mut self, i: usize, p: Field<M>) {
        if self.0.len() < i + 1 {
            self.0.resize(i + 1, Field::<M>::zero());
        }
        self.0[i] = p;
        self.normalize();
    }

    /// Returns the `i`-th coefficient
    pub fn get(&mut self, i: usize) -> Option<&Field<M>> {
        self.0.get(i)
    }
}

impl<const M: u64> Display for Poly<M> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let mut first: bool = true;
        for i in 0..=self.degree() {
            if self.0[i].is_zero() && self.degree() > 1 {
                continue;
            }
            let v = format!("{}", self.0[i]);
            if !first {
                write!(f, "+")?;
            }
            if i == 0 || v != "1" {
                write!(f, "{}", v)?;
            }
            if i >= 1 {
                write!(f, "x")?;
            }
            if i >= 2 {
                write!(f, "^{}", i)?;
            }
            first = false;
        }
        Ok(())
    }
}

/* main arithmetics */
impl<const M: u64> AddAssign<&Poly<M>> for Poly<M> {
    fn add_assign(&mut self, rhs: &Poly<M>) {
        for n in 0..max(self.0.len(), rhs.0.len()) {
            if n >= self.0.len() {
                self.0.push(rhs.0[n]);
            } else if n < self.0.len() && n < rhs.0.len() {
                self.0[n] += rhs.0[n];
            }
        }
        self.normalize();
    }
}

impl<const M: u64> AddAssign<&Field<M>> for Poly<M> {
    fn add_assign(&mut self, rhs: &Field<M>) {
        self.0[0] += rhs;
        self.normalize();
    }
}

impl<const M: u64> SubAssign<&Field<M>> for Poly<M> {
    fn sub_assign(&mut self, rhs: &Field<M>) {
        self.0[0] -= rhs;
        self.normalize();
    }
}

impl<const M: u64> SubAssign<&Poly<M>> for Poly<M> {
    fn sub_assign(&mut self, rhs: &Poly<M>) {
        for n in 0..max(self.0.len(), rhs.0.len()) {
            if n >= self.0.len() {
                self.0.push(rhs.0[n]);
            } else if n < self.0.len() && n < rhs.0.len() {
                self.0[n] -= &rhs.0[n];
            }
        }
        self.normalize();
    }
}

impl<const M: u64> Mul<&Poly<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: &Poly<M>) -> Self::Output {
        let mut mul: Vec<Field<M>> = iter::repeat(Field::<M>::zero())
            .take(self.0.len() + rhs.0.len() - 1)
            .collect();
        for n in 0..self.0.len() {
            for m in 0..rhs.0.len() {
                mul[n + m] += self.0[n] * rhs.0[m];
            }
        }
        let mut m = Poly(mul);
        m.normalize();
        m
    }
}

impl<const M: u64> MulAssign<&Field<M>> for Poly<M> {
    fn mul_assign(&mut self, rhs: &Field<M>) {
        if rhs.is_zero() {
            *self = Poly::zero()
        } else {
            self.0.iter_mut().for_each(|v| *v = *v * rhs);
        }
    }
}

impl<const M: u64> Div for Poly<M> {
    type Output = (Poly<M>, Poly<M>);

    fn div(self, rhs: Poly<M>) -> Self::Output {
        let (mut q, mut r) = (Poly::zero(), self);
        while !r.is_zero() && r.degree() >= rhs.degree() {
            let lead_r = r.0[r.0.len() - 1];
            let lead_d = rhs.0[rhs.0.len() - 1];
            let mut t = Poly::zero();
            t.set(r.0.len() - rhs.0.len(), lead_r * lead_d.inv().unwrap());
            q += &t;
            r -= &(&rhs * &t);
        }
        q.normalize();
        r.normalize();
        (q, r)
    }
}

/* helpers */
impl<const M: u64> Add for Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: Poly<M>) -> Self::Output {
        let mut v = self.clone();
        v += &rhs;
        v
    }
}

impl<const M: u64> Add<Poly<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: Poly<M>) -> Self::Output {
        let mut v = self.clone();
        v += &rhs;
        v
    }
}

impl<const M: u64> Add<&Poly<M>> for Poly<M> {
    type Output = Poly<M>;
    fn add(mut self, rhs: &Poly<M>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const M: u64> Add<&Poly<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: &Poly<M>) -> Self::Output {
        let mut v = self.clone();
        v += rhs;
        v
    }
}

impl<const M: u64> Add<&Field<M>> for Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: &Field<M>) -> Self::Output {
        let mut cloned = self.clone();
        cloned += rhs;
        cloned
    }
}

impl<const M: u64> Add<Field<M>> for Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: Field<M>) -> Self::Output {
        let mut cloned = self.clone();
        cloned += &rhs;
        cloned
    }
}

impl<const M: u64> Sub<Field<M>> for Poly<M> {
    type Output = Poly<M>;
    fn sub(self, rhs: Field<M>) -> Self::Output {
        let mut cloned = self.clone();
        cloned -= &rhs;
        cloned
    }
}

impl<const M: u64> Add<Field<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn add(self, rhs: Field<M>) -> Self::Output {
        let mut cloned = self.clone();
        cloned += &rhs;
        cloned
    }
}

impl<const M: u64> Sub<Poly<M>> for Poly<M> {
    type Output = Poly<M>;
    fn sub(mut self, rhs: Poly<M>) -> Self::Output {
        self -= &rhs;
        self
    }
}

impl<const M: u64> Mul<&Poly<M>> for Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: &Poly<M>) -> Self::Output {
        &self * rhs
    }
}

impl<const M: u64> Mul<Poly<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Poly<M>) -> Self::Output {
        self * &rhs
    }
}

impl<const M: u64> Mul<Poly<M>> for Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Poly<M>) -> Self::Output {
        self * &rhs
    }
}

impl<const M: u64> Mul<&Field<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: &Field<M>) -> Self::Output {
        let mut m = self.clone();
        m *= rhs;
        m
    }
}

impl<const M: u64> Mul<&Field<M>> for Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: &Field<M>) -> Self::Output {
        &self * rhs
    }
}

impl<const M: u64> Mul<Field<M>> for Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Field<M>) -> Self::Output {
        self * &rhs
    }
}

impl<const M: u64> Mul<Field<M>> for &Poly<M> {
    type Output = Poly<M>;
    fn mul(self, rhs: Field<M>) -> Self::Output {
        self * &rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = Field<15485863>;
    type P = Poly<15485863>;

    #[test]
    fn test_poly_add() {
        let mut p246 = P::from(&[1, 2, 3]);
        p246 += &P::from(&[1, 2, 3]);
        assert_eq!(p246, P::from(&[2, 4, 6]));

        let mut p24645 = P::from(&[1, 2, 3]);
        p24645 += &P::from(&[1, 2, 3, 4, 5]);
        assert_eq!(p24645, P::from(&[2, 4, 6, 4, 5]));

        let mut p24646 = P::from(&[1, 2, 3, 4, 6]);
        p24646 += &P::from(&[1, 2, 3]);
        assert_eq!(p24646, Poly::from(&[2, 4, 6, 4, 6]));
    }

    #[test]
    fn test_poly_sub() {
        let mut p0 = P::from(&[1, 2, 3]);
        p0 -= &P::from(&[1, 2, 3]);
        assert_eq!(p0, P::from(&[0]));

        let mut p003 = P::from(&[1, 2, 3]);
        p003 -= &P::from(&[1, 2]);
        assert_eq!(p003, P::from(&[0, 0, 3]));
    }

    #[test]
    fn test_poly_mul() {
        assert_eq!(
            &P::from(&[5, 0, 10, 6]) * &P::from(&[1, 2, 4]),
            P::from(&[5, 10, 30, 26, 52, 24])
        );
    }

    #[test]
    fn test_poly_div() {
        fn do_test(n: P, d: P) {
            let (q, r) = n.clone() / d.clone();
            let mut n2 = &q * &d;
            n2 += &r;
            assert_eq!(n, n2);
        }

        do_test(P::from(&[1]), P::from(&[1, 1]));
        do_test(P::from(&[1, 1]), P::from(&[1, 1]));
        do_test(P::from(&[1, 2, 1]), P::from(&[1, 1]));
        do_test(P::from(&[1, 2, 1, 2, 5, 8, 1, 9]), P::from(&[1, 1, 5, 4]));
    }

    #[test]
    fn test_poly_print() {
        assert_eq!("1+2x+x^2", format!("{}", P::from(&[1, 2, 1])));
        assert_eq!("1+x^2", format!("{}", P::from(&[1, 0, 1])));
        assert_eq!("x^2", format!("{}", P::from(&[0, 0, 1])));
        assert_eq!("2x^2", format!("{}", P::from(&[0, 0, 2])));
    }

    #[test]
    fn test_poly_lagrange_multi() {
        let points = vec![
            (F::from(1), F::from(2)),
            (F::from(5), F::from(7)),
            (F::from(7), F::from(9)),
            (F::from(3), F::from(1)),
        ];
        let l = Poly::lagrange(&points);
        points.iter().for_each(|p| assert_eq!(l.eval(&p.0), p.1));
    }
    #[test]
    fn test_poly_z() {
        assert_eq!(
            P::z(&vec![F::from(1), F::from(5)]),
            P::from(&[5, -6, 1]) // f(x) = (x-1)(x-5) = x^2-6x+5
        );
    }
    #[test]
    fn test_poly_eval() {
        // check that (x^2+2x+1)(2) = 9
        assert_eq!(P::from(&[1, 2, 1]).eval(&F::from(2)), F::from(9));
    }
    #[test]
    fn test_poly_normalize() {
        let mut p1 = P::from(&[1, 0, 0, 0]);
        p1.normalize();
        assert_eq!(p1, P::from(&[1]));
    }
}
