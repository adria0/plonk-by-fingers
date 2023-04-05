pub use crate::ec::Field;
pub use crate::matrix::Matrix;
pub use anyhow::{anyhow, bail, Result};
use std::convert::TryFrom;
use std::{
    cmp::max,
    fmt::{Display, Formatter},
    ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, SubAssign},
};

#[derive(Clone, Debug, PartialEq)]
pub struct Poly<F: Field>(Vec<F>);

impl<F: Field> Poly<F> {
    /// Creates a new Poly from its `coeffs`icients, first element the coefficient for x^0
    /// for safetly, input value is normalized (trailing zeroes are removed)
    pub fn new<I: IntoIterator<Item=F>>(coeffs: I) -> Self {
        let mut poly = Poly(coeffs.into_iter().collect());
        poly.normalize();
        poly
    }


    /// Creates a new polinomial where the `coeffs` fits in u64 values
    pub fn from<I: IntoIterator<Item=i64>>(coeffs: I) -> Self {
        Poly::new(coeffs.into_iter().map(|n| F::from(n)).collect::<Vec<F>>())
    }
    // Returns the x polinomial
    pub fn x() -> Self {
        Poly::new(vec![F::zero(), F::one()])
    }

    // Returns a constant polinomial
    pub fn n(n: F) -> Self {
        Poly::new(vec![n])
    }

    /// Parses an expression
    pub fn parse(mut s: &str) -> Result<Poly<F>> {
        let orig = s;
        let mut coefficients = Vec::new();
        while !s.is_empty() {
            let sign: i64 = if let Some(next) = s.strip_prefix("+") {
                s = next;
                1
            } else if let Some(next) = s.strip_prefix("-") {
                s = next;
                -1
            } else {
                1
            };
            let expr = if let Some(next_pos) = s.find(&['-', '+']) {
                let expr = &s[..next_pos];
                s = &s[next_pos..];
                expr
            } else {
                let expr = s;
                s = &s[s.len()..];
                expr
            };
            let (mut coeff, exp) = if let Some(xexpr) = expr.find("x^") {
                (&expr[..xexpr], &expr[xexpr + 2..])
            } else if let Some(xexpr) = expr.find("x") {
                (&expr[..xexpr], "1")
            } else {
                (expr, "0")
            };
            if coeff.len() == 0 {
                coeff = "1";
            }
            let coeff = sign * i64::from_str_radix(coeff, 10).unwrap();
            let exp = usize::from_str_radix(exp, 10).unwrap();
            while coefficients.len() <= exp {
                coefficients.push(0i64);
            }
            coefficients[exp] += coeff;
        }
        let poly = Self::from(coefficients);
        if poly.to_string() != orig {
            bail!(
                "Not canonical poly orig:'{}' parsed:'{}'",
                orig,
                poly.to_string()
            );
        }
        Ok(poly)
    }

    pub fn coeffs(&self) -> &[F] {
        &self.0
    }
    pub fn into_coeffs(self) -> Vec<F> {
        self.0
    }
    /// Returns p(x)=0
    pub fn zero() -> Self {
        Poly(vec![F::zero()])
    }

    /// Returns p(x)=1
    pub fn one() -> Self {
        Poly(vec![F::one()])
    }

    /// Creates a polinomial that contains a set of `p` points, by using lagrange
    /// see <https://en.wikipedia.org/wiki/Lagrange_polynomial>
    pub fn lagrange(p: &[(F, F)]) -> Self {
        let k = p.len();
        let mut l = Poly::zero();
        for j in 0..k {
            let mut l_j = Poly::one();
            for i in 0..k {
                if i != j {
                    let c = (p[j].0 - p[i].0).inv();
                    assert!(c.is_some(), "lagrange polinomial x points must be unique");
                    let c = c.unwrap();
                    l_j = &l_j * &Poly::new(vec![-(c * p[i].0), c]);
                }
            }
            l += &(&l_j * p[j].1);
        }
        l
    }

    /// Creates a polinomial that has roots at the selected points (x-p_1)(x-p_2)...(x-p_n)
    pub fn z<I: IntoIterator<Item=F>>(points: I) -> Self {
        points
            .into_iter()
            .fold(Poly::one(), |acc, x| &acc * &Poly::new(vec![-x, F::one()]))
    }

    /// Evals the polinomial at the desired point
    pub fn eval(&self, x: &F) -> F {
        let mut x_pow = F::one();
        let mut y = self.0[0];
        for (i, _) in self.0.iter().enumerate().skip(1) {
            x_pow *= x;
            y += x_pow * self.0[i];
        }
        y
    }

    /// Substitutes x with p(x)
    pub fn subst_x(&self, p: &Poly<F>) -> Poly<F> {
        let mut res = Poly::new(vec![self.0[0]]);
        let mut p_pow = p.clone();
        for coeff in self.0.iter().skip(1) {
            res = res + &p_pow * coeff;
            p_pow = p_pow * p;
        }
        res
    }

    /// Pow
    pub fn pow(&self, mut exp: u64) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base.clone();
            }
            exp >>= 1;
            base = base.clone() * base;
        }
        result
    }

    /// Evals the polinomial suplying the `x_pows` x^0, x^1, x^2
    pub fn eval_with_pows(&self, x_pow: &[F]) -> F {
        let mut y = self.0[0];
        for (i, _) in self.0.iter().enumerate() {
            y += x_pow[i] * self.0[i];
        }
        y
    }

    /// Returns the degree of the polinominal, degree(x+1) = 1
    pub fn degree(&self) -> u64 {
        (self.0.len() - 1) as u64
    }

    /// Normalizes the coefficients, removing ending zeroes
    pub fn normalize(&mut self) {
        if self.0.len() > 1 && self.0[self.0.len() - 1].is_zero() {
            let first_non_zero = self.0.iter().rev().position(|p| !p.is_zero());
            if let Some(first_non_zero) = first_non_zero {
                self.0.resize(self.0.len() - first_non_zero, F::zero());
            } else {
                self.0.resize(1, F::zero());
            }
        }
    }

    /// Returns if p(x)=0
    pub fn is_zero(&self) -> bool {
        self.0.len() == 1 && self.0[0].is_zero()
    }

    /// Sets the `i`-th coefficient to the selected `p` value
    pub fn set(&mut self, i: usize, p: F) {
        if self.0.len() < i + 1 {
            self.0.resize(i + 1, F::zero());
        }
        self.0[i] = p;
        self.normalize();
    }

    /// Returns the `i`-th coefficient
    pub fn get(&mut self, i: usize) -> Option<&F> {
        self.0.get(i)
    }
}

impl<F: Field> Display for Poly<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut first: bool = true;
        for i in 0..=self.degree() as usize {
            if self.0[i].is_zero() && self.degree() > 0 {
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

impl<F: Field> TryFrom<Matrix<F>> for Poly<F> {
    type Error = anyhow::Error;
    fn try_from(matrix: Matrix<F>) -> Result<Self, Self::Error> {
        if matrix.cols() == 1 {
            Ok(Poly::new(Into::<Vec<_>>::into(matrix)))
        } else {
            Err(anyhow!("only one row"))
        }
    }
}

/* main arithmetics */
impl<F: Field> AddAssign<&Poly<F>> for Poly<F> {
    fn add_assign(&mut self, rhs: &Poly<F>) {
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

impl<F: Field> AddAssign<&F> for Poly<F> {
    fn add_assign(&mut self, rhs: &F) {
        self.0[0] += rhs;
        self.normalize();
    }
}

impl<F: Field> SubAssign<&F> for Poly<F> {
    fn sub_assign(&mut self, rhs: &F) {
        self.0[0] -= rhs;
        self.normalize();
    }
}

impl<F: Field> SubAssign<&Poly<F>> for Poly<F> {
    fn sub_assign(&mut self, rhs: &Poly<F>) {
        for n in 0..max(self.0.len(), rhs.0.len()) {
            if n >= self.0.len() {
                self.0.push(-rhs.0[n]);
            } else if n < self.0.len() && n < rhs.0.len() {
                self.0[n] -= &rhs.0[n];
            }
        }
        self.normalize();
    }
}

impl<F: Field> Mul<&Poly<F>> for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: &Poly<F>) -> Self::Output {
        let mut mul = vec![F::zero(); self.0.len() + rhs.0.len()];
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

impl<F: Field> MulAssign<&F> for Poly<F> {
    fn mul_assign(&mut self, rhs: &F) {
        if rhs.is_zero() {
            *self = Poly::zero()
        } else {
            self.0.iter_mut().for_each(|v| *v = *v * rhs);
        }
    }
}

impl<F: Field> Div for Poly<F> {
    type Output = (Poly<F>, Poly<F>);

    fn div(self, rhs: Poly<F>) -> Self::Output {
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
impl<F: Field> Add for Poly<F> {
    type Output = Poly<F>;
    fn add(mut self, rhs: Poly<F>) -> Self::Output {
        self += &rhs;
        self
    }
}

impl<F: Field> Add<Poly<F>> for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, rhs: Poly<F>) -> Self::Output {
        let mut v = self.clone();
        v += &rhs;
        v
    }
}

impl<F: Field> Add<&Poly<F>> for Poly<F> {
    type Output = Poly<F>;
    fn add(mut self, rhs: &Poly<F>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<F: Field> Add<&Poly<F>> for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, rhs: &Poly<F>) -> Self::Output {
        let mut v = self.clone();
        v += rhs;
        v
    }
}

impl<F: Field> Add<&F> for Poly<F> {
    type Output = Poly<F>;
    fn add(mut self, rhs: &F) -> Self::Output {
        self += rhs;
        self
    }
}

impl<F: Field> Add<F> for Poly<F> {
    type Output = Poly<F>;
    fn add(mut self, rhs: F) -> Self::Output {
        self += &rhs;
        self
    }
}

impl<F: Field> Sub<F> for Poly<F> {
    type Output = Poly<F>;
    fn sub(mut self, rhs: F) -> Self::Output {
        self -= &rhs;
        self
    }
}

impl<F: Field> Add<F> for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, rhs: F) -> Self::Output {
        let mut cloned = self.clone();
        cloned += &rhs;
        cloned
    }
}

impl<F: Field> Sub<Poly<F>> for Poly<F> {
    type Output = Poly<F>;
    fn sub(mut self, rhs: Poly<F>) -> Self::Output {
        self -= &rhs;
        self
    }
}

impl<F: Field> Mul<&Poly<F>> for Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: &Poly<F>) -> Self::Output {
        &self * rhs
    }
}

impl<F: Field> Mul<Poly<F>> for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: Poly<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F: Field> Mul<Poly<F>> for Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: Poly<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F: Field> Mul<&F> for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: &F) -> Self::Output {
        let mut m = self.clone();
        m *= rhs;
        m
    }
}

impl<F: Field> Mul<&F> for Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: &F) -> Self::Output {
        &self * rhs
    }
}

#[allow(clippy::op_ref)]
impl<F: Field> Mul<F> for Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: F) -> Self::Output {
        self * &rhs
    }
}

#[allow(clippy::op_ref)]
impl<F: Field> Mul<F> for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, rhs: F) -> Self::Output {
        self * &rhs
    }
}

struct PolyFreq<F: Field>(Vec<F>);

impl<F: Field> Mul for PolyFreq<F> {
    type Output = PolyFreq<F>;
    fn mul(self, rhs: PolyFreq<F>) -> Self::Output {
        let zero = F::zero();
        let mut c_freq = Vec::new();

        for n in 0..self.0.len() {
            let l = self.0.get(n).unwrap_or(&zero);
            let r = rhs.0.get(n).unwrap_or(&zero);
            c_freq.push(*l * r);
        }
        Self(c_freq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::U64Field;
    type F = U64Field<15485863>;
    type P = Poly<F>;

    #[test]
    fn test_poly_add() {
        let mut p246 = P::from([1, 2, 3]);
        p246 += &P::from([1, 2, 3]);
        assert_eq!(p246, P::from([2, 4, 6]));

        let mut p24645 = P::from([1, 2, 3]);
        p24645 += &P::from([1, 2, 3, 4, 5]);
        assert_eq!(p24645, P::from([2, 4, 6, 4, 5]));

        let mut p24646 = P::from([1, 2, 3, 4, 6]);
        p24646 += &P::from([1, 2, 3]);
        assert_eq!(p24646, P::from([2, 4, 6, 4, 6]));
    }

    #[test]
    fn test_poly_sub() {
        let mut p0 = P::from([1, 2, 3]);
        p0 -= &P::from([1, 2, 3]);
        assert_eq!(p0, P::from([0]));

        let mut p003 = P::from([1, 2, 3]);
        p003 -= &P::from([1, 2]);
        assert_eq!(p003, P::from([0, 0, 3]));
    }

    #[test]
    fn test_poly_mul() {
        assert_eq!(
            &P::from([5, 0, 10, 6]) * &P::from([1, 2, 4]),
            P::from([5, 10, 30, 26, 52, 24])
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

        do_test(P::from([1]), P::from([1, 1]));
        do_test(P::from([1, 1]), P::from([1, 1]));
        do_test(P::from([1, 2, 1]), P::from([1, 1]));
        do_test(P::from([1, 2, 1, 2, 5, 8, 1, 9]), P::from([1, 1, 5, 4]));
    }

    #[test]
    fn test_poly_print() {
        assert_eq!("1+2x+x^2", format!("{}", P::from([1, 2, 1])));
        assert_eq!("1+x^2", format!("{}", P::from([1, 0, 1])));
        assert_eq!("x^2", format!("{}", P::from([0, 0, 1])));
        assert_eq!("2x^2", format!("{}", P::from([0, 0, 2])));
    }

    #[test]
    fn test_poly_lagrange_multi() {
        let points = vec![
            (F::from(1u64), F::from(2u64)),
            (F::from(5u64), F::from(7u64)),
            (F::from(7u64), F::from(9u64)),
            (F::from(3u64), F::from(1u64)),
        ];
        let l = Poly::lagrange(&points);
        points.iter().for_each(|p| assert_eq!(l.eval(&p.0), p.1));
    }
    #[test]
    fn test_poly_z() {
        assert_eq!(
            P::z([F::from(1u64), F::from(5u64)]),
            P::from([5, -6, 1]) // f(x) = (x-1)(x-5) = x^2-6x+5
        );
    }
    #[test]
    fn test_poly_eval() {
        // check that (x^2+2x+1)(2) = 9
        assert_eq!(P::from([1, 2, 1]).eval(&F::from(2u64)), F::from(9u64));
    }
    #[test]
    fn test_poly_normalize() {
        let mut p1 = P::from([1, 0, 0, 0]);
        p1.normalize();
        assert_eq!(p1, P::from([1]));
    }
    #[test]
    fn test_parse() {
        Poly::<F>::parse("1").unwrap();
        Poly::<F>::parse("x").unwrap();
        Poly::<F>::parse("2x").unwrap();
        Poly::<F>::parse("x^2").unwrap();
        Poly::<F>::parse("1+x+2x^12").unwrap();
    }

    #[test]
    fn test_pow() {
        let check = |base: &str, exp: u64, res: &str| {
            let base = Poly::<F>::parse(&base).unwrap();
            let res = Poly::<F>::parse(&res).unwrap();
            assert_eq!(base.pow(exp), res);
        };
        check("2", 3, "8");
        check("x", 0, "1");
        check("x", 6, "x^6");
        check("1+x", 2, "1+2x+x^2");
    }

    #[test]
    fn test_subst() {
        let check = |f_str: &str, subst_str: &str| {
            let val = F::from(13331u64);
            let f = Poly::<F>::parse(&f_str).unwrap();
            let subst = Poly::<F>::parse(&subst_str).unwrap();
            let ev1 = f.eval(&subst.eval(&val));
            let ev2 = f.subst_x(&subst).eval(&val);
            assert_eq!(ev1, ev2, "{} {} => {}", f, subst, f.subst_x(&subst));
        };

        check("2", "2x");
        check("2x", "2x");
        check("1+2x^3", "1+x+2x^2");
        check("x^2", "6+16x+2x^2+13x^3");
    }
}
