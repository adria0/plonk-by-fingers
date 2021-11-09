#!allow[(unused_imports, dead_code)]

use std::{
    fmt::Display,
    ops::{Add, Index, IndexMut, Mul},
};

use super::field::Field;
use super::poly::Poly;

#[derive(Debug, PartialEq, Clone)]
pub struct Matrix<const M: u64> {
    m: usize, // rows
    n: usize, // cols
    v: Vec<Field<M>>,
}

impl<const M: u64> Matrix<M> {
    pub fn zero(m: usize, n: usize) -> Self {
        let mut v = Vec::new();
        v.resize(n * m, Field::<M>::from(0));
        Self { m, n, v }
    }
    pub fn new(v: &[Field<M>], m: usize, n: usize) -> Self {
        assert_eq!(v.len(), m * n);
        Matrix {
            m,
            n,
            v: v.to_owned(),
        }
    }

    pub fn from(v: &[u64], m: usize, n: usize) -> Self {
        assert_eq!(v.len(), m * n);
        Matrix {
            m,
            n,
            v: v.iter().map(|x| Field::<M>::from(*x)).collect(),
        }
    }
    pub fn cols(&self) -> usize {
        self.n
    }
    pub fn rows(&self) -> usize {
        self.m
    }
    pub fn inv(&self) -> Matrix<M> {
        let len = self.n;
        let mut aug = Matrix::zero(len, len * 2);
        for i in 0..len {
            for j in 0..len {
                aug[(i, j)] = self[(i, j)];
            }
            aug[(i, i + len)] = Field::one();
        }

        aug.gauss_jordan_general();

        let mut unaug = Matrix::zero(len, len);
        for i in 0..len {
            for j in 0..len {
                unaug[(i, j)] = aug[(i, j + len)];
            }
        }
        unaug
    }

    //Begin Generalised Reduced Row Echelon Form
    fn gauss_jordan_general(&mut self) {
        let mut lead = 0;
        let row_count = self.m;
        let col_count = self.n;

        for r in 0..row_count {
            if col_count <= lead {
                break;
            }
            let mut i = r;
            while self[(i, lead)] == Field::zero() {
                i = i + 1;
                if row_count == i {
                    i = r;
                    lead = lead + 1;
                    if col_count == lead {
                        break;
                    }
                }
            }

            // swap i,r rows
            for col in 0..self.n {
                self.v.swap(self.n * i + col, self.n * r + col);
            }

            if self[(r, lead)] != Field::zero() {
                let div = self[(r, lead)];
                for j in 0..col_count {
                    self[(r, j)] = (self[(r, j)] / div).unwrap();
                }
            }

            for k in 0..row_count {
                if k != r {
                    let mult = self[(k, lead)];
                    for j in 0..col_count {
                        self[(k, j)] = self[(k, j)] - self[(r, j)] * mult;
                    }
                }
            }
            lead = lead + 1;
        }
    }
    pub fn into_poly(self) -> Poly<M> {
        assert_eq!(1, self.n);
        Poly::new(self.v)
    }
}

impl<const M: u64> Index<(usize, usize)> for Matrix<M> {
    type Output = Field<M>;
    // row, column
    fn index(&self, p: (usize, usize)) -> &Self::Output {
        assert!(p.0 < self.m && p.1 < self.n);
        &self.v[p.1 + p.0 * self.n]
    }
}

impl<const M: u64> IndexMut<(usize, usize)> for Matrix<M> {
    fn index_mut(&mut self, p: (usize, usize)) -> &mut Self::Output {
        assert!(p.0 < self.m && p.1 < self.n);
        &mut self.v[p.1 + p.0 * self.n]
    }
}

impl<const M: u64> Mul<Matrix<M>> for Matrix<M> {
    type Output = Matrix<M>;
    fn mul(self, rhs: Matrix<M>) -> Self::Output {
        &self * rhs
    }
}

impl<const M: u64> Mul<Matrix<M>> for &Matrix<M> {
    type Output = Matrix<M>;
    fn mul(self, rhs: Matrix<M>) -> Self::Output {
        assert!(self.n == rhs.m);
        let mut c = Matrix::zero(self.m, rhs.n);
        for i in 0..self.m {
            for j in 0..rhs.n {
                for k in 0..self.n {
                    c[(i, j)] = c[(i, j)] + self[(i, k)] * rhs[(k, j)];
                }
            }
        }
        c
    }
}

impl<const M: u64> Add<Matrix<M>> for Matrix<M> {
    type Output = Matrix<M>;
    fn add(self, rhs: Matrix<M>) -> Self::Output {
        assert!(self.m == rhs.m);
        assert!(self.n == rhs.n);
        let mut c = Matrix::zero(self.m, self.n);
        for i in 0..self.v.len() {
            c.v[i] = self.v[i] + rhs.v[i];
        }
        c
    }
}

impl<const M: u64> Display for Matrix<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{}x{}", self.m, self.n)?;
        for r in 0..self.m {
            write!(f, "[")?;
            for c in 0..self.n - 1 {
                write!(f, "{} ", self[(r, c)])?;
            }
            writeln!(f, "{}]", self[(r, self.n - 1)])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = Field<104729>;
    type M = Matrix<104729>;
    #[test]
    fn test_sum_matrix() {
        assert_eq!(
            M::from(&[1, 2], 2, 1) + M::from(&[3, 4], 2, 1),
            M::from(&[4, 6], 2, 1)
        );
    }
    #[test]
    fn test_index_matrix() {
        let m = M::from(&[1, 2, 3, 4], 2, 2);
        // row 0 column 1
        assert_eq!(m[(0, 1)], F::from(2));
    }
    #[test]
    fn test_mul_matrix() {
        assert_eq!(
            M::from(&[1, 2, 3, 4, 5, 6], 2, 3) * M::from(&[10, 11, 20, 21, 30, 31], 3, 2),
            M::from(&[140, 146, 320, 335], 2, 2)
        );
    }
    #[test]
    fn test_inv_matrix() {
        let m = M::from(&[1, 2, 3, 4, 1, 6, 7, 8, 9], 3, 3);
        assert_ne!(m, m.inv());
        assert_eq!(m, m.inv().inv());
    }
}
