#![allow(dead_code, unused_imports)]

use std::fmt::Display;

use super::ec::{g1f, g2f, G1Point, G2Point};
use super::field::Field;
use super::matrix::Matrix;
use super::plonk::{f101, F17, P17};
use super::poly::Poly;

// (q_l * a) + (q_r * b) + (q_o * c) + (q_m * a * b) + q_c = 0
// where a,b,c are the left, right and output wires of the gate
pub struct Gate {
    pub q_l: F17,
    pub q_r: F17,
    pub q_o: F17,
    pub q_m: F17,
    pub q_c: F17,
}

impl Gate {
    pub fn new(q_l: F17, q_r: F17, q_o: F17, q_m: F17, q_c: F17) -> Self {
        Gate {
            q_l,
            q_r,
            q_o,
            q_m,
            q_c,
        }
    }
    pub fn sum_a_b() -> Self {
        Gate {
            q_l: F17::one(),
            q_r: F17::one(),
            q_o: -F17::one(),
            q_m: F17::zero(),
            q_c: F17::zero(),
        }
    }
    pub fn mul_a_b() -> Self {
        Gate {
            q_l: F17::zero(),
            q_r: F17::zero(),
            q_o: -F17::one(),
            q_m: F17::one(),
            q_c: F17::zero(),
        }
    }
    pub fn bind_a(value: F17) -> Self {
        Gate {
            q_l: F17::one(),
            q_r: F17::zero(),
            q_o: F17::zero(),
            q_m: F17::one(),
            q_c: value,
        }
    }
}

#[derive(Debug)]
pub enum CopyOf {
    A(usize),
    B(usize),
    C(usize),
}

impl Display for Gate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}a+{}b+{}ab+{}c+{}=0",
            self.q_l, self.q_r, self.q_o, self.q_m, self.q_c
        )
    }
}

#[derive(Debug)]
pub struct Constrains {
    pub q_l: Vec<F17>,
    pub q_r: Vec<F17>,
    pub q_o: Vec<F17>,
    pub q_m: Vec<F17>,
    pub q_c: Vec<F17>,
    pub c_a: Vec<CopyOf>,
    pub c_b: Vec<CopyOf>,
    pub c_c: Vec<CopyOf>,
}

pub struct Assigment {
    pub a: F17,
    pub b: F17,
    pub c: F17,
}

impl Assigment {
    pub fn new(a: F17, b: F17, c: F17) -> Self {
        Self { a, b, c }
    }
}

pub struct Assigments {
    pub a: Vec<F17>,
    pub b: Vec<F17>,
    pub c: Vec<F17>,
}

impl Constrains {
    pub fn new(gates: &[Gate], copy_constrains: (Vec<CopyOf>, Vec<CopyOf>, Vec<CopyOf>)) -> Self {
        Self {
            q_l: gates.iter().map(|g| g.q_l).collect(),
            q_r: gates.iter().map(|g| g.q_r).collect(),
            q_o: gates.iter().map(|g| g.q_o).collect(),
            q_m: gates.iter().map(|g| g.q_m).collect(),
            q_c: gates.iter().map(|g| g.q_c).collect(),
            c_a: copy_constrains.0,
            c_b: copy_constrains.1,
            c_c: copy_constrains.2,
        }
    }

    pub fn satisfies(&self, v: &Assigments) -> bool {
        // check gates (q_l * a) + (q_r * b) + (q_o * c) + (q_m * a * b) + q_c = 0
        assert_eq!(v.a.len(), self.q_l.len());
        for n in 0..v.a.len() {
            let r = self.q_l[n] * v.a[n]
                + self.q_l[n] * v.b[n]
                + self.q_o[n] * v.c[n]
                + self.q_m[n] * v.a[n] * v.b[n]
                + self.q_c[n];
            if r != Field::zero() {
                return false;
            }
        }

        // check copy constrains
        assert_eq!(v.a.len(), self.c_a.len());
        assert_eq!(v.a.len(), self.c_b.len());
        assert_eq!(v.a.len(), self.c_c.len());
        for n in 0..self.c_a.len() {
            let value = |c: &CopyOf| match c {
                CopyOf::A(n) => &v.a[*n - 1],
                CopyOf::B(n) => &v.b[*n - 1],
                CopyOf::C(n) => &v.c[*n - 1],
            };
            if &v.a[n] != value(&self.c_a[n])
                || &v.b[n] != value(&self.c_b[n])
                || &v.c[n] != value(&self.c_c[n])
            {
                return false;
            }
        }
        true
    }
}

impl Assigments {
    pub fn new(assigments: &[Assigment]) -> Self {
        Self {
            a: assigments.iter().map(|v| v.a).collect(),
            b: assigments.iter().map(|v| v.b).collect(),
            c: assigments.iter().map(|v| v.c).collect(),
        }
    }
    pub fn len(&self) -> usize {
        self.a.len()
    }
}
