#![allow(clippy::len_without_is_empty)]

use std::fmt::Display;

use super::ec::Field;

// (q_l * a) + (q_r * b) + (q_o * c) + (q_m * a * b) + q_c = 0
// where a,b,c are the left, right and output wires of the gate
pub struct Gate<F: Field> {
    pub q_l: F,
    pub q_r: F,
    pub q_o: F,
    pub q_m: F,
    pub q_c: F,
}

impl<F: Field> Gate<F> {
    pub fn new(q_l: F, q_r: F, q_o: F, q_m: F, q_c: F) -> Self {
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
            q_l: F::one(),
            q_r: F::one(),
            q_o: -F::one(),
            q_m: F::zero(),
            q_c: F::zero(),
        }
    }
    pub fn mul_a_b() -> Self {
        Gate {
            q_l: F::zero(),
            q_r: F::zero(),
            q_o: -F::one(),
            q_m: F::one(),
            q_c: F::zero(),
        }
    }
    pub fn bind_a(value: F) -> Self {
        Gate {
            q_l: F::one(),
            q_r: F::zero(),
            q_o: F::zero(),
            q_m: F::one(),
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

impl<F: Field> Display for Gate<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}a+{}b+{}ab+{}c+{}=0",
            self.q_l, self.q_r, self.q_o, self.q_m, self.q_c
        )
    }
}

#[derive(Debug)]
pub struct Constrains<F: Field> {
    pub q_l: Vec<F>,
    pub q_r: Vec<F>,
    pub q_o: Vec<F>,
    pub q_m: Vec<F>,
    pub q_c: Vec<F>,
    pub c_a: Vec<CopyOf>,
    pub c_b: Vec<CopyOf>,
    pub c_c: Vec<CopyOf>,
}

pub struct Assigment<F: Field> {
    pub a: F,
    pub b: F,
    pub c: F,
}

impl<F: Field> Assigment<F> {
    pub fn new(a: F, b: F, c: F) -> Self {
        Self { a, b, c }
    }
}

pub struct Assigments<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}

impl<F: Field> Constrains<F> {
    pub fn new(
        gates: &[Gate<F>],
        copy_constraints: (Vec<CopyOf>, Vec<CopyOf>, Vec<CopyOf>),
    ) -> Self {
        Self {
            q_l: gates.iter().map(|g| g.q_l).collect(),
            q_r: gates.iter().map(|g| g.q_r).collect(),
            q_o: gates.iter().map(|g| g.q_o).collect(),
            q_m: gates.iter().map(|g| g.q_m).collect(),
            q_c: gates.iter().map(|g| g.q_c).collect(),
            c_a: copy_constraints.0,
            c_b: copy_constraints.1,
            c_c: copy_constraints.2,
        }
    }

    pub fn satisfies(&self, v: &Assigments<F>) -> bool {
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

        // check copy constraints
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

impl<F: Field> Assigments<F> {
    pub fn new(assigments: &[Assigment<F>]) -> Self {
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
