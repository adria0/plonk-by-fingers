use crate::poly::Poly;

use super::PlonkTypes;
use crate::pairing::{G1, G2};

pub struct Srs<P: PlonkTypes> {
    pub g1s: Vec<P::G1>,
    pub g2_1: P::G2,
    pub g2_s: P::G2,
}

impl<P: PlonkTypes> Srs<P> {
    pub fn create(s: P::GF, n: usize) -> Self {
        let mut g1s = Vec::new();
        let mut s_pow = s;
        g1s.push(P::G1::generator());
        for _ in 0..n {
            g1s.push(P::G1::generator() * s_pow);
            s_pow = s_pow * s;
        }
        Self {
            g1s,
            g2_1: P::G2::generator(),
            g2_s: P::G2::generator() * s,
        }
    }

    // evaluate a polinomil at secret point s using SRS G1s
    pub fn eval_at_s(&self, vs: &Poly<P::HF>) -> P::G1 {
        vs.coeffs()
            .iter()
            .enumerate()
            .fold(G1::identity(), |acc, (n, v)| acc + self.g1s[n] * P::gf(*v))
    }
}
