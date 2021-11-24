#![allow(clippy::many_single_char_names)]

use super::{
    f101,
    g1::G1P,
    g2::{g2f, G2P},
};
use crate::ec::{Field, G1Point, G2Point, GTPoint, Pairing};

pub struct PBHPairing {}
impl Pairing for PBHPairing {
    type G1 = G1P;
    type G2 = G2P;
    type GT = G2P;

    fn pairing(g1: Self::G1, g2: Self::G2) -> Self::GT {
        let p = <G1P as G1Point>::F::order();
        let r = G1P::generator_subgroup_size().as_u64();
        let k = G2P::embeeding_degree();

        let exp = (p.pow(k as u32) - 1) / r;

        pairing_f(r, g1, g2).pow(exp)
    }
}

fn pairing_f(r: u64, p: G1P, q: G2P) -> G2P {
    // line equation from a to b point
    let l = |a: G1P, b: G1P| {
        let m = b.x - a.x;
        let n = b.y - a.y;

        let x = n;
        let y = -m;
        let c = m * a.y - n * a.x;

        (x, y, c)
    };

    if r == 1 {
        g2f(1, 0)
    } else if r % 2 == 1 {
        let r = r - 1;
        let (x, y, c) = l(p * f101(r), p);
        pairing_f(r, p, q) * G2P::new(q.a * x + c, q.b * y)
    } else {
        let r = r / 2;
        let (x, y, c) = l(p * f101(r), -p * f101(r) * f101(2));
        pairing_f(r, p, q).pow(2) * G2P::new(q.a * x + c, q.b * y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::ops::Mul;

    #[test]
    fn test_parings() {
        let p = G1P::generator().mul(f101(1));
        let q = G2P::generator().mul(f101(3));
        let a = f101(5);

        let ê = |g1, g2| PBHPairing::pairing(g1, g2);

        // ê(aP, Q) = ê(P,aQ)

        assert_eq!(ê(p * a, q), ê(p, q * a));

        // ê(aP,Q) = ê(P,Q)^a

        assert_eq!(ê(p * a, q), ê(p, q).pow(a.as_u64()));

        // ê(aP,Q) = ê((a-1)p,Q) * ê(P,Q)

        assert_eq!(ê(p * a, q), ê(p * (a - f101(1)), q) * ê(p, q));
    }
}
