#![allow(clippy::many_single_char_names)]

use crate::{
    field::U64Field,
    pairing::{Pairing, G1, G2, GT},
    poly::Field,
};

use super::{g1::G1P, g2::G2P, gt::GTP, types::f101};

pub struct PBHPairing {}
impl Pairing for PBHPairing {
    type G1 = G1P;
    type G2 = G2P;
    type GT = GTP;

    fn pairing(g1: Self::G1, g2: Self::G2) -> Self::GT {
        let p = <G1P as G1>::F::order();
        let r = G1P::generator_subgroup_size().as_u64();
        let k = G2P::embeeding_degree();

        let exp = (p.pow(k as u32) - 1) / r;

        pairing_f(r, g1, g2).pow(exp)
    }
}

fn pairing_f(r: u64, p: G1P, q: G2P) -> GTP {
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
        GTP::new(f101(1), U64Field::<101>(0))
    } else if r % 2 == 1 {
        let r = r - 1;
        let (x, y, c) = l(p * f101(r), p);
        pairing_f(r, p, q) * GTP::new(q.a * x + c, q.b * y)
    } else {
        let r = r / 2;
        let (x, y, c) = l(p * f101(r), -p * f101(r) * f101(2));
        pairing_f(r, p, q).pow(2) * GTP::new(q.a * x + c, q.b * y)
    }
}

#[cfg(test)]
mod tests {
    use crate::plonk::plonk_by_hand::pairing::PBHPairing;
    use crate::{
        pairing::{Pairing, G1, G2, GT},
        plonk::plonk_by_hand::{g1::G1P, g2::G2P},
    };

    use super::*;

    use std::ops::Mul;

    #[test]
    fn test_parings() {
        let p = G1P::generator().mul(f101(1));
        let r = G1P::generator().mul(f101(4));
        let q = G2P::generator().mul(f101(3));
        let a = f101(5);

        let ê = |g1, g2| PBHPairing::pairing(g1, g2);

        // ê(aP, Q) = ê(P,aQ)

        assert_eq!(ê(p * a, q), ê(p, q * a));

        // ê(aP,Q) = ê(P,Q)^a

        assert_eq!(ê(p * a, q), ê(p, q).pow(a.as_u64()));

        // ê(P+R,Q) = ê(P,Q) * ê(R,Q)

        assert_eq!(ê(p + r, q), ê(p, q) * ê(r, q));
    }
}
