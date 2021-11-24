#![allow(dead_code)]

use std::ops::Mul;

use super::{f101, g2::G2P, F101};
use crate::ec::{Field, G2Point, GTPoint};

impl GTPoint for G2P {
    type S = u64;
    fn pow(&self, mut n: Self::S) -> Self {
        // frobenious map reduction: p^101 = -p
        let (mut p, mut base) = if n >= 101 {
            let base = -self.pow(n / 101);
            n %= 101;
            (base, *self)
        } else {
            (G2P::new(F101::one(), F101::zero()), *self)
        };

        // montgomery reduction
        while n > 0 {
            if n % 2 == 1 {
                p = p * base;
            }
            n >>= 1;
            base = base * base;
        }

        p
    }
}
impl Mul<G2P> for G2P {
    type Output = G2P;
    fn mul(self, rhs: G2P) -> Self::Output {
        G2P::new(
            self.a * rhs.a - f101(2) * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{super::g2::g2f, *};
    #[test]
    fn test_gt_vectors() {
        // check GT multiplication
        assert_eq!(g2f(26, 97) * g2f(93, 76), g2f(97, 89));

        // check GT exp
        assert_eq!(g2f(42, 49).pow(6), g2f(97, 89));
        assert_eq!(g2f(93, 76).pow(101), -g2f(93, 76));
        assert_eq!(g2f(93, 76).pow(102), (-g2f(93, 76)) * g2f(93, 76));
        assert_eq!(g2f(68, 47).pow(600), g2f(97, 89));
    }
}
