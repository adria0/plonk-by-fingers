use std::{
    fmt::Display,
    ops::{Add, Mul, Neg},
};

use super::{f101, F101};
use crate::ec::{Field, G2Point};

#[allow(non_snake_case)]
pub fn g2f(a: u64, b: u64) -> G2P {
    G2P::new(F101::from(a), F101::from(b))
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct G2P {
    pub a: F101,
    pub b: F101,
}

impl G2Point for G2P {
    type F = F101;
    type S = F101;

    // a + b · u
    fn new(a: Self::F, b: Self::F) -> Self {
        G2P { a, b }
    }
    fn generator() -> Self {
        G2P {
            a: f101(36u64),
            b: f101(31u64),
        }
    }
    fn embeeding_degree() -> u64 {
        2
    }
    fn x(&self) -> &Self::F {
        &self.a
    }
    fn y(&self) -> &Self::F {
        &self.b
    }
}

impl Display for G2P {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}+{}u", self.a, self.b)
    }
}

impl Neg for G2P {
    type Output = G2P;
    fn neg(self) -> Self::Output {
        G2P::new(self.a, -self.b)
    }
}

impl Add for G2P {
    type Output = G2P;
    fn add(self, rhs: G2P) -> Self {
        if self == rhs {
            let two = f101(2);
            let three = f101(3);
            let m_u = ((three * self.a.pow(2)) / (two * self.b)).unwrap(); // in u units
            let u_pow_2_inv = (-f101(2)).inv().unwrap(); // 1/(u^2) = -1/2
            let m_pow_2 = m_u.pow(2) * u_pow_2_inv;
            G2P::new(
                m_pow_2 - two * self.a,
                u_pow_2_inv * m_u * (three * self.a - m_pow_2) - self.b,
            )
        } else {
            let lambda_u = ((rhs.b - self.b) / (rhs.a - self.a)).unwrap();
            let lambda_pow_2 = lambda_u.pow(2) * -f101(2);
            let a = lambda_pow_2 - self.a - rhs.a;
            let b = lambda_u * (self.a - a) - self.b;

            G2P::new(a, b)
        }
    }
}

impl Mul<F101> for G2P {
    type Output = G2P;
    fn mul(self, rhs: F101) -> Self::Output {
        let mut rhs = rhs.as_u64();
        let mut result = None;
        let mut base = self;
        while rhs > 0 {
            if rhs % 2 == 1 {
                result = Some(if let Some(result) = result {
                    result + base
                } else {
                    base
                })
            }
            rhs >>= 1;
            base = base + base;
        }
        result.unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g2_vectors() {
        let g = G2P::generator();

        // check point doubling
        assert_eq!(g2f(90, 82), g + g);

        // check point addition
        assert_eq!((g + g) + (g + g), g + g + g + g);

        // check point multiplication
        assert_eq!(g * f101(6), g + g + g + g + g + g);

        // check G2 multiplication
        assert_eq!(g2f(26, 97) * g2f(93, 76), g2f(97, 89));
    }
}
