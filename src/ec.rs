#![allow(dead_code)]

use std::{
    fmt::Display,
    ops::{Add, Mul, Neg},
};

use super::plonk::{f101, F101};

#[allow(non_snake_case)]
pub fn g1f(x: u64, y: u64) -> G1Point {
    G1Point::new(f101(x), f101(y))
}

#[allow(non_snake_case)]
pub fn g2f(a: u64, b: u64) -> G2Point {
    G2Point::new(f101(a), f101(b))
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct G1Point {
    pub x: F101,
    pub y: F101,
    pub inf: bool,
}
impl G1Point {
    pub fn new(x: F101, y: F101) -> Self {
        G1Point { x, y, inf: false }
    }
    pub fn in_curve(&self) -> bool {
        self.y.pow(2) == self.x.pow(3) + f101(3)
    }
    pub fn is_inf(&self) -> bool {
        self.inf
    }
    pub fn inf() -> Self {
        G1Point {
            x: f101(0),
            y: f101(0),
            inf: true,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct G2Point {
    pub a: F101,
    pub b: F101,
}
impl G2Point {
    // a + b · u
    pub fn new(a: F101, b: F101) -> Self {
        G2Point { a, b }
    }

    pub fn pow(&self, mut n: u64) -> Self {
        // frobenious map reduction: p^101 = -p
        let (mut p, mut base) = if n >= 101 {
            let base = -self.pow(n / 101);
            n = n % 101;
            (base, *self)
        } else {
            (G2Point::new(f101(1), f101(0)), *self)
        };

        // montgomery reduction
        while n > 0 {
            if n % 2 == 1 {
                p = p * base;
            }
            n = n >> 1;
            base = base * base;
        }

        p
    }
}

impl Display for G1Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.inf {
            write!(f, "inf")
        } else {
            write!(f, "({},{})", self.x, self.y)
        }
    }
}

impl Display for G2Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}+{}u", self.a, self.b)
    }
}

impl Neg for G1Point {
    type Output = G1Point;
    fn neg(self) -> Self::Output {
        if self.inf {
            self
        } else {
            G1Point::new(self.x, -self.y)
        }
    }
}

impl Neg for G2Point {
    type Output = G2Point;
    fn neg(self) -> Self::Output {
        G2Point::new(self.a, -self.b)
    }
}

impl Add for G1Point {
    type Output = G1Point;
    fn add(self, rhs: G1Point) -> Self {
        if self.inf {
            rhs
        } else if rhs.inf {
            self
        } else if self == -rhs {
            G1Point::inf()
        } else if self == rhs {
            let two = F101::from(2);
            let three = F101::from(3);
            let m = ((three * self.x.pow(2)) / (two * self.y)).unwrap();
            G1Point::new(
                m * m - two * self.x,
                m * (three * self.x - m.pow(2)) - self.y,
            )
        } else {
            // https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#G1Point_addition
            let lambda = ((rhs.y - self.y) / (rhs.x - self.x))
                .expect(&format!("cannot add {}+{}", self, rhs));
            let x = lambda.pow(2) - self.x - rhs.x;
            G1Point::new(x, lambda * (self.x - x) - self.y)
        }
    }
}

impl Add for G2Point {
    type Output = G2Point;
    fn add(self, rhs: G2Point) -> Self {
        if self == rhs {
            let two = F101::from(2);
            let three = F101::from(3);
            let m_u = ((three * self.a.pow(2)) / (two * self.b)).unwrap(); // in u units
            let u_pow_2_inv = (-F101::from(2)).inv().unwrap(); // 1/(u^2) = -1/2
            let m_pow_2 = m_u.pow(2) * u_pow_2_inv;
            G2Point::new(
                m_pow_2 - two * self.a,
                u_pow_2_inv * m_u * (three * self.a - m_pow_2) - self.b,
            )
        } else {
            let lambda_u = ((rhs.b - self.b) / (rhs.a - self.a)).unwrap();
            let lambda_pow_2 = lambda_u.pow(2) * -f101(2);
            let a = lambda_pow_2 - self.a - rhs.a;
            let b = lambda_u * (self.a - a) - self.b;

            G2Point::new(a, b)
        }
    }
}

impl Mul<F101> for G1Point {
    type Output = G1Point;
    fn mul(self, rhs: F101) -> Self::Output {
        let mut rhs = rhs.as_u64();
        if rhs == 0 || self.is_inf() {
            return G1Point::inf();
        }
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
            rhs = rhs >> 1;
            base = base + base;
        }
        result.unwrap()
    }
}

impl Mul<F101> for G2Point {
    type Output = G2Point;
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
            rhs = rhs >> 1;
            base = base + base;
        }
        result.unwrap()
    }
}

impl Mul<G2Point> for G2Point {
    type Output = G2Point;
    fn mul(self, rhs: G2Point) -> Self::Output {
        G2Point::new(
            self.a * rhs.a - f101(2) * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )
    }
}

fn pairing_f(r: u64, p: G1Point, q: G2Point) -> G2Point {
    // line equation from a to b point
    let l = |a: G1Point, b: G1Point| {
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
        pairing_f(r, p, q) * G2Point::new(q.a * x + c, q.b * y)
    } else {
        let r = r / 2;
        let (x, y, c) = l(p * f101(r), -p * f101(r) * f101(2));
        pairing_f(r, p, q).pow(2) * G2Point::new(q.a * x + c, q.b * y)
    }
}

pub fn pairing(g1: G1Point, g2: G2Point) -> G2Point {
    let p = 101u64;
    let r = 17u64;
    let k = 2u32;

    let exp = (p.pow(k) - 1) / r;

    pairing_f(17, g1, g2).pow(exp)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_g1_vectors() {
        let g = g1f(1, 2);
        let two_g = g + g;
        let four_g = two_g + two_g;
        let eight_g = four_g + four_g;
        let sixteen_g = eight_g + eight_g;

        assert_eq!(g1f(1, 99), -g);
        assert_eq!(g1f(68, 74), two_g);
        assert_eq!(g1f(68, 27), -two_g);
        assert_eq!(g1f(65, 98), four_g);
        assert_eq!(g1f(65, 3), -four_g);
        assert_eq!(g1f(18, 49), eight_g);
        assert_eq!(g1f(18, 52), -eight_g);
        assert_eq!(g1f(1, 99), sixteen_g);
        assert_eq!(g1f(1, 2), -sixteen_g);

        // since g = -16 g, this subgroup has order 17

        assert_eq!(g1f(26, 45), two_g + g);
        assert_eq!(g1f(12, 32), four_g + g);
        assert_eq!(g1f(18, 52), eight_g + g);
        assert_eq!(four_g + two_g, two_g + four_g);

        assert_eq!(g * f101(1), g);
        assert_eq!(g * f101(2), g + g);
        assert_eq!(g * f101(6), g + g + g + g + g + g);
    }
    #[test]
    fn test_g2_vectors() {
        let g = g2f(36, 31);

        // check point doubling
        assert_eq!(g2f(90, 82), g + g);

        // check point addition
        assert_eq!((g + g) + (g + g), g + g + g + g);

        // check point multiplication
        assert_eq!(g * f101(6), g + g + g + g + g + g);

        // // check frobenious map
        // assert_eq!(-g, g * 100);

        // check G2 multiplication
        assert_eq!(g2f(26, 97) * g2f(93, 76), g2f(97, 89));

        // check G2 exp
        assert_eq!(g2f(42, 49).pow(6), g2f(97, 89));
        assert_eq!(g2f(93, 76).pow(101), -g2f(93, 76));
        assert_eq!(g2f(93, 76).pow(102), (-g2f(93, 76)) * g2f(93, 76));
        assert_eq!(g2f(68, 47).pow(600), g2f(97, 89));
    }
}
