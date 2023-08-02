#![allow(dead_code)]

use crate::field::Field;
use crate::pairing::GT;
use std::fmt::Display;
use std::ops::Mul;
use std::ops::Neg;

use super::types::f101;
use super::types::F101;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GTP {
    a: F101,
    b: F101,
}

impl Display for GTP {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}+{}u", self.a, self.b)
    }
}

impl Neg for GTP {
    type Output = GTP;
    fn neg(self) -> Self::Output {
        GTP {
            a: self.a,
            b: -self.b,
        }
    }
}

impl GT for GTP {
    type S = u64;
    fn pow(&self, mut n: Self::S) -> Self {
        // frobenious map reduction: p^101 = -p
        let (mut p, mut base) = if n >= 101 {
            let base = -self.pow(n / 101);
            n %= 101;
            (base, *self)
        } else {
            (
                GTP {
                    a: F101::one(),
                    b: F101::zero(),
                },
                *self,
            )
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
impl Mul<GTP> for GTP {
    type Output = GTP;
    fn mul(self, rhs: GTP) -> Self::Output {
        GTP {
            a: self.a * rhs.a - f101(2) * self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a,
        }
    }
}

impl GTP {
    pub fn new(a: F101, b: F101) -> Self {
        Self { a, b }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn gtf(a: u64, b: u64) -> GTP {
        GTP {
            a: f101(a),
            b: f101(b),
        }
    }

    #[test]
    fn test_gt_vectors() {
        // check GT multiplication
        assert_eq!(gtf(26, 97) * gtf(93, 76), gtf(97, 89));

        // check GT exp
        assert_eq!(gtf(42, 49).pow(6), gtf(97, 89));
        assert_eq!(gtf(93, 76).pow(101), -gtf(93, 76));
        assert_eq!(gtf(93, 76).pow(102), (-gtf(93, 76)) * gtf(93, 76));
        assert_eq!(gtf(68, 47).pow(600), gtf(97, 89));
    }
}
