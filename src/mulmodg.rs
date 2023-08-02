use crate::field::Field;
use crate::poly::Poly;

#[derive(Clone, Debug)]
pub struct MulGroupMod<F: Field> {
    g: F,
    coset: F,
}

#[derive(Clone, Debug)]
pub struct MulGroupModIterator<F: Field> {
    g: F,
    n: F,
    coset: F,
}

impl<F: Field> MulGroupMod<F> {
    pub fn new(g: F) -> Self {
        Self { g, coset: F::one() }
    }
    pub fn lagrange<I: IntoIterator<Item = F>>(&self, values: I) -> Poly<F> {
        let y = values.into_iter();
        let points: Vec<(F, F)> = self.iter().zip(y).collect::<Vec<_>>();
        Poly::lagrange(&points)
    }
    pub fn coset(&self, coset: F) -> MulGroupMod<F> {
        MulGroupMod { g: self.g, coset }
    }
    pub fn at(&self, p: u64) -> F {
        self.coset * self.g.pow(p)
    }
    pub fn iter(&self) -> MulGroupModIterator<F> {
        MulGroupModIterator {
            g: self.g,
            n: F::zero(),
            coset: self.coset,
        }
    }
}

impl<F: Field> Iterator for MulGroupModIterator<F> {
    type Item = F;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.n == F::zero() {
            self.n = self.g;
            Some(self.coset)
        } else {
            let current = self.n;
            self.n = self.n * self.g;

            if self.n == self.g {
                None
            } else {
                Some(self.coset * current)
            }
        }
    }
}
