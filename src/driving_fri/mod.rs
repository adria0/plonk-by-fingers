#![allow(non_snake_case, unused)]

use crate::field::Field;
use crate::poly::Poly;
use anyhow::Result;
use std::ops::Index;

#[derive(Clone, Debug)]
struct MulGroupMod<F: Field> {
    g: F,
    coset: F,
}

#[derive(Clone, Debug)]
struct MulGroupModIterator<F: Field> {
    g: F,
    n: F,
    coset: F,
}

// https://blog.lambdaclass.com/diving-deep-fri/

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

#[test]
fn fri_test() -> Result<()> {
    type F17 = crate::field::U64Field<17>;
    let f17 = |n: u64| F17::from(n);

    let trace_len = 4;

    // Problem statement & interpolation ---------------------------------

    let g = f17(13u64);
    let a = std::iter::successors(Some(f17(3)), |n| Some(n * n));
    let D_t = MulGroupMod::new(g);
    let t_x = D_t.lagrange(a.take(trace_len));

    println!("{}", t_x);

    assert_eq!(t_x.eval(&f17(1)), f17(3));
    assert_eq!(t_x.eval(&f17(13)), f17(9));
    assert_eq!(t_x.eval(&f17(16)), f17(13));
    assert_eq!(t_x.eval(&f17(4)), f17(16));

    // Committing to the trace polynomial -------------------------------

    let D_0 = MulGroupMod::new(f17(9)).coset(f17(3));

    let x_t_x = D_0.iter().map(|x| (x, t_x.eval(&x))).collect::<Vec<_>>();

    // Enter the constraints --------------------------------------------

    // Bundary constraint

    let p_1 = &t_x + Poly::from([-3i64]);
    let p_1_vanishing = Poly::z([f17(1)]);
    let C_1 = p_1 / p_1_vanishing;

    // Transmission verification constraint

    let t_gx = t_x.compose(&Poly::new(vec![F17::zero(), g]));

    let P_x = Poly::parse("x^2").unwrap().compose(&t_x);
    let P_y = Poly::parse("x").unwrap().compose(&t_gx);

    let p_2 = P_y - P_x;
    let p_2_vanishing = Poly::z(D_t.iter().take(trace_len - 1));
    let C_2 = p_2 / p_2_vanishing;

    // we want total degree must me a power of 2
    let max_degree = std::cmp::max(C_1.degree(), C_2.degree());
    let D = 2u64.pow((max_degree as f32).log2().ceil() as u32);

    let (alpha_1, beta_1, alpha_2, beta_2) = (f17(1), f17(3), f17(2), f17(4));

    let H_x = &C_1 * (Poly::y(alpha_1) * Poly::x().pow(D - C_1.degree()) + Poly::y(beta_1))
        + &C_2 * (Poly::y(alpha_2) * Poly::x().pow(D - C_2.degree()) + Poly::y(beta_2));

    // split H_x into even and odds, check that H(x) = H_1(x^2) + x·H_2(x^2)

    let (H_1_x, H_2_x) = {
        let H_1_coeffs =
            H_x.coeffs()
                .iter()
                .enumerate()
                .filter_map(|(i, n)| if i % 2 == 0 { Some(*n) } else { None });

        let H_2_coeffs =
            H_x.coeffs()
                .iter()
                .enumerate()
                .filter_map(|(i, n)| if i % 2 == 1 { Some(*n) } else { None });

        let H_1_x = Poly::new(H_1_coeffs);
        let H_2_x = Poly::new(H_2_coeffs);

        // check that H(x) = H_1(x^2) + x·H_2(x^2)
        assert_eq!(
            H_x,
            H_1_x.compose(&Poly::x().pow(2)) + Poly::x() * H_2_x.compose(&Poly::x().pow(2))
        );

        (H_1_x, H_2_x)
    };

    // Evaluate D0 in H1, H2
    let (H_1_eval, H_2_eval) = {
        let x_pow2 = D_0.iter().take(trace_len).map(|v| v * v);

        let H_1_eval = x_pow2.clone().map(|v| H_1_x.eval(&v));
        let H_2_eval = x_pow2.map(|v| H_2_x.eval(&v));
        (H_1_eval, H_2_eval)
    };

    assert_eq!(
        H_1_eval.collect::<Vec<_>>(),
        vec![f17(12), f17(13), f17(12), f17(13)]
    );
    assert_eq!(
        H_2_eval.collect::<Vec<_>>(),
        vec![f17(7), f17(10), f17(15), f17(12)]
    );

    // Sampling outside the original region ----------------------------------------

    // The prover sends to the verifier z
    let z = f17(8);
    let H_z = H_x.eval(&z);

    let pp = |s| Poly::<F17>::parse(s).unwrap();

    let H_1_z = H_1_x.eval(&(z * z));
    let H_2_z = H_2_x.eval(&(z * z));

    // Why does the verifier needs this? --------------------------------------------
    // The verifier has to check that the provided H_1(z^2) and H_2(z^2) matches the
    // calculation of H(z) from the trace elements

    let deep_1 = (t_x.clone() - Poly::y(t_x.eval(&z))) / Poly::new([-z, f17(1)]);
    let deep_2 = (t_x.clone() - Poly::y(t_x.eval(&(g * z)))) / Poly::new([-g * z, f17(1)]);

    // !! There is a mistake here: divided by x-z not x-z^2
    let deep_3 =
        (H_1_x.compose(&Poly::x().pow(2)) - H_1_x.eval(&(z * z))) / Poly::new([-z, f17(1)]);

    let deep_4 =
        (H_2_x.compose(&Poly::x().pow(2)) - H_2_x.eval(&(z * z))) / Poly::new([-z, f17(1)]);

    let gamma = [f17(1); 4];

    assert_eq!(deep_1, pp("13") * pp("5+16x+x^2"));
    assert_eq!(deep_2, pp("13") * pp("16+10x+x^2"));
    assert_eq!(deep_3, pp("15") * pp("15+x") * pp("8+x") * pp("2+x"));
    assert_eq!(deep_4, pp("9") * pp("8+x"));

    let P_0_x = deep_1 * gamma[0] + deep_2 * gamma[1] + deep_3 * gamma[2] + deep_4 * gamma[3];

    println!("P_0_x={}", P_0_x);

    println!("-----------------------------------");

    unreachable!()
}
