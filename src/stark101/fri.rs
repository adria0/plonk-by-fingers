use crate::{field::{Field, U64Field}, poly::Poly};

struct FriProof<F: Field> {
    evals: Vec<(F, F)>,
}


impl<F: Field> FriProof<F> {
    // Prove what we now the polinomial p
    fn prove<RNDS>(mut poly: Poly<F>, domain: &[F], rands: RNDS) -> Self
    where
        RNDS: IntoIterator<Item = F>,
    {
        let mut rands = rands.into_iter();
        let query = rands.next().unwrap().as_u64() as usize;

        let mut evals = Vec::new();
        let mut domain = domain.to_vec();

        while poly.degree() > 0 {
            let x = domain[query % domain.len()];
            evals.push((poly.eval(&x), poly.eval(&-x)));

            let mut g_x2 = vec![];
            let mut h_x2 = vec![];

            for (i, coeff) in poly.coeffs().iter().enumerate() {
                if i % 2 == 0 {
                    g_x2.push(*coeff);
                } else {
                    h_x2.push(*coeff);
                }
            }

            let g_x2 = Poly::new(g_x2);
            let h_x2 = Poly::new(h_x2);

            poly = g_x2 + h_x2 * Poly::new([rands.next().unwrap()]);

            domain = (0..domain.len() / 2)
                .map(|i| domain[i * 2])
                .collect::<Vec<_>>();
        }
        let x = domain[query % domain.len()];
        evals.push((poly.eval(&x), F::zero()));

        Self { evals }
    }

    // Verify low degree
    fn verify<RNDS>(&self, domain: &[F], rands: RNDS) -> bool
    where
        RNDS: IntoIterator<Item = F>,
    {
        let mut rands = rands.into_iter();
        let query = rands.next().unwrap().as_u64() as usize;
        let mut domain = domain.to_vec();
        for layer in 0..self.evals.len() - 1 {
            let x = domain[query % domain.len()];
            let rand = rands.next().unwrap();
            let check = ((rand + x) / (x * F::from(2))).unwrap() * self.evals[layer].0
                + ((rand - x) / (-x * F::from(2))).unwrap() * self.evals[layer].1;

            if self.evals[layer + 1].0 != check {
                return false;
            }

            domain = (0..domain.len() / 2)
                .map(|i| domain[i * 2])
                .collect::<Vec<_>>();
        }

        true
    }
}

#[test]
fn fri() {
    type FF = U64Field<3221225473>;
    let w: FF = 5.into();
    let g = w.pow(3 * 2u64.pow(20));
    let domain: Vec<_> = crate::mulmodg::MulGroupMod::<FF>::new(g).iter().take(1024).collect();

    let rands = std::iter::repeat(w);

    let proof = FriProof::prove(Poly::from([1,2,3,4,5,6,7,8]), &domain, rands.clone());
    assert!(proof.verify(&domain, rands));
}