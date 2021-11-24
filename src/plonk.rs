#![allow(clippy::many_single_char_names)]

use crate::{
    ec::{GTPoint, Pairing},
    poly::Field,
};

use super::{
    constraints::*,
    ec::{G1Point, G2Point},
    matrix::Matrix,
    poly::Poly,
};

pub trait PlonkTypes: PartialEq {
    type GF: Field;
    type SF: Field;
    type G1: G1Point<S = Self::GF>;
    type G2: G2Point<S = Self::GF>;
    type GT: GTPoint;
    type E: Pairing<G1 = Self::G1, G2 = Self::G2, GT = Self::GT>;
    fn gf(sf: Self::SF) -> Self::GF;
}

pub struct SRS<P: PlonkTypes> {
    pub g1s: Vec<P::G1>,
    pub g2_1: P::G2,
    pub g2_s: P::G2,
}

impl<P: PlonkTypes> SRS<P> {
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
    pub fn eval_at_s(&self, vs: &Poly<P::SF>) -> P::G1 {
        vs.coeffs()
            .iter()
            .enumerate()
            .fold(G1Point::identity(), |acc, (n, v)| {
                acc + self.g1s[n] * P::gf(*v)
            })
    }
}

#[derive(Debug, PartialEq)]
pub struct Proof<P: PlonkTypes> {
    pub a_s: P::G1,
    pub b_s: P::G1,
    pub c_s: P::G1,
    pub z_s: P::G1,
    pub t_lo_s: P::G1,
    pub t_mid_s: P::G1,
    pub t_hi_s: P::G1,
    pub w_z_s: P::G1,
    pub w_z_omega_s: P::G1,
    pub a_z: P::SF,
    pub b_z: P::SF,
    pub c_z: P::SF,
    pub s_sigma_1_z: P::SF,
    pub s_sigma_2_z: P::SF,
    pub r_z: P::SF,
    pub z_omega_z: P::SF,
}

pub struct Challange<P: PlonkTypes> {
    pub alpha: P::SF,
    pub beta: P::SF,
    pub gamma: P::SF,
    pub z: P::SF,
    pub v: P::SF,
}

pub struct Plonk<P: PlonkTypes> {
    srs: SRS<P>,
    omega: P::SF,
    h: Vec<P::SF>,
    h_pows_inv: Matrix<P::SF>,
    k1: P::SF,
    k2: P::SF,
    k1_h: Vec<P::SF>,
    k2_h: Vec<P::SF>,
    z_h_x: Poly<P::SF>,
}

impl<P: PlonkTypes> Plonk<P> {
    pub fn new(srs: SRS<P>, omega: P::SF, omega_pows: usize, k1: P::SF, k2: P::SF) -> Self {
        // This roots of unity should be able to be generated through a generator
        // So the generator (called omega) creates these roots of unity (H)

        let h: Vec<_> = (0..omega_pows).map(|n| omega.pow(n as u64)).collect();

        // We need to label all of the values in our assignment with different field elements.
        // To do this, will use the roots of unity H along with two cosets of H.
        // We can get these cosets by multiplying each element of H by each of two constants: $k_1$ and $k_2$.
        // The constant $k_1$ is chosen so that it is not an element of $H$, and $k_2$ is chosen so that it is
        // neither an element of H nor $k_1H$. This ensures we have all the field elements to use as labels.
        // $H$ will be used to index $a$ values, $k_1H$ for $b$ values, $k_2H$ for $cC values

        assert!(!h.contains(&k1));
        assert!(!h.contains(&k2));

        let k1_h: Vec<_> = h.iter().map(|r| *r * k1).collect(); // k1_h is a coset of H

        assert!(!k1_h.contains(&k2));
        let k2_h: Vec<_> = h.iter().map(|r| *r * k2).collect(); // k2_h is a coset of H

        // In some point we will want to build polinomials that evaluates with specific values
        // at the roots of unity, for instance at the $a$ values:
        // $f_a(x)$ = polinomial that evaluates
        //   $f_a(w^0)$=a[0]
        //   $f_a(w^1)$=a[1]
        //   ...
        //   $f_a(w^n)$=a[n]
        //
        // we create here the h_pows_inv, an interpolation matrix, where the values
        // that we want to interpolate are multiplied (as a row matrix) with the
        // interpolation matrix to get the interpolation polinomial.

        let mut h_pows = Matrix::<P::SF>::zero(h.len(), h.len());
        for c in 0..h_pows.cols() {
            for r in 0..h_pows.rows() {
                h_pows[(r, c)] = h[r].pow(c as u64);
            }
        }

        let h_pows_inv = h_pows.inv();

        // The polynomial $Z_H$ is the polynomial that is zero on all the elements of
        // our subgroup $H$.

        let z_h_x = Poly::z(&h);

        Plonk {
            srs,
            omega,
            h,
            h_pows_inv,
            k1,
            k2,
            k1_h,
            k2_h,
            z_h_x,
        }
    }
    fn interpolate(&self, vv: &[P::SF]) -> Poly<P::SF> {
        (&self.h_pows_inv * Matrix::<_>::new(vv, vv.len(), 1)).into_poly()
    }
    fn copy_constraints_to_roots(&self, c: &[CopyOf]) -> Vec<P::SF> {
        c.iter()
            .map(|c| match c {
                CopyOf::A(n) => self.h[n - 1],
                CopyOf::B(n) => self.k1_h[n - 1],
                CopyOf::C(n) => self.k2_h[n - 1],
            })
            .collect::<Vec<_>>()
    }

    pub fn prove(
        &self,
        constraints: &Constrains<P::SF>,
        assigments: &Assigments<P::SF>,
        challange: &Challange<P>,
        rand: [P::SF; 9],
    ) -> Proof<P> {
        // check that the constraints satisfies the assigments
        assert!(constraints.satisfies(assigments));

        // In a non-interactive protocol, these values are derived from a
        // cryptographic hash of a transcript of the Prover's process. These values
        // would be unpredictable by the either the Prover or Verifier, but each party
        // is be able to compute the same challenge values using the same hash function
        // on the transcript. This is known as the Fiat-Shamir transform which can turn
        // an interactive protocol non-interactive.

        let Challange {
            alpha,
            beta,
            gamma,
            z,
            v,
        } = challange;

        let n = constraints.c_a.len() as u64;

        // create sigmas, mappings from the copy constraints to roots of unity elements
        // ---------------------------------------------------------------------------

        let sigma_1 = self.copy_constraints_to_roots(&constraints.c_a);
        let sigma_2 = self.copy_constraints_to_roots(&constraints.c_b);
        let sigma_3 = self.copy_constraints_to_roots(&constraints.c_c);

        // create a set of polinomials, those polinomials evaluates at roots of unity
        // for all the components of the plonk circuit
        //
        //   (a,b,c)                  : assigments / values
        //   (o,m,l,r,c)              : gate constraints
        //   (sigma1, sigma2, sigma3) : copy constraints
        // ---------------------------------------------------------------------------

        let f_a_x = self.interpolate(&assigments.a);
        let f_b_x = self.interpolate(&assigments.b);
        let f_c_x = self.interpolate(&assigments.c);
        let q_o_x = self.interpolate(&constraints.q_o);
        let q_m_x = self.interpolate(&constraints.q_m);
        let q_l_x = self.interpolate(&constraints.q_l);
        let q_r_x = self.interpolate(&constraints.q_r);
        let q_c_x = self.interpolate(&constraints.q_c);
        let s_sigma_1 = self.interpolate(&sigma_1);
        let s_sigma_2 = self.interpolate(&sigma_2);
        let s_sigma_3 = self.interpolate(&sigma_3);

        // round 1 - eval a(x), b(x), c(x) at s
        // ---------------------------------------------------------------------------

        let (b1, b2, b3, b4, b5, b6) = (rand[0], rand[1], rand[2], rand[3], rand[4], rand[5]);

        let a_x = Poly::new(vec![b2, b1]) * &self.z_h_x + f_a_x;
        let b_x = Poly::new(vec![b4, b3]) * &self.z_h_x + f_b_x;
        let c_x = Poly::new(vec![b6, b5]) * &self.z_h_x + f_c_x;

        // ouput of first step
        let a_s = self.srs.eval_at_s(&a_x);
        let b_s = self.srs.eval_at_s(&b_x);
        let c_s = self.srs.eval_at_s(&c_x);

        // round 2 - eval acummulator vector polinomial at s
        // ---------------------------------------------------------------------------
        // check https://vitalik.ca/general/2019/09/22/plonk.html
        //
        // alpha, beta and gamma are "random" that comes from the H(transcript)
        // b7, b8, b9 are random
        //

        let (b7, b8, b9) = (rand[6], rand[7], rand[8]);

        // create the accumulator vector, it has the property that is going to have
        // the same value than
        //

        let mut acc = vec![P::SF::one()];
        for i in 1..n as usize {
            let a = assigments.a[i - 1];
            let b = assigments.b[i - 1];
            let c = assigments.c[i - 1];

            // omega_pow is the root of unity
            let omega_pow = self.omega.pow((i as u64) - 1);

            // combine each wire value with its index position
            let dend = (a + *beta * omega_pow + gamma)
                * (b + *beta * self.k1 * omega_pow + gamma)
                * (c + *beta * self.k2 * omega_pow + gamma);

            // combine each wire value with its *permuted* index position
            let dsor = (a + *beta * s_sigma_1.eval(&omega_pow) + gamma)
                * (b + *beta * s_sigma_2.eval(&omega_pow) + gamma)
                * (c + *beta * s_sigma_3.eval(&omega_pow) + gamma);

            let v = acc[i - 1] * (dend / dsor).unwrap();
            acc.push(v);
        }

        // the accumulator vector is interpolated into a polynomial acc(x)
        // and this is used to create the polynomial z.

        let acc_x = self.interpolate(&acc);

        let z_x = Poly::new(vec![b9, b8, b7]) * &self.z_h_x + acc_x;

        // Then z is evaluated at the secret number s using the SRS from the setup.
        // output of second step
        let z_s = self.srs.eval_at_s(&z_x);

        // round 3 - compute the quotient polynomial
        // ---------------------------------------------------------------------------
        // Next comes the most massive computation of the entire protocol.
        // Our goal is to compute the polynomial t, which will be of degree 3n+5 for n gates.
        // The polynomial t encodes the majority of the information contained in our circuit
        // and assignments all at once.

        // L1 refers to a Lagrange basis polynomial over our roots of unity H.
        // Specifically L1(1)=1, but takes the value 0 on each of the other roots of unity.
        // The coefficients for L1 can be found by interpolating the vector

        let lagrange_vector: Vec<_> = std::iter::once(P::SF::one())
            .chain(std::iter::repeat(P::SF::zero()))
            .take(self.h.len())
            .collect();

        let l_1_x = self.interpolate(&lagrange_vector);

        // public inputs
        let p_i_x = P::SF::zero().as_poly(); // Public input

        // helper variables
        let a_x_b_x_q_m_x = (&a_x * &b_x) * &q_m_x;
        let a_x_q_l_x = &a_x * &q_l_x;
        let b_x_q_r_x = &b_x * &q_r_x;
        let c_x_q_o_x = &c_x * &q_o_x;
        let alpha_a_x_beta_x_gamma = &(&a_x + &Poly::new(vec![*gamma, *beta])) * alpha;
        let b_x_beta_k1_x_gamma = &b_x + &Poly::new(vec![*gamma, *beta * self.k1]);
        let c_x_beta_k2_x_gamma = &c_x + &Poly::new(vec![*gamma, *beta * self.k2]);
        let z_omega_x = Poly::new(
            z_x.coeffs()
                .iter()
                .enumerate()
                .map(|(n, c)| *c * self.omega.pow(n as u64))
                .collect::<Vec<_>>(),
        );
        let alpha_a_x_beta_s_sigma1_x_gamma = (&a_x + &s_sigma_1 * beta + gamma) * alpha;
        let b_x_beta_s_sigma2_x_gamma = &b_x + &s_sigma_2 * beta + gamma;
        let c_x_beta_s_sigma3_x_gamma = &c_x + &s_sigma_3 * beta + gamma;
        let alpha_2_z_x_1_l_1_x = ((&z_x + (-P::SF::one()).as_poly()) * alpha.pow(2)) * &l_1_x;

        // compute t(x)

        let t_1_z_h = a_x_b_x_q_m_x + a_x_q_l_x + b_x_q_r_x + c_x_q_o_x + p_i_x + &q_c_x;
        let t_2_z_h = alpha_a_x_beta_x_gamma * b_x_beta_k1_x_gamma * c_x_beta_k2_x_gamma * &z_x;
        let t_3_z_h = alpha_a_x_beta_s_sigma1_x_gamma
            * b_x_beta_s_sigma2_x_gamma
            * c_x_beta_s_sigma3_x_gamma
            * &z_omega_x;

        let t_4_z_h = alpha_2_z_x_1_l_1_x;

        let (t_x, rem) = (t_1_z_h + t_2_z_h - t_3_z_h + t_4_z_h) / self.z_h_x.clone();
        assert_eq!(rem, Poly::zero());

        // It turns out that for n constraints, t will have degree 3n+5, which is too large to use the SRS
        // from the setup phase in Part 1. However we can break t into three parts of degree n+1 each.
        // Each of these polynomials will use n+2 coefficients from t(x)

        let t_hi_x = Poly::new(t_x.coeffs()[12..18].to_owned());
        let t_mid_x = Poly::new(t_x.coeffs()[6..12].to_owned());
        let t_lo_x = Poly::new(t_x.coeffs()[0..6].to_owned());

        // After dividing into three parts, we evaluate each part at s using the SRS.
        // output of the third step

        let t_hi_s = self.srs.eval_at_s(&t_hi_x);
        let t_mid_s = self.srs.eval_at_s(&t_mid_x);
        let t_lo_s = self.srs.eval_at_s(&t_lo_x);

        // round 4 - compute linearization polynomial
        // ---------------------------------------------------------------------------
        // we create a polynomial r that is a kind of partner to t, where many of the polynomials that were
        // included in t are replaced by field elements that are evaluations of those polynomials at a
        // challenge value 𝔷.

        let a_z = a_x.eval(z);
        let b_z = b_x.eval(z);
        let c_z = c_x.eval(z);
        let s_sigma_1_z = s_sigma_1.eval(z);
        let s_sigma_2_z = s_sigma_2.eval(z);
        let t_z = t_x.eval(z);
        let z_omega_z = z_omega_x.eval(z);

        let a_z_b_z_q_m_x = &q_m_x * a_z * b_z;

        let a_z_q_l_x = &q_l_x * a_z;
        let b_z_q_r_x = &q_r_x * b_z;
        let c_z_q_o_x = &q_o_x * c_z;

        let r_1_x = a_z_b_z_q_m_x + a_z_q_l_x + b_z_q_r_x + c_z_q_o_x + q_c_x;
        let r_2_x = &z_x
            * ((a_z + *beta * z + gamma)
                * (b_z + *beta * self.k1 * z + gamma)
                * (c_z + *beta * self.k2 * z + gamma)
                * alpha);

        let r_3_x = &z_x
            * (s_sigma_3 * *beta * z_omega_z)
            * ((a_z + *beta * s_sigma_1_z + gamma) * (b_z + *beta * s_sigma_2_z + gamma) * alpha);
        let r_4_x = &z_x * l_1_x.eval(z) * alpha.pow(2);

        let r_x = r_1_x + r_2_x + r_3_x + r_4_x;

        // compute linearization evaluation
        let r_z = r_x.eval(z);

        // round 5
        // ---------------------------------------------------------------------------
        // we create two large polynomials that combine all the polynomials we've been
        // using so far and we output commitments to them.

        // compute opening proof polynomial w_z_x
        let w_z_x = (t_lo_x + t_mid_x * z.pow(n + 2) + t_hi_x * z.pow(2 * n + 4) - t_z)
            + (r_x - r_z) * v
            + (a_x - a_z) * v.pow(2)
            + (b_x - b_z) * v.pow(3)
            + (c_x - c_z) * v.pow(4)
            + (s_sigma_1 - s_sigma_1_z) * v.pow(5)
            + (s_sigma_2 - s_sigma_2_z) * v.pow(6);
        let (w_z_x, rem) = w_z_x / Poly::new(vec![-*z, P::SF::one()]);
        assert_eq!(rem, Poly::zero());

        // compute opening proof polinomial w_zw_x
        let (w_z_omega_x, rem) =
            (z_x - z_omega_z) / Poly::new(vec![-*z * self.omega, P::SF::one()]);
        assert_eq!(rem, Poly::zero());

        // compute opening proof polinomials at s
        let w_z_s = self.srs.eval_at_s(&w_z_x);
        let w_z_omega_s = self.srs.eval_at_s(&w_z_omega_x);

        Proof {
            a_s,
            b_s,
            c_s,
            z_s,
            t_lo_s,
            t_mid_s,
            t_hi_s,
            w_z_s,
            w_z_omega_s,
            a_z,
            b_z,
            c_z,
            s_sigma_1_z,
            s_sigma_2_z,
            r_z,
            z_omega_z,
        }
    }

    pub fn verify(
        &self,
        constraints: &Constrains<P::SF>,
        proof: &Proof<P>,
        challange: &Challange<P>,
        rand: [P::SF; 1],
    ) -> bool {
        let Proof {
            a_s,
            b_s,
            c_s,
            z_s,
            t_lo_s,
            t_mid_s,
            t_hi_s,
            w_z_s,
            w_z_omega_s,
            a_z,
            b_z,
            c_z,
            s_sigma_1_z,
            s_sigma_2_z,
            r_z,
            z_omega_z,
        } = proof;

        let Challange {
            alpha,
            beta,
            gamma,
            z,
            v,
        } = challange;

        // verifier preprocessing
        // ---------------------------------------------------------------------------
        let sigma_1 = self.copy_constraints_to_roots(&constraints.c_a);
        let sigma_2 = self.copy_constraints_to_roots(&constraints.c_b);
        let sigma_3 = self.copy_constraints_to_roots(&constraints.c_c);

        let q_m_s = self.srs.eval_at_s(&self.interpolate(&constraints.q_m));
        let q_l_s = self.srs.eval_at_s(&self.interpolate(&constraints.q_l));
        let q_r_s = self.srs.eval_at_s(&self.interpolate(&constraints.q_r));
        let q_o_s = self.srs.eval_at_s(&self.interpolate(&constraints.q_o));
        let q_c_s = self.srs.eval_at_s(&self.interpolate(&constraints.q_c));
        let sigma_1_s = self.srs.eval_at_s(&self.interpolate(&sigma_1));
        let sigma_2_s = self.srs.eval_at_s(&self.interpolate(&sigma_2));
        let sigma_3_s = self.srs.eval_at_s(&self.interpolate(&sigma_3));

        let u = rand[0];

        // Step 1. Validate proof points in G1

        if !a_s.in_curve()
            || !b_s.in_curve()
            || !c_s.in_curve()
            || !z_s.in_curve()
            || !t_lo_s.in_curve()
            || !t_mid_s.in_curve()
            || !t_hi_s.in_curve()
            || !w_z_s.in_curve()
            || !w_z_omega_s.in_curve()
        {
            return false;
        }

        // Step 2. Validate proof fields in SF

        if !a_z.in_field()
            || !b_z.in_field()
            || !c_z.in_field()
            || !s_sigma_1_z.in_field()
            || !s_sigma_2_z.in_field()
            || !r_z.in_field()
            || !z_omega_z.in_field()
        {
            return false;
        }

        // Step 3. We do not have public inputs, nothing to do

        // Step 4. Evaluate z_h at z

        let z_h_z = self.z_h_x.eval(z);

        // Step 5. Evaluate lagrange on z

        let lagrange_vector: Vec<_> = std::iter::once(P::SF::one())
            .chain(std::iter::repeat(P::SF::zero()))
            .take(self.h.len())
            .collect();

        let l_1_z = self.interpolate(&lagrange_vector).eval(z);

        // Step 6. We do not have public inputs, nothing to do

        let p_i_z = P::SF::zero();

        // Step 7. Compute quotient polinomial evaluation

        let a_z_beta_s_sigma_1_z_gamma = *beta * s_sigma_1_z + gamma + a_z;

        let b_z_beta_s_sigma_2_z_gamma = *beta * s_sigma_2_z + gamma + b_z;
        let c_z_gamma = *c_z + gamma;
        let l_1_z_alpha_2 = l_1_z * alpha.pow(2);
        let t_z = ((*r_z + p_i_z
            - (a_z_beta_s_sigma_1_z_gamma * b_z_beta_s_sigma_2_z_gamma * c_z_gamma * z_omega_z)
            - l_1_z_alpha_2)
            / z_h_z)
            .unwrap();

        // Step 8. Compute the first part of batched polinomial commitment

        let d_1_s = q_m_s * P::gf(*a_z * b_z * v)
            + q_l_s * P::gf(*a_z * v)
            + q_r_s * P::gf(*b_z * v)
            + q_o_s * P::gf(*c_z * v)
            + q_c_s * P::gf(*v);

        let d_2_s = *z_s
            * P::gf(
                (*a_z + *beta * z + gamma)
                    * (*b_z + *beta * self.k1 * z + gamma)
                    * (*c_z + *beta * self.k2 * z + gamma)
                    * alpha
                    * v
                    + l_1_z * alpha.pow(2) * v
                    + u,
            );

        let d_3_s = sigma_3_s
            * P::gf(
                (*a_z + *beta * s_sigma_1_z + gamma)
                    * (*b_z + *beta * s_sigma_2_z + gamma)
                    * alpha
                    * v
                    * beta
                    * z_omega_z,
            );

        let d_s = d_1_s + d_2_s + -d_3_s;

        // Step 9. Compute full batched polinomial commitment

        let n = constraints.c_a.len() as u64;

        let f_s = t_lo_s.to_owned()
            + t_mid_s.to_owned() * P::gf(z.pow(n + 2))
            + t_hi_s.to_owned() * P::gf(z.pow(2 * n + 4))
            + d_s
            + a_s.to_owned() * P::gf(v.pow(2))
            + b_s.to_owned() * P::gf(v.pow(3))
            + c_s.to_owned() * P::gf(v.pow(4))
            + sigma_1_s.to_owned() * P::gf(v.pow(5))
            + sigma_2_s.to_owned() * P::gf(v.pow(6));

        // Step 10. Compute group encoded batch evaluation [E]

        let e_s = self.srs.eval_at_s(&Poly::from(&[1]))
            * P::gf(
                t_z + *v * r_z
                    + v.pow(2) * a_z
                    + v.pow(3) * b_z
                    + v.pow(4) * c_z
                    + v.pow(5) * s_sigma_1_z
                    + v.pow(6) * s_sigma_2_z
                    + u * z_omega_z,
            );

        // Step 11. Batch validate all equations

        let e_1_q1 = *w_z_s + *w_z_omega_s * P::gf(u);
        let e_1_q2 = self.srs.g2_s;
        let e_2_q1 = *w_z_s * P::gf(*z) + *w_z_omega_s * P::gf(u * z * self.omega) + f_s + -e_s;
        let e_2_q2 = self.srs.g2_1;

        let e_1 = P::E::pairing(e_1_q1, e_1_q2);
        let e_2 = P::E::pairing(e_2_q1, e_2_q2);

        e_1 == e_2
    }
}
