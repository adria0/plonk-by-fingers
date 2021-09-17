#![allow(dead_code, unused_imports)]

use super::constrains::*;
use super::ec::{g1f, g2f, pairing, G1Point, G2Point};
use super::field::Field;
use super::matrix::Matrix;
use super::plonk::{f101, f17, p17, F101, F17, P17};
use super::poly::Poly;

pub struct SRS {
    pub g1s: Vec<G1Point>,
    pub g2_1: G2Point,
    pub g2_s: G2Point,
}

impl SRS {
    pub fn create(s: F101, n: usize, g1: G1Point, g2: G2Point) -> SRS {
        let mut g1s = Vec::new();
        let mut s_pow = s;
        g1s.push(g1);
        for _ in 0..n {
            g1s.push(g1 * s_pow);
            s_pow = s_pow * s;
        }
        Self {
            g1s,
            g2_1: g2,
            g2_s: g2 * s,
        }
    }

    // evaluate a polinomial at secret point s using SRS G1s
    pub fn eval_at_s(&self, vs: &Poly<17>) -> G1Point {
        vs.coeffs()
            .iter()
            .enumerate()
            .fold(G1Point::inf(), |acc, (n, v)| {
                acc + self.g1s[n] * v.rebase::<101>()
            })
    }
}

#[derive(Debug, PartialEq)]
pub struct Proof {
    a_s: G1Point,
    b_s: G1Point,
    c_s: G1Point,
    z_s: G1Point,
    t_lo_s: G1Point,
    t_mid_s: G1Point,
    t_hi_s: G1Point,
    w_z_s: G1Point,
    w_z_omega_s: G1Point,
    a_z: Field<17>,
    b_z: Field<17>,
    c_z: Field<17>,
    s_sigma_1_z: Field<17>,
    s_sigma_2_z: Field<17>,
    r_z: Field<17>,
    z_omega_z: Field<17>,
}

pub struct Challange {
    pub alpha: F17,
    pub beta: F17,
    pub gamma: F17,
    pub z: F17,
    pub v: F17,
}

struct Plonk {
    srs: SRS,
    omega: F17,
    h: [F17; 4],
    h_pows_inv: Matrix<17>,
    k1: F17,
    k2: F17,
    k1_h: Vec<F17>,
    k2_h: Vec<F17>,
    z_h_x: P17, // Zh(x) = x^4 - 1
}

impl Plonk {
    pub fn new(srs: SRS) -> Plonk {
        let omega = f17(4);

        // roots of unity for x^omega = 1
        let h = [f17(1), f17(4), f17(16), f17(13)];

        let mut h_pows = Matrix::<17>::zero(4, 4);
        for c in 0..h_pows.cols() {
            for r in 0..h_pows.rows() {
                h_pows[(r, c)] = h[r].pow(c as u64);
            }
        }
        let h_pows_inv = h_pows.inv();

        let (k1, k2) = (f17(2), f17(3));
        assert!(!h.contains(&k1));
        assert!(!h.contains(&k2));

        let k1_h: Vec<F17> = h.iter().map(|r| *r * k1).collect();
        let k2_h: Vec<F17> = h.iter().map(|r| *r * k2).collect();

        let z_h_x = p17(&[-1, 0, 0, 0, 1]); // Zh(x) = x^4 - 1

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
    fn interpolate(&self, vv: &[Field<17>]) -> P17 {
        (&self.h_pows_inv * Matrix::<17>::new(&vv, vv.len(), 1)).into_poly()
    }
    fn copyc_to_roots(&self, c: &[CopyOf]) -> Vec<F17> {
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
        constrains: &Constrains,
        assigments: &Assigments,
        challange: &Challange,
        rand: [F17; 9],
    ) -> Proof {
        // check that the constrains satisfies the assigments
        assert!(constrains.satisfies(&assigments));

        let Challange {
            alpha,
            beta,
            gamma,
            z,
            v,
        } = challange;

        let n = constrains.c_a.len() as u64;

        // now I have a generator g=(1,2) on Field::<M>, that generates a subgroup F17
        // we need at least <gates_len> nth-roots-of-unity from F17
        // fortunately <gates_len>=4 and there are 4 elements in F17 that x^4 = 1
        // roots_of_unity are H={1,4,16,13}
        // create sigmas, mappings from the copy constrains to roots of unity elements
        // ---------------------------------------------------------------------------

        let sigma_1 = self.copyc_to_roots(&constrains.c_a);
        let sigma_2 = self.copyc_to_roots(&constrains.c_b);
        let sigma_3 = self.copyc_to_roots(&constrains.c_c);

        // prepare vectors as polys
        // ---------------------------------------------------------------------------

        let f_a_x = self.interpolate(&assigments.a);
        let f_b_x = self.interpolate(&assigments.b);
        let f_c_x = self.interpolate(&assigments.c);
        let q_o_x = self.interpolate(&constrains.q_o);
        let q_m_x = self.interpolate(&constrains.q_m);
        let q_l_x = self.interpolate(&constrains.q_l);
        let q_r_x = self.interpolate(&constrains.q_r);
        let q_c_x = self.interpolate(&constrains.q_c);
        let s_sigma_1 = self.interpolate(&sigma_1);
        let s_sigma_2 = self.interpolate(&sigma_2);
        let s_sigma_3 = self.interpolate(&sigma_3);

        // round 1
        // ---------------------------------------------------------------------------

        let (b1, b2, b3, b4, b5, b6) = (rand[0], rand[1], rand[2], rand[3], rand[4], rand[5]);

        let a_x = P17::new(vec![b2, b1]) * &self.z_h_x + f_a_x;
        let b_x = P17::new(vec![b4, b3]) * &self.z_h_x + f_b_x;
        let c_x = P17::new(vec![b6, b5]) * &self.z_h_x + f_c_x;

        let a_s = self.srs.eval_at_s(&a_x);
        let b_s = self.srs.eval_at_s(&b_x);
        let c_s = self.srs.eval_at_s(&c_x);

        // round 2
        // ---------------------------------------------------------------------------

        let (b7, b8, b9) = (rand[6], rand[7], rand[8]);

        let mut acc = vec![f17(1)];
        for i in 1..n as usize {
            let a = assigments.a[i - 1];
            let b = assigments.b[i - 1];
            let c = assigments.c[i - 1];
            let omega_pow = self.omega.pow((i as u64) - 1);
            let dend = (a + beta * omega_pow + gamma)
                * (b + beta * self.k1 * omega_pow + gamma)
                * (c + beta * self.k2 * omega_pow + gamma);
            let dsor = (a + beta * s_sigma_1.eval(&omega_pow) + gamma)
                * (b + beta * s_sigma_2.eval(&omega_pow) + gamma)
                * (c + beta * s_sigma_3.eval(&omega_pow) + gamma);

            let v = acc[i - 1] * (dend / dsor).unwrap();
            acc.push(v);
        }

        let acc_x = self.interpolate(&acc);
        let z_x = P17::new(vec![b9, b8, b7]) * &self.z_h_x + acc_x;
        let z_s = self.srs.eval_at_s(&z_x);

        // round 3
        // ---------------------------------------------------------------------------

        let l_1_x = self.interpolate(&[f17(1), f17(0), f17(0), f17(0)]);
        let p_i_x = f17(0).as_poly(); // Public input

        let a_x_b_x_q_m_x = (&a_x * &b_x) * &q_m_x;
        let a_x_q_l_x = &a_x * &q_l_x;
        let b_x_q_r_x = &b_x * &q_r_x;
        let c_x_q_o_x = &c_x * &q_o_x;
        let alpha_a_x_beta_x_gamma = &(&a_x + &Poly::new(vec![*gamma, *beta])) * alpha;
        let b_x_beta_k1_x_gamma = &b_x + &Poly::new(vec![*gamma, beta * self.k1]);
        let c_x_beta_k2_x_gamma = &c_x + &Poly::new(vec![*gamma, beta * self.k2]);
        let z_omega_x = Poly::new(
            z_x.coeffs()
                .iter()
                .enumerate()
                .map(|(n, c)| c * &self.omega.pow(n as u64))
                .collect::<Vec<_>>(),
        );
        let alpha_a_x_beta_s_sigma1_x_gamma = (&a_x + &s_sigma_1 * beta + gamma) * alpha;
        let b_x_beta_s_sigma2_x_gamma = &b_x + &s_sigma_2 * beta + gamma;
        let c_x_beta_s_sigma3_x_gamma = &c_x + &s_sigma_3 * beta + gamma;
        let alpha_2_z_x_1_l_1_x = ((&z_x + (-f17(1)).as_poly()) * alpha.pow(2)) * &l_1_x;

        let t_1_z_h = a_x_b_x_q_m_x + a_x_q_l_x + b_x_q_r_x + c_x_q_o_x + p_i_x + &q_c_x;
        let t_2_z_h = alpha_a_x_beta_x_gamma * b_x_beta_k1_x_gamma * c_x_beta_k2_x_gamma * &z_x;
        let t_3_z_h = alpha_a_x_beta_s_sigma1_x_gamma
            * b_x_beta_s_sigma2_x_gamma
            * c_x_beta_s_sigma3_x_gamma
            * &z_omega_x;

        let t_4_z_h = alpha_2_z_x_1_l_1_x;

        let (t_x, rem) = (t_1_z_h + t_2_z_h - t_3_z_h + t_4_z_h) / self.z_h_x.clone();
        assert_eq!(rem, Poly::zero());

        let t_hi_x = Poly::new(t_x.coeffs()[2 * 6..3 * 6].to_owned());
        let t_mid_x = Poly::new(t_x.coeffs()[1 * 6..2 * 6].to_owned());
        let t_lo_x = Poly::new(t_x.coeffs()[0 * 6..1 * 6].to_owned());

        let t_hi_s = self.srs.eval_at_s(&t_hi_x);
        let t_mid_s = self.srs.eval_at_s(&t_mid_x);
        let t_lo_s = self.srs.eval_at_s(&t_lo_x);

        // round 4
        // ---------------------------------------------------------------------------

        let a_z = a_x.eval(&z);
        let b_z = b_x.eval(&z);
        let c_z = c_x.eval(&z);
        let s_sigma_1_z = s_sigma_1.eval(&z);
        let s_sigma_2_z = s_sigma_2.eval(&z);
        let t_z = t_x.eval(&z);
        let z_omega_z = z_omega_x.eval(&z);

        let a_z_b_z_q_m_x = &q_m_x * &a_z * &b_z;

        let a_z_q_l_x = &q_l_x * &a_z;
        let b_z_q_r_x = &q_r_x * &b_z;
        let c_z_q_o_x = &q_o_x * &c_z;

        let r_1_x = a_z_b_z_q_m_x + a_z_q_l_x + b_z_q_r_x + c_z_q_o_x + q_c_x;
        let r_2_x = &z_x
            * ((a_z + beta * z + gamma)
                * (b_z + beta * self.k1 * z + gamma)
                * (c_z + beta * self.k2 * z + gamma)
                * alpha);
        let r_3_x = &z_x
            * ((a_z + beta * s_sigma_1_z + gamma)
                * (b_z + beta * s_sigma_2_z + gamma)
                * (beta * z_omega_z * s_sigma_3)
                * alpha);
        let r_4_x = &z_x * l_1_x.eval(&z) * alpha.pow(2);

        let r_x = r_1_x + r_2_x + r_3_x + r_4_x;
        let r_z = r_x.eval(&z);

        // round 5
        // ---------------------------------------------------------------------------

        let v = v;
        let w_z_x = (t_lo_x + t_mid_x * z.pow(n + 2) + t_hi_x * z.pow(2 * n + 4) - t_z)
            + v * (r_x - r_z)
            + v.pow(2) * (a_x - a_z)
            + v.pow(3) * (b_x - b_z)
            + v.pow(4) * (c_x - c_z)
            + v.pow(5) * (s_sigma_1 - s_sigma_1_z)
            + v.pow(6) * (s_sigma_2 - s_sigma_2_z);
        let (w_z_x, rem) = w_z_x / Poly::new(vec![-z, f17(1)]);
        assert_eq!(rem, Poly::zero());

        let (w_z_omega_x, rem) = (z_x - z_omega_z) / Poly::new(vec![-z * self.omega, f17(1)]);
        assert_eq!(rem, Poly::zero());

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

    fn verify(
        &self,
        constrains: &Constrains,
        proof: &Proof,
        challange: &Challange,
        rand: [F17; 1],
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
        let sigma_1 = self.copyc_to_roots(&constrains.c_a);
        let sigma_2 = self.copyc_to_roots(&constrains.c_b);
        let sigma_3 = self.copyc_to_roots(&constrains.c_c);

        let q_m_s = self.srs.eval_at_s(&self.interpolate(&constrains.q_m));
        let q_l_s = self.srs.eval_at_s(&self.interpolate(&constrains.q_l));
        let q_r_s = self.srs.eval_at_s(&self.interpolate(&constrains.q_r));
        let q_o_s = self.srs.eval_at_s(&self.interpolate(&constrains.q_o));
        let q_c_s = self.srs.eval_at_s(&self.interpolate(&constrains.q_c));
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

        // Step 2. Validate proof fields in F17

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

        let z_h_z = self.z_h_x.eval(&z);

        // Step 5. Evaluate lagrange on z

        let l_1_z = self.interpolate(&[f17(1), f17(0), f17(0), f17(0)]).eval(&z);

        // Step 6. We do not have public inputs, nothing to do

        let p_i_z = F17::zero();

        // Step 7. Compute quotient polinomial evaluation

        let a_z_beta_s_sigma_1_z_gamma = a_z + beta * s_sigma_1_z + gamma;
        let b_z_beta_s_sigma_2_z_gamma = b_z + beta * s_sigma_2_z + gamma;
        let c_z_gamma = c_z + gamma;
        let l_1_z_alpha_2 = l_1_z * alpha.pow(2);
        let t_z = ((r_z + p_i_z
            - (a_z_beta_s_sigma_1_z_gamma * b_z_beta_s_sigma_2_z_gamma * c_z_gamma * z_omega_z)
            - l_1_z_alpha_2)
            / z_h_z)
            .unwrap();

        // Step 8. Compute the first part of batched polinomial commitment

        let d_1_s = q_m_s * (a_z * b_z * v).rebase::<101>()
            + q_l_s * (a_z * v).rebase::<101>()
            + q_r_s * (b_z * v).rebase::<101>()
            + q_o_s * (c_z * v).rebase::<101>()
            + q_c_s * (v).rebase::<101>();

        let d_2_s = *z_s
            * ((a_z + beta * z + gamma)
                * (b_z + beta * self.k1 * z + gamma)
                * (c_z + beta * self.k2 * z + gamma)
                * alpha
                * v
                + l_1_z * alpha.pow(2) * v
                + u)
                .rebase::<101>();

        let d_3_s = sigma_3_s
            * ((a_z + beta * s_sigma_1_z + gamma)
                * (b_z + beta * s_sigma_2_z + gamma)
                * alpha
                * v
                * beta
                * z_omega_z)
                .rebase::<101>();

        let d_s = d_1_s + d_2_s + -d_3_s;

        // Step 9. Compute full batched polinomial commitment

        let n = constrains.c_a.len() as u64;

        let f_s = t_lo_s.to_owned()
            + t_mid_s.to_owned() * z.pow(n + 2).rebase::<101>()
            + t_hi_s.to_owned() * z.pow(2 * n + 4).rebase::<101>()
            + d_s
            + a_s.to_owned() * v.pow(2).rebase::<101>()
            + b_s.to_owned() * v.pow(3).rebase::<101>()
            + c_s.to_owned() * v.pow(4).rebase::<101>()
            + sigma_1_s.to_owned() * v.pow(5).rebase::<101>()
            + sigma_2_s.to_owned() * v.pow(6).rebase::<101>();

        // Step 10. Compute group encoded batch evaluation [E]

        let e_s = self.srs.eval_at_s(&p17(&[1]))
            * (t_z
                + v * r_z
                + v.pow(2) * a_z
                + v.pow(3) * b_z
                + v.pow(4) * c_z
                + v.pow(5) * s_sigma_1_z
                + v.pow(6) * s_sigma_2_z
                + u * z_omega_z)
                .rebase::<101>();

        // Step 11. Batch validate all equations

        let e_1_q1 = *w_z_s + *w_z_omega_s * u.rebase::<101>();
        let e_1_q2 = self.srs.g2_s;
        let e_2_q1 = *w_z_s * z.rebase::<101>()
            + *w_z_omega_s * (u * z * self.omega).rebase::<101>()
            + f_s
            + -e_s;
        let e_2_q2 = self.srs.g2_1;

        let e_1 = pairing(e_1_q1, e_1_q2);
        let e_2 = pairing(e_2_q1, e_2_q2);

        println!("e1({},{}) = {}", e_1_q1, e_1_q2, e_1);
        println!("e2({},{}) = {}", e_2_q1, e_2_q2, e_2);

        e_1 == e_2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_plonk_gen_proof() {
        // create the trusted setup
        let s = 2; // the toxic waste
        let srs = SRS::create(f101(s), 6, g1f(1, 2), g2f(36, 31));

        let plonk = Plonk::new(srs);

        // constrains and assigments
        let constrains = Constrains::new(
            &[
                Gate::mul_a_b(),
                Gate::mul_a_b(),
                Gate::mul_a_b(),
                Gate::sum_a_b(),
            ],
            (
                vec![CopyOf::B(1), CopyOf::B(2), CopyOf::B(3), CopyOf::C(1)],
                vec![CopyOf::A(1), CopyOf::A(2), CopyOf::A(3), CopyOf::C(2)],
                vec![CopyOf::A(4), CopyOf::B(4), CopyOf::C(4), CopyOf::C(3)],
            ),
        );

        let assigments = Assigments::new(&[
            Assigment::new(Field::from(3), Field::from(3), Field::from(9)),
            Assigment::new(Field::from(4), Field::from(4), Field::from(16)),
            Assigment::new(Field::from(5), Field::from(5), Field::from(25)),
            Assigment::new(Field::from(9), Field::from(16), Field::from(25)),
        ]);

        // random numbers (the b's)
        let rand = [
            f17(7),
            f17(4),
            f17(11),
            f17(12),
            f17(16),
            f17(2),
            f17(14),
            f17(11),
            f17(7),
        ];

        // values that are sent from the verifier to the prover
        let challange = Challange {
            alpha: f17(15),
            beta: f17(12),
            gamma: f17(13),
            z: f17(5),
            v: f17(12),
        };

        let proof = plonk.prove(&constrains, &assigments, &challange, rand);

        let expected = Proof {
            a_s: g1f(91, 66),
            b_s: g1f(26, 45),
            c_s: g1f(91, 35),
            z_s: g1f(32, 59),
            t_lo_s: g1f(12, 32),
            t_mid_s: g1f(26, 45),
            t_hi_s: g1f(91, 66),
            w_z_s: g1f(91, 35),
            w_z_omega_s: g1f(65, 98),
            a_z: f17(15),
            b_z: f17(13),
            c_z: f17(5),
            s_sigma_1_z: f17(1),
            s_sigma_2_z: f17(12),
            r_z: f17(15),
            z_omega_z: f17(15),
        };
        assert_eq!(proof, expected);

        let rand = [f17(4)];
        assert!(plonk.verify(&constrains, &proof, &challange, rand));
    }
}
