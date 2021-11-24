pub mod g1;
pub mod g2;
pub mod gt;
pub mod pairing;

use crate::{ec::Field, plonk::PlonkTypes, utils::U64Field};

pub type F101 = U64Field<101>;
pub fn f101(x: u64) -> F101 {
    F101::from(x)
}

pub type F17 = U64Field<17>;
pub fn f17(x: u64) -> F17 {
    F17::from(x)
}

#[derive(Debug, PartialEq)]
pub struct PlonkByHandTypes {}
impl PlonkTypes for PlonkByHandTypes {
    type G1 = g1::G1P;
    type G2 = g2::G2P;
    type GT = gt::GTP;
    type E = pairing::PBHPairing;
    type GF = F101;
    type SF = F17;

    fn gf(sg: F17) -> F101 {
        F101::from(sg.as_u64())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        constraints::{Assigment, Assigments, Constrains, CopyOf, Gate},
        pbh::g1::g1f,
        plonk::{Challange, Plonk, Proof, SRS},
    };

    #[test]
    fn test_plonk_gen_proof() {
        // create the trusted setup
        let s = f101(2); // the toxic waste
        let srs = SRS::<PlonkByHandTypes>::create(s, 6);

        let plonk = Plonk::new(
            srs,
            f17(4), // omega
            4,      // omega pows
            f17(2), // k1
            f17(3), // k1
        );

        // constraints and assigments
        let constraints = Constrains::new(
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
            Assigment::new(f17(3), f17(3), f17(9)),
            Assigment::new(f17(4), f17(4), f17(16)),
            Assigment::new(f17(5), f17(5), f17(25)),
            Assigment::new(f17(9), f17(16), f17(25)),
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

        let proof = plonk.prove(&constraints, &assigments, &challange, rand);

        let expected = Proof::<PlonkByHandTypes> {
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
        assert!(plonk.verify(&constraints, &proof, &challange, rand));
    }
}
