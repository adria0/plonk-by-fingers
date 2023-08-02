mod g1;
mod g2;
mod gt;
mod pairing;
mod types;

#[cfg(test)]
mod tests {
    use crate::plonk::plonk_by_hand::{
        g1::g1f,
        types::{f101, f17, PlonkByHandTypes},
    };

    use super::super::super::plonk::*;

    #[test]
    fn test_plonk_gen_proof() {
        // create the trusted setup
        let s = f101(2); // the toxic waste
        let srs = Srs::<PlonkByHandTypes>::create(s, 6);

        let plonk = Plonk::new(
            srs,
            f17(4), // omega pows
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
