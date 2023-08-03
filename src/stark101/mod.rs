#![allow(non_snake_case)]
#![allow(clippy::module_inception)]

use crate::field::{Field, U64Field};
use crate::poly::Poly;

mod channel;
mod stark;
mod format;
mod mt;
mod fri;

// https://starkware.co/stark-101/

pub type FF = crate::field::U64Field<3221225473>;

#[cfg(test)]
mod stark101 {

    use itertools::Itertools;

    use crate::{stark101::stark::Program, field::Field, poly::Poly, stark101::stark::StarkDomain};
    use super::{FF, mt::MerkleTree, channel::Channel, stark::FriCommit};

    struct FibonacciSq {}

    impl Program for FibonacciSq {
        fn witness(public: Vec<FF>) -> Vec<FF> {
            let n = 1022;

            let x = public[0];

            let mut out = Vec::new();
            let mut a0 = FF::one();
            let mut a1 = x;

            out.push(a0);
            out.push(a1);

            for _ in 0..n - 1 {
                let a3 = a0 * a0 + a1 * a1;
                out.push(a3);

                a0 = a1;
                a1 = a3;
            }

            out
        }
        fn constrains(domain: &[FF], witness_on_domain: &Poly<FF>) -> Vec<Poly<FF>> {
            // First constraint, first point exists in f, so we remove this point without remainder
            let (x0, y0) = (domain[0], witness_on_domain.eval(&domain[0]));
            let p0 = (witness_on_domain.clone() - y0) / (Poly::x() - x0);

            // First constraint, last point exists in f, so we remove this point without remainder
            let (x1, y1) = (domain[1022], witness_on_domain.eval(&domain[1022]));
            let p1 = (witness_on_domain.clone() - y1) / (Poly::x() - x1);

            // Third constraint
            let numer2 = witness_on_domain.compose(&(Poly::x() * domain[2]))
                - witness_on_domain.compose(&(Poly::x() * domain[1])).pow(2)
                - witness_on_domain.pow(2);
            let denom2_numer: Poly<FF> = Poly::mset([(0, FF::from(-1)), (1024, FF::one())]);
            let denom2_denom = Poly::new([-domain[1021], FF::one()])
                * Poly::new([-domain[1022], FF::one()])
                * Poly::new([-domain[1023], FF::one()]);
            let p2 = numer2 / (denom2_numer / denom2_denom);

            vec![p0, p1, p2]
        }
        fn rotations() -> Vec<usize> {
            vec![
                0,  // f(x)
                8,  // f(gx)
                16, // f(g^2x)
            ]
        }
    }


#[test]
fn stark101() {
    /*
        Overview

        0)  Program is a(x) = a(x-1)^2 + a(x-2)^2, a(0)=1, a(1022)=2338775057
        1)  Compute the trace of the execution (witness) 1024 points and we
            interpolate the results `f`=interpolate(g^i,witness(i))
        2)  Create the Low Degree Extension with 8192 points (`eval_domain`=wg^i and
            interpolate f in this domain `f_eval`=interpolate(eval_domain[i],f(eval_domain[i]))
            note: it's a coset
        3)  Generate three polinomial constraints
            - p0(x) = (f(x)-1) / (x-g^0)
            - p1(x) = (f(x)-2338775057) / (x-g^1022)
            - p2(x) = (f(xg^2) - f(gx)^2 - f(x)^2)) / ( (x-g0)*...*(x-g^20) )
        4)  Create Composition polinomial using random linear combination
            - `CP` =  α0 p0(x) + α1 p1(x) + α2 p2(x)
        5)  Trick: proving CP is a polynomial <=> CP is close to a polinomal of low degree
        6)  FRI: in each iteration reduce the degree of polinomial CP by splitting in odd
            and even coeffs with the help of a binding value from the transcript.
        7)  The prover sends to the verifier
            - The MT root of the execution trace and the MT path of the public witness
            - The MT root of the composition polinomial
            -
    */

    let x = 3141592.into();
    let witness = FibonacciSq::witness(vec![x]);

    let dom = StarkDomain::new();
    let witness_at_domain = dom.lagrange(witness);
    let witness_at_coset: Vec<_> = dom.eval_at_coset(&witness_at_domain);

    // constraints
    let constraints = FibonacciSq::constrains(&dom.domain, &witness_at_domain);

    // composition polynomial
    let witness_at_coset_mt = MerkleTree::new(&witness_at_coset);
    let mut channel = Channel::new();
    channel.send(&witness_at_coset_mt.root());

    let mut cp = Poly::zero();
    let alphas: Vec<FF> = constraints
        .iter()
        .map(|_| channel.receive_random_field_element())
        .collect();

    for (constrain, alpha) in constraints.iter().zip_eq(alphas.iter()) {
        cp = cp + constrain * alpha;
    }

    // Prove that cp is close to a polynomial of low degree

    let commit = FriCommit::commit(cp, dom.coset.clone(), &mut channel);
    let queries = commit.decommit(
        &witness_at_coset_mt,
        witness_at_coset,
        FibonacciSq::rotations(),
        &mut channel,
    );

    for q in queries {
        assert!(q.verify(
            &dom.coset,
            witness_at_coset_mt.root(),
            &commit.cps_at_domains_root,
            FibonacciSq::rotations(),
            &alphas,
            &commit.cps_betas,
            &constraints
        ));
    }
}

}