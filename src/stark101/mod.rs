#![allow(non_snake_case)]

use crate::poly::Poly;
use crate::field::Field;

mod channel;
mod example;
mod format;
mod mt;

// https://starkware.co/stark-101/

pub type FF = crate::field::U64Field<3221225473>;

struct FriProof<F: Field> {
    evals : Vec<(F,F)>
}
impl<F : Field> FriProof<F> {
    // Prove what we now the polinomial p
    fn prove<RNDS>(mut poly: Poly<F>, domain: &[F], rands: RNDS) -> Self
    where RNDS : IntoIterator<Item=F> {
        let mut rands = rands.into_iter();
        let idx = rands.next().unwrap().as_u64() as usize;
        let mut evals = Vec::new();
        let mut domain = domain.to_vec();

        while poly.degree() > 0 {
            let x = domain[idx % domain.len()];
            evals.push((poly.eval(&x),poly.eval(&-x)));

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

        Self { evals }

    }

    // Verify low degree
    fn verify<RNDS>(&self, domain: &[F], rands: RNDS) -> bool
    where RNDS : IntoIterator<Item=F> {
        let mut rands = rands.into_iter();
        let idx = rands.next().unwrap().as_u64() as usize;
        let mut domain = domain.to_vec();
        for layer in 0..self.evals.len() - 1 {
            let x = domain[idx % domain.len()];
            let rand = rands.next().unwrap();
            let check = ((rand + x) / (x * F::from(2))).unwrap() * self.evals[layer].0
                + ((rand - x) / (-x * F::from(2))).unwrap() * self.evals[layer].1;

            if self.evals[layer+1].0 != check  {
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
fn fri_base() {

}

#[test]
fn fri_hand() {
    use crate::field::Field;

    type F = crate::field::U64Field<41>;
    let w = F::from(27);

    let domain0 = (0..8).map(|n| w.pow(n)).collect::<Vec<_>>();
    assert_eq!(
        domain0,
        vec![1, 27, 32, 3, 40, 14, 9, 38]
            .into_iter()
            .map(F::from)
            .collect::<Vec<_>>()
    );

    // trace polynomial is 7+8x+9x^2
    // cp0 polinomial is the same
    let trace = Poly::<F>::from([7, 8, 9]); // 7+8x+9x^2

    let cp0 = Poly::lagrange(
        &domain0
            .iter()
            .map(|x| (x.clone(), trace.eval(x)))
            .collect::<Vec<_>>(),
    );

    // split poly in p(x) = g(x^2)+ xh(x^2)
    // g(x^2)=7+9x
    // h(x^2)=8

    let mut g_x2 = vec![];
    let mut h_x2 = vec![];

    for (i, coeff) in cp0.coeffs().iter().enumerate() {
        if i % 2 == 0 {
            g_x2.push(*coeff);
        } else {
            h_x2.push(*coeff);
        }
    }

    let g_x2 = Poly::new(g_x2);
    let h_x2 = Poly::new(h_x2);

    {
        let x2 = Poly::from([0, 0, 1]);
        let check = g_x2.compose(&x2) + Poly::x() * h_x2.compose(&x2);
        assert_eq!(check, cp0);
    }

    // use random beta, and random linear combine the two polinomials
    // p1 is in the domain x^2
    // betq = 7, y = x^2
    // cp1(y) = g(y)+7·h(y)
    // cp1(x) = 7+9x+7·9 = 9x + 22
    // new domain is [1,32,40,9]

    let beta1 = F::from(11);
    let cp1 = g_x2.clone() + h_x2.clone() * Poly::new([beta1]);

    let domain1 = (0..domain0.len() / 2)
        .map(|i| domain0[i * 2])
        .collect::<Vec<_>>();

    assert_eq!(
        domain1,
        vec![1, 32, 40, 9]
            .into_iter()
            .map(F::from)
            .collect::<Vec<_>>()
    );

    // generate proof, index comes from transript, gives idx_{pos,neg}_cp{0,1}
    // ---------------------------------------------------------------------------------

    let idx = 5;

    // x domain: 1,27(+),32,3,40,14,9,38(-)
    let idx_pos_cp0 = cp0.eval(&domain0[idx % domain0.len()]);
    let idx_neg_cp0 = cp0.eval(&-domain0[idx % domain0.len()]);

    // x^2 domain: 1 32(+) 40 9(-)
    let idx_pos_cp1 = cp1.eval(&domain1[idx % domain1.len()]);
    let idx_neg_cp1 = cp1.eval(&-domain1[idx % domain1.len()]);

    // verify transition from cp0 to cp1

    let x = domain0[idx % domain0.len()];
    let check = ((beta1 + x) / (x * F::from(2))).unwrap() * idx_pos_cp0
        + ((beta1 - x) / (-x * F::from(2))).unwrap() * idx_neg_cp0;

    assert_eq!(check, idx_pos_cp1);
    let rands = [2,6,7,3,1,3,4].map(F::from);

    let fri = Fri::prove(Poly::from([1,2,3,4,5,6,8]), &domain0, rands);
    assert!(fri.verify(&domain0, rands));

}

#[test]
fn stark101_test_1() {
    use itertools::Itertools;

    use crate::{
        field::Field,
        mulmodg::MulGroupMod,
        poly::Poly,
        stark101::{
            channel::Channel,
            format::{as_neg_str, sha256hex},
            mt::MerkleTree,
        },
    };

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

    // PART I : Trace and Low-Degree Extension =============================
    // =====================================================================

    // FibonacciSq Trace ------------------------------

    let fibonacci_sq = |x: FF, n: usize| -> Vec<FF> {
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
    };

    let a = fibonacci_sq(3141592.into(), 1022);
    assert_eq!(a[1022], 2338775057u64.into());

    // Thinking of Polynomials ------------------------------

    let w: FF = 5.into();

    // G mul mod, g is 5 ^ (3 * 2^20)

    let g = w.pow(3 * 2u64.pow(20));
    let mm = MulGroupMod::<FF>::new(g);
    let G: Vec<_> = mm.iter().take(1024).collect();

    // check ciclic group

    assert_eq!(G[1023] * g, 1.into());

    // create interpolation poly

    let f = mm.lagrange(a);
    assert_eq!(f.eval(&2.into()), 1302089273.into());

    // Evaluating on a Larger Domain --------------------------

    // We say that f is a low degree polinomial (1024 points)
    // Creating a larger domain polinomial ( Low Degeee Extension ) with the same
    // points over a Coset ( with 8k points )

    // Cosets

    let h = w.pow((3 * 2u64.pow(30)) / 8192);
    let H = MulGroupMod::<FF>::new(h);
    let eval_domain: Vec<_> = H.coset(w).iter().collect();

    assert_eq!(H.at(1), h);

    let w_inv = w.inv().unwrap();
    for i in 0..8192 {
        assert_eq!((w_inv * eval_domain[1]).pow(i as u64) * w, eval_domain[i]);
    }

    // Evaluate on a Coset
    // so `f_eval` this is the LDE (low degree externsion) of `f`
    let f_eval: Vec<_> = eval_domain.iter().map(|v| f.eval(v)).collect();

    // Commitments --------------------------------------------

    let mt = MerkleTree::new(&f_eval);
    assert_eq!(
        mt.root(),
        "6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04"
    );

    // Channel -------------------------------------------------------------

    let mut channel = Channel::new();

    channel.send(&mt.root());

    // PART II : Constraints ===============================================
    // =====================================================================

    // Step 1 - FibonacciSq Constraints
    // Step 2 - Polynomial Constraints
    // Step 3 - Rational Functions (That are in Fact Polynomials)

    // First constraint, first point exists in f, so we remove this point without remainder
    let (x0, y0) = (G[0], f.eval(&G[0]));
    let p0 = (f.clone() - y0) / (Poly::x() - x0);

    // First constraint, last point exists in f, so we remove this point without remainder
    let (x1, y1) = (G[1022], f.eval(&G[1022]));
    let p1 = (f.clone() - y1) / (Poly::x() - x1);

    // Third constraint
    let numer2 = f.compose(&(Poly::x() * g.pow(2))) - f.compose(&(Poly::x() * g)).pow(2) - f.pow(2);
    let denom2_numer: Poly<FF> = Poly::mset([(0, FF::from(-1)), (1024, FF::one())]);
    let denom2_denom = Poly::new([-g.pow(1021), FF::one()])
        * Poly::new([-g.pow(1022), FF::one()])
        * Poly::new([-g.pow(1023), FF::one()]);
    let p2 = numer2 / (denom2_numer / denom2_denom);

    assert_eq!(p2.degree(), 1023);
    assert_eq!(p2.eval(&FF::from(31415)), 2090051528.into());

    // Step 4 - Composition Polynomial

    // Commit on the Composition Polynomial
    {
        let mut test_channel = Channel::new();
        let alpha0 = test_channel.receive_random_field_element();
        let alpha1 = test_channel.receive_random_field_element();
        let alpha2 = test_channel.receive_random_field_element();
        let CP = alpha0 * p0.clone() + alpha1 * p1.clone() + alpha2 * p2.clone();

        assert_eq!(CP.degree(), 1023);
        assert_eq!(CP.eval(&FF::from(2439804)), 838767343.into());

        let CP_eval: Vec<_> = eval_domain.iter().map(|v| CP.eval(v)).collect();
        let CP_root = MerkleTree::new(&CP_eval).root();

        assert_eq!(
            CP_root,
            "a8c87ef9764af3fa005a1a2cf3ec8db50e754ccb655be7597ead15ed4a9110f1"
        );
    }

    let alpha0 = channel.receive_random_field_element();
    let alpha1 = channel.receive_random_field_element();
    let alpha2 = channel.receive_random_field_element();

    let cp = alpha0 * p0 + alpha1 * p1 + alpha2 * p2;
    let cp_eval: Vec<_> = eval_domain.iter().map(|v| cp.eval(v)).collect();

    // PART 3: FRI Commitments ============================================
    // ====================================================================

    // Domain Generation

    let next_fri_domain = |fri_domain: &[FF]| {
        fri_domain
            .iter()
            .map(|v| v * v)
            .take(fri_domain.len() / 2)
            .collect::<Vec<_>>()
    };

    {
        // test
        let s = sha256hex(
            next_fri_domain(&eval_domain)
                .iter()
                .map(|v| as_neg_str(*v))
                .join(","),
        );
        assert_eq!(
            s,
            "5446c90d6ed23ea961513d4ae38fc6585f6614a3d392cb087e837754bfd32797"
        );
    }

    let next_fri_polynomial = |poly: &Poly<FF>, beta: FF| -> Poly<FF> {
        let mut even_coeffs = vec![];
        let mut odd_coeffs = vec![];

        for (i, coeff) in poly.coeffs().iter().enumerate() {
            if i % 2 == 0 {
                even_coeffs.push(*coeff);
            } else {
                odd_coeffs.push(*coeff);
            }
        }

        Poly::new(odd_coeffs) * beta + Poly::new(even_coeffs)
    };

    {
        // vector test changed due bug, see https://github.com/starkware-industries/stark101/issues/8
        let next_p = next_fri_polynomial(&cp, FF::from(987654321));
        let s = sha256hex(next_p.coeffs().iter().map(|v| as_neg_str(*v)).join(","));
        assert_eq!(
            "242f36b1d7d5b3e19948e892459774f14c038bc5864ba8884817112aa1405257",
            s
        );
    }

    let next_fri_layer =
        |poly: &Poly<FF>, domain: &[FF], beta: FF| -> (Poly<FF>, Vec<FF>, Vec<FF>) {
            let next_poly = next_fri_polynomial(poly, beta);
            let next_domain = next_fri_domain(domain);
            let next_layer: Vec<_> = next_domain.iter().map(|v| next_poly.eval(v)).collect();
            (next_poly, next_domain, next_layer)
        };

    {
        // test
        let test_poly = Poly::from([2, 3, 0, 1]);
        let test_domain = vec![FF::from(3), FF::from(5)];
        let beta = FF::from(7);
        let (next_p, next_d, next_l) = next_fri_layer(&test_poly, &test_domain, beta);
        assert_eq!(next_p, Poly::from([23, 7]));
        assert_eq!(next_d, vec![FF::from(9)]);
        assert_eq!(next_l, vec![FF::from(86)]);
    }

    let fri_commit = |cp: Poly<_>,
                      domain: Vec<_>,
                      cp_eval: Vec<_>,
                      channel: &mut Channel|
     -> (Vec<_>, Vec<_>, Vec<_>, Vec<_>) {
        let mut fri_polys = vec![cp];
        let mut fri_domains = vec![domain];
        let mut fri_merkles = vec![MerkleTree::new(&cp_eval)];
        let mut fri_layers = vec![cp_eval];

        while fri_polys[fri_polys.len() - 1].degree() > 0 {
            let beta = channel.receive_random_field_element();
            let (next_poly, next_domain, next_layer) = next_fri_layer(
                &fri_polys[fri_polys.len() - 1],
                &fri_domains[fri_domains.len() - 1],
                beta,
            );

            fri_polys.push(next_poly);
            fri_domains.push(next_domain);
            fri_merkles.push(MerkleTree::new(&next_layer));
            fri_layers.push(next_layer);
            channel.send(&fri_merkles[fri_merkles.len() - 1].root());
        }

        let n = fri_polys.iter_mut().last().and_then(|v| v.get(0)).unwrap();
        channel.send_ff(*n);
        (fri_polys, fri_domains, fri_layers, fri_merkles)
    };

    let (fri_polys, fri_domains, fri_layers, fri_merkles) =
        fri_commit(cp, eval_domain, cp_eval, &mut channel);

    assert_eq!(fri_layers.len(), 11, "Expected number of FRI layers is 11");
    assert_eq!(
        fri_layers[fri_layers.len() - 1].len(),
        8,
        "Expected last layer to contain exactly 8 elements"
    );

    assert!(
        fri_layers[fri_layers.len() - 1]
            .iter()
            .all(|v| v == &FF::from(1443223587)),
        "Expected last layer to be constant"
    );

    assert_eq!(
        fri_polys[fri_polys.len() - 1].degree(),
        0,
        "Expected last polynomial to be constant (degree 0)"
    );
    assert_eq!(
        fri_merkles[fri_merkles.len() - 1].root(),
        "38fe69bdc218d0327beba3197aaee3c0ba707912d0dec81301b8a7fca6a022bf"
    );

    assert_eq!(
        channel.state,
        "eba2039aca50b08d50c5f4775863e1adccc9d133dd52e74b624efb9937780ae7"
    );

    // PART 4: FRI Commitments ============================================
    // ====================================================================

    let decommit_on_fri_layers = |idx: usize, channel: &mut Channel| {
        for (layer, merkle) in fri_layers
            .iter()
            .dropping_back(1)
            .zip_eq(fri_merkles.iter().dropping_back(1))
        {
            let length = layer.len();
            let idx = idx % length;
            let sib_idx = (idx + length / 2) % length;
            channel.send_ff(layer[idx]);
            channel.send_path(&merkle.get_authentication_path(idx));
            channel.send_ff(layer[sib_idx]);
            channel.send_path(&merkle.get_authentication_path(sib_idx));
        }
        channel.send_ff(*fri_layers.iter().last().and_then(|l| l.get(0)).unwrap());
    };

    {
        // test
        let mut test_channel = Channel::new();
        for query in [7527, 8168, 1190, 2668, 1262, 1889, 3828, 5798, 396, 2518] {
            decommit_on_fri_layers(query, &mut test_channel);
        }
        assert_eq!(
            test_channel.state,
            "0c7931382b6a846b4d91b485dcb51f8511b971f94513f72870392bfe7641ca36"
        );
    }

    let decommit_on_query = |idx: usize, channel: &mut Channel| {
        assert!(idx + 16 < f_eval.len());
        channel.send_ff(f_eval[idx]); // f(x).
        channel.send_path(&mt.get_authentication_path(idx)); // auth path for f(x).
        channel.send_ff(f_eval[idx + 8]); // f(gx).
        channel.send_path(&mt.get_authentication_path(idx + 8)); // auth path for f(gx).
        channel.send_ff(f_eval[idx + 16]); // f(g^2x).
        channel.send_path(&mt.get_authentication_path(idx + 16)); // auth path for f(g^2x).
        decommit_on_fri_layers(idx, channel)
    };

    {
        // test
        let mut test_channel = Channel::new();
        for query in [8134, 1110, 1134, 6106, 7149, 4796, 144, 4738, 957] {
            decommit_on_query(query, &mut test_channel);
        }

        assert_eq!(
            test_channel.state,
            "856682d2782140a10a371eb258ee85c05618f34ee0c695ae836aec6e4fc3f9b1"
        );
    }

    let decommit_fri = |channel: &mut Channel| {
        for _ in 0..3 {
            // Get a random index from the verifier and send the corresponding decommitment.
            let rnd = channel.receive_random_int(0u64.into(), (8191u64 - 16).into());
            let idx = if let Some(n) = rnd.to_u64_digits().first() {
                *n as usize
            } else {
                0
            };
            decommit_on_query(idx, channel);
        }
    };

    {
        let mut test_channel = Channel::new();
        decommit_fri(&mut test_channel);
        assert_eq!(
            test_channel.state,
            "3424b78347d24261f921cd1a3b6dec864c90729b8d7c2c678cde327c68ab7c5b"
        );
    }
}
