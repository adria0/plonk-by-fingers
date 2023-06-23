use crate::field::Field;
use crate::mulmodg::MulGroupMod;
use crate::poly::Poly;
use itertools::Itertools;
use num_bigint::BigUint;
use sha2::{Digest, Sha256};
use std::ops::Index;
use std::{convert::TryInto, marker::PhantomData};

// https://starkware.co/stark-101/

type FF = crate::field::U64Field<3221225473>;

fn sha256hex(data: String) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data.as_bytes());
    hex::encode(hasher.finalize())
}
fn as_neg_str(v: FF) -> String {
    let m1 = FF::zero() - FF::one();
    let limit = (m1 / 2.into()).unwrap() + FF::one();
    if v > limit {
        let mv = m1 - v + FF::one();
        format!("-{}", mv)
    } else {
        format!("{}", v)
    }
}

fn merkle_tree(data: &[FF]) -> String {
    fn recursive_build_tree(data: &[FF], node_id: usize) -> String {
        if node_id >= data.len() {
            // A leaf.
            let id_in_data = node_id - data.len();
            let leaf_data = as_neg_str(data[id_in_data]);
            let h = sha256hex(leaf_data);
            h
        } else {
            // An internal node.
            let left = recursive_build_tree(data, node_id * 2);
            let right = recursive_build_tree(data, node_id * 2 + 1);
            let h = sha256hex(left + &right);
            h
        }
    }
    recursive_build_tree(data, 1)
}

struct Channel {
    state: String,
}

impl Channel {
    pub fn new() -> Self {
        Self {
            state: String::from("0"),
        }
    }
    pub fn send(&mut self, s: &str) {
        self.state = sha256hex(format!("{}{}", self.state, s));
    }

    pub fn receive_random_int(&mut self, min: BigUint, max: BigUint) -> BigUint {
        let bytes = if self.state.len() % 2 == 0 {
            hex::decode(self.state.clone())
        } else {
            hex::decode(String::from("0") + &self.state)
        }
        .unwrap();
        let state = BigUint::from_bytes_be(&bytes);
        let num = min.clone() + state % (max - min + 1u32);
        self.state = sha256hex(self.state.clone());
        num
    }
    pub fn receive_random_field_element(&mut self) -> FF {
        let num = self.receive_random_int(0u64.into(), (FF::order() - 1).into());
        if let Some(n) = num.iter_u64_digits().next() {
            FF::from(n)
        } else {
            FF::zero()
        }
    }
}

#[test]
fn stark101_test_1() {
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

    let N = 1022;

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
    let eval_domain : Vec<_> = H.coset(w).iter().collect();

    assert_eq!(H.at(1), h);

    let w_inv = w.inv().unwrap();
    for i in 0..8192 {
        assert_eq!((w_inv * eval_domain[1]).pow(i as u64) * w, eval_domain[i]);
    }

    // Evaluate on a Coset
    // so `f_eval` this is the LDE (low degree externsion) of `f`
    let f_eval: Vec<_> = eval_domain.iter().map(|v| f.eval(&v)).collect();

    // Commitments --------------------------------------------

    let root = merkle_tree(&f_eval);
    assert_eq!(
        root,
        "6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04"
    );

    // Channel -------------------------------------------------------------

    let mut channel = Channel::new();

    channel.send(&root);

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

        let CP_eval: Vec<_> = eval_domain.iter().map(|v| CP.eval(&v)).collect();
        let CP_root = merkle_tree(&CP_eval);

        assert_eq!(
            CP_root,
            "a8c87ef9764af3fa005a1a2cf3ec8db50e754ccb655be7597ead15ed4a9110f1"
        );
    }

    let alpha0 = channel.receive_random_field_element();
    let alpha1 = channel.receive_random_field_element();
    let alpha2 = channel.receive_random_field_element();

    let cp = alpha0 * p0.clone() + alpha1 * p1.clone() + alpha2 * p2.clone();
    let cp_eval: Vec<_> = eval_domain.iter().map(|v| cp.eval(&v)).collect();

    // PART 3: FRI Commitments ==?==========================================
    // =====================================================================

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

    let fri_commit =
        |cp: Poly<_>, domain: Vec<_>, cp_eval: Vec<_>, channel: &mut Channel| -> (Vec<_>, Vec<_>, Vec<_>, Vec<_>) {
            let mut fri_polys = vec![cp];
            let mut fri_domains = vec![domain];
            let mut fri_merkles = vec![merkle_tree(&cp_eval)];
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
                fri_merkles.push(merkle_tree(&next_layer));
                fri_layers.push(next_layer);
                channel.send(&fri_merkles[fri_merkles.len() - 1]);
            }
            channel.send(&as_neg_str(
                *fri_polys[fri_merkles.len() - 1].get(0).unwrap(),
            ));
            (fri_polys, fri_domains, fri_layers, fri_merkles)
        };

    {
        // test
        let mut test_channel = Channel::new();
        let (fri_polys, fri_domains, fri_layers, fri_merkles) = fri_commit(cp, eval_domain, cp_eval,  &mut test_channel);

        assert_eq!(fri_layers.len() ,11);
        assert_eq!(fri_layers[fri_layers.len()-1].len(), 8);
        assert!(fri_layers[fri_layers.len()-1].iter().all(|v| v == &fri_layers[fri_layers.len()-1][0]));
        assert_eq!(fri_polys[fri_polys.len()-1].degree(),0);

        /*
        assert all([x == FieldElement(-1138734538) for x in fri_layers[-1]]), f'Expected last layer to be constant.'
        assert fri_merkles[-1].root == '1c033312a4df82248bda518b319479c22ea87bd6e15a150db400eeff653ee2ee', 'Last layer Merkle root is wrong.'
        assert test_channel.state == '61452c72d8f4279b86fa49e9fb0fdef0246b396a4230a2bfb24e2d5d6bf79c2e', 'The channel state is not as expected.'
        */
    }

    unreachable!();
}
