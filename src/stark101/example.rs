#![allow(non_snake_case, dead_code)]

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
// https://starkware.co/stark-101/

pub type FF = crate::field::U64Field<3221225473>;

trait Program {
    fn witness(public: Vec<FF>) -> Vec<FF>;
    fn constrains(G: &[FF], f: &Poly<FF>) -> Vec<Poly<FF>>;
    fn rotations() -> Vec<usize>;
}

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

struct StarkDomain {
    domain: Vec<FF>,
    coset: Vec<FF>,
}

impl StarkDomain {
    pub fn new() -> StarkDomain {
        // initial group
        let w: FF = 5.into();
        let g = w.pow(3 * 2u64.pow(20));
        let initial_group = MulGroupMod::<FF>::new(g);
        let initial: Vec<_> = initial_group.iter().take(1024).collect();

        // coset
        let h = w.pow((3 * 2u64.pow(30)) / 8192);
        let larger_group = MulGroupMod::<FF>::new(h).coset(w);
        let larger: Vec<_> = larger_group.iter().collect();

        Self {
            domain: initial,
            coset: larger,
        }
    }
    pub fn lagrange<I: IntoIterator<Item = FF>>(&self, values: I) -> Poly<FF> {
        let y = values.into_iter();
        let points: Vec<(FF, FF)> = self.domain.clone().into_iter().zip(y).collect::<Vec<_>>();
        Poly::lagrange(&points)
    }

    pub fn eval_at_coset(&self, p: &Poly<FF>) -> Vec<FF> {
        self.coset.iter().map(|v| p.eval(v)).collect()
    }
}

struct FriCommit {
    cps_at_domains: Vec<Vec<FF>>,
    cps_at_domains_root: Vec<String>,
    cps_at_domains_leaf: FF,
    cps_betas: Vec<FF>,
}

struct FriLayer {
    length: usize,
    root: String,
    pos_x: FF,
    neg_x: FF,
    pos_x_path: Vec<String>,
    neg_x_path: Vec<String>,
}

struct CpDecommit {
    layers: Vec<FriLayer>,
    leaf: FF,
}

struct WitnessDecomit {
    path: Vec<String>,
    leaf: FF,
}

struct FriQuery {
    idx: usize,
    cp: CpDecommit,
    witness: Vec<WitnessDecomit>,
}

impl FriQuery {
    fn verify(
        &self,
        witness_root: String,
        cp_roots: &[String],
        rotations: Vec<usize>,
        alphas: &[FF],
        betas: &[FF],
        constraints: &[Poly<FF>],
    ) -> bool {
        // Check that witness proofs are ok
        for (wd, rot) in self.witness.iter().zip_eq(rotations.iter()) {
            if !MerkleTree::verify_decommitment(self.idx + rot, wd.leaf, &wd.path, &witness_root) {
                return false;
            }
        }

        // Recompute CP from constraints
        let mut cp = Poly::zero();
        for (constrain, alpha) in constraints.iter().zip_eq(alphas.iter()) {
            cp = cp + constrain * alpha;
        }

        println!("verify::cp={}", cp.hash());

        let l = 0;
        let eval_g_x2 =
            ((self.cp.layers[l].pos_x - self.cp.layers[l].neg_x) / FF::from(2)).unwrap();
        let eval_h_x2 = ((self.cp.layers[l].pos_x - self.cp.layers[l].neg_x)
            / (FF::from(2) * FF::from(self.idx as u64)))
        .unwrap();
        let eval_cp_next = eval_g_x2 + betas[l] * eval_h_x2;

        assert_eq!(eval_cp_next, self.cp.layers[l + 1].pos_x);

        for (root, layer) in cp_roots
            .iter()
            .dropping_back(1)
            .zip_eq(self.cp.layers.iter().dropping(1))
        {
            let idx = std::cmp::max(self.idx % layer.length, 1);
            let sib_idx = (idx + layer.length / 2) % layer.length;

            if !MerkleTree::verify_decommitment(idx, layer.pos_x, &layer.pos_x_path, root) {
                return false;
            }
            if !MerkleTree::verify_decommitment(sib_idx, layer.neg_x, &layer.neg_x_path, root) {
                return false;
            }
        }
        true
    }
}

struct FriDecommit(Vec<FriQuery>);

impl FriCommit {
    pub fn commit(cp: Poly<FF>, domain: Vec<FF>, channel: &mut Channel) -> FriCommit {
        println!("commit::cp={}", cp.hash());

        let next_domain = |domain: &[FF]| {
            domain
                .iter()
                .map(|v| v * v)
                .take(domain.len() / 2)
                .collect::<Vec<_>>()
        };

        let next_cp = |poly: &Poly<FF>, beta: FF| -> Poly<FF> {
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

        let poly_at_domain = domain.iter().map(|v| cp.eval(v)).collect::<Vec<_>>();

        let mut cps = vec![cp];
        let mut domains = vec![domain];
        let mut cps_at_domains = vec![poly_at_domain];

        let mut roots = vec![];
        let mut betas = vec![];

        while cps[cps.len() - 1].degree() > 0 {
            let beta = channel.receive_random_field_element();

            let next_cp = next_cp(&cps[cps.len() - 1], beta);
            let next_domain = next_domain(&domains[domains.len() - 1]);
            let next_cp_at_domain: Vec<FF> = next_domain.iter().map(|v| next_cp.eval(v)).collect();

            println!("commit::cp_next={}", next_cp.hash());
            cps.push(next_cp);

            domains.push(next_domain);
            let root = MerkleTree::new(&next_cp_at_domain).root();
            channel.send(&root);
            roots.push(root);
            betas.push(beta);
            cps_at_domains.push(next_cp_at_domain);
        }

        let n = cps.iter_mut().last().and_then(|v| v.get(0)).unwrap();
        channel.send_ff(*n);

        FriCommit {
            cps_at_domains,
            cps_at_domains_root: roots,
            cps_at_domains_leaf: *n,
            cps_betas: betas,
        }
    }

    fn decommit(
        &self,
        mt: &MerkleTree,
        witness_at_coset: Vec<FF>,
        rotations: Vec<usize>,
        channel: &mut Channel,
    ) -> Vec<FriQuery> {
        let decommit_on_layers = |idx: usize, channel: &mut Channel| -> CpDecommit {
            let mut layers = vec![];

            let nodes = &self.cps_at_domains[..self.cps_at_domains.len() - 1];
            let leaf = &self.cps_at_domains[self.cps_at_domains.len() - 1];

            for poly_at_domain in nodes {
                let mt = MerkleTree::new(poly_at_domain);

                let length = poly_at_domain.len();
                let idx = idx % length;
                let sib_idx = (idx + length / 2) % length;

                println!("sampling at index = {}", idx);

                let layer = FriLayer {
                    length,
                    root: mt.root(),
                    pos_x: poly_at_domain[idx],
                    pos_x_path: mt.get_authentication_path(idx),
                    neg_x: poly_at_domain[sib_idx],
                    neg_x_path: mt.get_authentication_path(sib_idx),
                };

                channel.send_ff(layer.pos_x);
                channel.send_path(&layer.pos_x_path);
                channel.send_ff(layer.neg_x);
                channel.send_path(&layer.neg_x_path);

                layers.push(layer);
            }
            // send leaf
            channel.send_ff(leaf[0]);
            CpDecommit {
                layers,
                leaf: leaf[0],
            }
        };

        let mut queries = vec![];

        for _ in 0..1 {
            // Get a random index from the verifier and send the corresponding decommitment.
            let rnd = channel.receive_random_int(0u64.into(), (8191u64 - 16).into());
            let mut idx = if let Some(n) = rnd.to_u64_digits().first() {
                *n as usize
            } else {
                0
            };
            let idx = 2;
            // = std::cmp::max(idx,2);

            let mut witness_decommit = vec![];

            for rotation in &rotations {
                assert!(idx + rotation < witness_at_coset.len());

                let leaf = witness_at_coset[idx + rotation];
                let path = mt.get_authentication_path(idx + rotation);

                channel.send_ff(leaf); // f(x).
                channel.send_path(&path); // auth path for f(x).

                witness_decommit.push(WitnessDecomit { leaf, path });
            }
            let cp_decommit = decommit_on_layers(idx, channel);

            queries.push(FriQuery {
                idx,
                witness: witness_decommit,
                cp: cp_decommit,
            });
        }

        queries
    }
}

#[test]
fn stark101_test_2() {
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

    let commit = FriCommit::commit(cp, dom.coset, &mut channel);
    let queries = commit.decommit(
        &witness_at_coset_mt,
        witness_at_coset,
        FibonacciSq::rotations(),
        &mut channel,
    );

    for q in queries {
        assert!(q.verify(
            witness_at_coset_mt.root(),
            &commit.cps_at_domains_root,
            FibonacciSq::rotations(),
            &alphas,
            &commit.cps_betas,
            &constraints
        ));
    }
}
