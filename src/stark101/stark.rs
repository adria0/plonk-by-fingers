#![allow(non_snake_case, dead_code)]
#![allow(clippy::needless_range_loop)]

use itertools::Itertools;

use crate::{
    field::Field,
    mulmodg::MulGroupMod,
    poly::Poly,
    stark101::{
        FF,
        channel::Channel,
        mt::MerkleTree,
    },
};
// https://starkware.co/stark-101/

pub trait Program {
    fn witness(public: Vec<FF>) -> Vec<FF>;
    fn constrains(G: &[FF], f: &Poly<FF>) -> Vec<Poly<FF>>;
    fn rotations() -> Vec<usize>;
}


pub struct StarkDomain {
    pub domain: Vec<FF>,
    pub coset: Vec<FF>,
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

pub struct FriCommit {
    pub cps_at_domains: Vec<Vec<FF>>,
    pub cps_at_domains_root: Vec<String>,
    pub cps_at_domains_leaf: FF,
    pub cps_betas: Vec<FF>,
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

pub struct FriQuery {
    idx: usize,
    cp: CpDecommit,
    witness: Vec<WitnessDecomit>,
}

impl FriQuery {
    pub fn verify(
        &self,
        domain: &[FF],
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

        // Check FRI layers
        let next_domain = |domain: &[FF]| {
            domain
                .iter()
                .map(|v| v * v)
                .take(domain.len() / 2)
                .collect::<Vec<_>>()
        };
        let mut domain = domain.to_vec();
        for layer in 0..self.cp.layers.len() - 1 {
            let rand = betas[layer];
            let x = domain[self.idx % domain.len()];

            let check = ((rand + x) / (x * FF::from(2))).unwrap() * self.cp.layers[layer].pos_x
            + ((rand - x) / (-x * FF::from(2))).unwrap() * self.cp.layers[layer].neg_x;

            assert_eq!(check, self.cp.layers[layer + 1].pos_x);
            domain = next_domain(&domain);
        }

        // Check MTs

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

    pub fn decommit(
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

        for _ in 0..3 {
            // Get a random index from the verifier and send the corresponding decommitment.
            let rnd = channel.receive_random_int(0u64.into(), (8191u64 - 16).into());
            let idx = if let Some(n) = rnd.to_u64_digits().first() {
                *n as usize
            } else {
                1
            };

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

