use std::collections::HashMap;

use crate::field::Field;

use super::format::*;
use super::FF;

pub struct MerkleTree {
    data: Vec<FF>,
    height: u32,
    root: String,
    facts: HashMap<String, (String, String)>,
}

impl MerkleTree {
    fn bin(mut n: usize) -> Vec<bool> {
        let mut r = Vec::new();
        while n > 0 {
            r.push(n % 2 == 1);
            n = n / 2;
        }
        r.reverse();
        r
    }

    fn recursive_build_tree(
        facts: &mut HashMap<String, (String, String)>,
        data: &[FF],
        node_id: usize,
    ) -> String {
        if node_id >= data.len() {
            // A leaf.
            let id_in_data = node_id - data.len();
            let leaf_data = as_neg_str(data[id_in_data]);
            let h = sha256hex(leaf_data);
            h
        } else {
            // An internal node.
            let left = Self::recursive_build_tree(facts, data, node_id * 2);
            let right = Self::recursive_build_tree(facts, data, node_id * 2 + 1);
            let h = sha256hex(format!("{left:}{right:}"));
            facts.insert(h.clone(), (left, right));
            h
        }
    }

    pub fn new(data: &[FF]) -> Self {
        let num_leaves = 2_usize.pow((data.len() as f32).log2().ceil() as u32);
        let mut data = data.to_vec();
        data.append(&mut vec![FF::zero(); num_leaves - data.len()]);
        let height = num_leaves.checked_ilog2().unwrap();
        let mut facts = HashMap::new();
        let root = Self::recursive_build_tree(&mut facts, &data, 1);

        Self {
            facts,
            data,
            height,
            root,
        }
    }

    pub fn get_authentication_path(&self, leaf_id: usize) -> Vec<String> {
        assert!(leaf_id < self.data.len());
        let node_id = leaf_id + self.data.len();
        let mut cur: String = self.root.clone();

        let mut decommitment = vec![];

        // In a Merkle Tree, the path from the root to a leaf, corresponds to the the leaf id's
        // binary representation, starting from the second-MSB, where '0' means 'left', and '1' means
        // 'right'.
        // We therefore iterate over the bits of the binary representation - skipping the '0b'
        // prefix, as well as the MSB.

        let path = Self::bin(node_id).into_iter().skip(1);
        for bit in path {
            let (mut left, mut right) = self.facts.get(&cur).unwrap().clone();
            if bit {
                std::mem::swap(&mut left, &mut right);
            }
            decommitment.push(right);
            cur = left;
        }

        decommitment
    }

    pub fn verify_decommitment(
        leaf_id: usize,
        leaf_data: FF,
        decommitment: &[String],
        root: &str,
    ) -> bool {
        let leaf_num = 2usize.pow(decommitment.len() as u32);
        let node_id = leaf_id + leaf_num;
        let mut cur = sha256hex(as_neg_str(leaf_data));

        let path = Self::bin(node_id).into_iter().skip(1);

        for (bit, auth) in path.zip(decommitment.iter()).rev() {
            let h = if bit {
                format!("{auth:}{cur:}")
            } else {
                format!("{cur:}{auth:}")
            };
            cur = sha256hex(h);
        }
        cur == root
    }
    pub fn root(&self) -> String {
        self.root.clone()
    }
}

#[test]
fn mt_witness() {
    let data = [1, 2, 3, 4, 7, 9, 11].map(FF::from);
    let mt = MerkleTree::new(&data);
    for (leaf_id, _) in data.iter().enumerate() {
        let decommitment = mt.get_authentication_path(leaf_id);
        assert!(MerkleTree::verify_decommitment(
            leaf_id,
            data[leaf_id],
            &decommitment,
            &mt.root()
        ))
    }
}
