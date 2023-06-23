use num_bigint::BigUint;

use crate::field::Field;

use super::format::*;
use super::FF;

pub struct Channel {
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
