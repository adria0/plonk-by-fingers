use crate::ec::Field;
use crate::matrix::Matrix;
use crate::poly::Poly;

trait FFT<F: Field> {
    fn fft(&self, values: &[F]) -> Vec<F>;
    fn fft_inv(&self, freq: &[F]) -> Vec<F>;
}

struct VandermondeMatrix<F: Field>(Matrix<F>);

impl<F: Field> VandermondeMatrix<F> {
    pub fn new(omega: F, len: usize) -> Self {
        let mut values = Vec::with_capacity(len.pow(2));
        for n in 0..len {
            for m in 0..len {
                let v = omega.pow((n * m) as u64);
                values.push(v);
            }
        }
        Self(crate::matrix::Matrix::new(values, len, len))
    }
}

impl<F: Field> FFT<F> for VandermondeMatrix<F> {
    fn fft(&self, values: &[F]) -> Vec<F> {
        (&self.0 * Poly::new(values.to_vec())).into_coeffs()
    }
    fn fft_inv(&self, freq: &[F]) -> Vec<F> {
        let vals = self.fft(freq);
        let vals_len_modinv = F::from(freq.len() as u64).inv().unwrap();
        std::iter::once(&vals[0])
            .chain(vals.iter().rev().take(vals.len() - 1))
            .map(|v| vals_len_modinv * v)
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::U64Field;

    #[test]
    fn test_fft() {
        type F = U64Field::<337>;
        let domain = VandermondeMatrix::new(F::from(85u64), 8);
        let values = [3, 1, 4, 1, 5, 9, 2, 6].map(|x| F::from(x as u64));
        let current_freq = domain.fft(&values);

        let expected_freq = [31, 70, 109, 74, 334, 181, 232, 4].map(|x| F::from(x as u64));
        assert_eq!(&current_freq, &expected_freq[..]);

        let current_values = domain.fft_inv(&current_freq);
        assert_eq!(&current_values, &values[..]);
    }
    /*
    #[test]
    fn test_ffn_poly_mul() {
        type F = U64Field<337>;
        let zero = F::zero();
        let domain = gen_domain(F::from(85u64));

        let d = |l: &[F]| {
            l.iter()
                .map(|x| format!("{}", x))
                .collect::<Vec<String>>()
                .join(",")
        };

        let mut a_vals = [1, 1].map(|x| F::from(x as u64)).to_vec();
        let mut b_vals = [5, 11].map(|x| F::from(x as u64)).to_vec();
        let sum = a_vals.len() + b_vals.len();

        let mp1 = crate::poly::Poly::new((&a_vals[..]).to_vec())
            * crate::poly::Poly::new((&b_vals[..]).to_vec());

        a_vals.extend(vec![F::zero(); sum - a_vals.len()]);
        b_vals.extend(vec![F::zero(); sum - b_vals.len()]);

        println!("a_vals={:?}", d(&a_vals));
        println!("b_vals={:?}", d(&b_vals));

        println!("domain = {:?}", d(&domain));

        let a_freq = fft(&a_vals, &domain);
        let b_freq = fft(&b_vals, &domain);

        let mut c_freq = Vec::new();

        println!("a_freq={:?}", d(&a_freq));
        println!("b_freq={:?}", d(&b_freq));

        for n in 0..a_freq.len() {
            let l = a_freq.get(n).unwrap_or(&zero);
            let r = b_freq.get(n).unwrap_or(&zero);

            let mut carry = F::zero();
            c_freq.push(l.carrying_mul(r, &mut carry));
        }

        println!("c_freq={:?}", d(&c_freq));
        //assert!(carry.is_zero());

        let c_vals = fft_inv(&c_freq, &domain);

        let mp2 = crate::poly::Poly::new(c_vals);

        let values = [3, 1, 4, 1, 5, 9, 2, 6].map(|x| F::from(x as u64));
        println!(
            "fft_with_matrix {}",
            d(&fft_matrix(&values, F::from(85u64)))
        );
        println!("fft_with_split {}", d(&fft(&values, &domain)));
    }
    */
}
