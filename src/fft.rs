use crate::ec::Field;
use crate::matrix::Matrix;
use crate::poly::Poly;

pub trait FFT<F: Field> {
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

// FE > = (a_vals.len + b_vals.len)  x ((F-1)^2)/2
// according https://stackoverflow.com/questions/52270320/implementing-fft-over-finite-fields
pub fn mul_ntt<F: Field, FE: Field, FFTI: FFT<FE>>(
    fft: FFTI,
    a_vals: Vec<F>,
    b_vals: Vec<F>,
    f_to_fe: &dyn Fn(F) -> FE,
    fe_to_f: &dyn Fn(FE) -> F,
) -> Vec<F> {
    let sum = a_vals.len() + b_vals.len();
    let zero = FE::zero();

    let mut a_vals: Vec<FE> = a_vals.into_iter().map(|v| f_to_fe(v)).collect();
    let mut b_vals: Vec<FE> = b_vals.into_iter().map(|v| f_to_fe(v)).collect();
    a_vals.extend(vec![FE::zero(); sum - a_vals.len()]);
    b_vals.extend(vec![FE::zero(); sum - b_vals.len()]);

    let a_freq = fft.fft(&a_vals);
    let b_freq = fft.fft(&b_vals);

    let mut c_freq = Vec::new();

    for n in 0..a_freq.len() {
        let l = a_freq.get(n).unwrap_or(&zero);
        let r = b_freq.get(n).unwrap_or(&zero);
        c_freq.push(*l * r);
    }

    fft.fft_inv(&c_freq).into_iter().map(|v| fe_to_f(v)).collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::U64Field;

    #[test]
    fn test_fft() {
        type F = U64Field<337>;
        let domain = VandermondeMatrix::new(F::from(85u64), 8);
        let values = [3, 1, 4, 1, 5, 9, 2, 6].map(|x| F::from(x as u64));
        let current_freq = domain.fft(&values);

        let expected_freq = [31, 70, 109, 74, 334, 181, 232, 4].map(|x| F::from(x as u64));
        assert_eq!(&current_freq, &expected_freq[..]);

        let current_values = domain.fft_inv(&current_freq);
        assert_eq!(&current_values, &values[..]);
    }

    #[test]
    fn test_ffn_poly_mul() {
        type F = U64Field<31>;
        type FE = U64Field<6737>;
        let f_to_fe = &|v: F| FE::from(v.as_u64());
        let fe_to_f = &|v: FE| F::from(v.as_u64() % 31);

        let a: Vec<F> = [24, 12, 28, 8].iter().map(|x| F::from(*x as u64)).collect();
        let b: Vec<F> = [4, 26, 29, 23].iter().map(|x| F::from(*x as u64)).collect();
        let poly_c = Poly::new(a.clone()) * Poly::new(b.clone());

        let domain = VandermondeMatrix::new(FE::from(5862u64), 8);
        let ntt_c  = Poly::new(mul_ntt(domain, a, b, &f_to_fe, &fe_to_f));

        assert_eq!(poly_c, ntt_c);
    }
}
