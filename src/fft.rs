use crate::ec::Field;
use crate::matrix::Matrix;
use crate::poly::Poly;

pub struct EvaluationDomainGenerator<F: Field> {
    omega: F,
    size: usize,
}

impl<F: Field> EvaluationDomainGenerator<F> {
    pub fn new(omega: F, size: usize) -> Self {
        EvaluationDomainGenerator { omega, size }
    }
}

pub trait FFT<F: Field> {
    fn new(domain: EvaluationDomainGenerator<F>) -> Self;
    fn fft(&self, values: &[F]) -> Vec<F>;
    fn fft_inv(&self, freq: &[F]) -> Vec<F>;
}

struct VandermondeMatrix<F: Field>(Matrix<F>);

impl<F: Field> VandermondeMatrix<F> {}

impl<F: Field> FFT<F> for VandermondeMatrix<F> {
    fn new(domain: EvaluationDomainGenerator<F>) -> Self {
        let mut values = Vec::with_capacity(domain.size.pow(2));
        for n in 0..domain.size {
            for m in 0..domain.size {
                let v = domain.omega.pow((n * m) as u64);
                values.push(v);
            }
        }
        Self(crate::matrix::Matrix::new(values, domain.size, domain.size))
    }
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

struct CooleyTurkey<F: Field> {
    pows: Vec<F>,
}

impl<F: Field> FFT<F> for CooleyTurkey<F> {
    fn new(domain: EvaluationDomainGenerator<F>) -> Self {
        let mut m = F::one();
        let mut pows = Vec::with_capacity(domain.size);
        pows.push(m);
        for _ in 1..domain.size {
            m = m * domain.omega;
            pows.push(m);
        }
        Self { pows }
    }
    fn fft(&self, values: &[F]) -> Vec<F> {
        let values: Vec<&F> = values.iter().collect::<Vec<_>>();
        let domain: Vec<&F> = self.pows.iter().collect::<Vec<_>>();
        cooley_tukey_fft(&values, &domain)
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

fn split<'a, F: Field>(v: &'a [&'a F], even: bool) -> Vec<&'a F> {
    v.iter()
        .enumerate()
        .filter(|(n, _)| (n % 2 == 0) == even)
        .map(|(_, v)| *v)
        .collect()
}

// https://www.algorithm-archive.org/contents/cooley_tukey/cooley_tukey.html
fn cooley_tukey_fft<F: Field>(vals: &[&F], domain: &[&F]) -> Vec<F> {
    if vals.len() == 1 {
        vec![*vals[0]]
    } else {
        let half_domain: Vec<_> = split(domain, true);
        let l = cooley_tukey_fft(&split(vals, true), &half_domain);
        let r = cooley_tukey_fft(&split(vals, false), &half_domain);

        let mut o = vec![F::zero(); vals.len()];
        for (i, (x, y)) in l.into_iter().zip(r.into_iter()).enumerate() {
            let y_times_root = y * domain[i];
            o[i] = x + y_times_root;
            o[i + vals.len() / 2] = x - y_times_root;
        }
        o
    }
}

// according https://stackoverflow.com/questions/52270320/implementing-fft-over-finite-fields
pub fn mul_ntt<F: Field, FFTI: FFT<F>>(
    fft: FFTI,
    mut a_vals: Vec<F>,
    mut b_vals: Vec<F>,
) -> Vec<F> {
    let sum = a_vals.len() + b_vals.len();
    let zero = F::zero();

    a_vals.extend(vec![F::zero(); sum - a_vals.len()]);
    b_vals.extend(vec![F::zero(); sum - b_vals.len()]);

    let a_freq = fft.fft(&a_vals);
    let b_freq = fft.fft(&b_vals);

    let mut c_freq = Vec::new();

    for n in 0..a_freq.len() {
        let l = a_freq.get(n).unwrap_or(&zero);
        let r = b_freq.get(n).unwrap_or(&zero);
        c_freq.push(*l * r);
    }

    fft.fft_inv(&c_freq)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::U64Field;

    #[test]
    fn test_fft_vandermonde() {
        type F = U64Field<337>;
        let ev = EvaluationDomainGenerator::new(F::from(85u64), 8);
        let domain = VandermondeMatrix::new(ev);
        let values = [3, 1, 4, 1, 5, 9, 2, 6].map(|x| F::from(x as u64));
        let current_freq = domain.fft(&values);

        let expected_freq = [31, 70, 109, 74, 334, 181, 232, 4].map(|x| F::from(x as u64));
        assert_eq!(&current_freq, &expected_freq[..]);

        let current_values = domain.fft_inv(&current_freq);
        assert_eq!(&current_values, &values[..]);
    }

    #[test]
    fn test_fft_cooley_turkey() {
        type F = U64Field<337>;
        let ev = EvaluationDomainGenerator::new(F::from(85u64), 8);
        let ct = CooleyTurkey::new(ev);

        let values = [3, 1, 4, 1, 5, 9, 2, 6].map(|x| F::from(x as u64));
        let current_freq = ct.fft(&values);

        let expected_freq = [31, 70, 109, 74, 334, 181, 232, 4].map(|x| F::from(x as u64));
        assert_eq!(&current_freq, &expected_freq[..]);

        let current_values = ct.fft_inv(&current_freq);
        assert_eq!(&current_values, &values[..]);
    }
    
    #[test]
    fn test_ntt_poly_mul() {
        type F = U64Field<337>;

        let a: Vec<F> = [24, 12, 28, 8].iter().map(|x| F::from(*x as u64)).collect();
        let b: Vec<F> = [4, 26, 29, 23].iter().map(|x| F::from(*x as u64)).collect();
        let poly_c = Poly::new(a.clone()) * Poly::new(b.clone());

        let ev = EvaluationDomainGenerator::new(F::from(85u64), 8);
        let domain = CooleyTurkey::new(ev);
        let ntt_c = Poly::new(mul_ntt(domain, a, b));

        assert_eq!(poly_c, ntt_c);
    }

}
