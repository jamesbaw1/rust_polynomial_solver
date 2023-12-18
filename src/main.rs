use num_complex::Complex;

#[derive(Debug)]
struct Polynomial(Vec<Complex<f64>>);

impl std::ops::Deref for Polynomial {
    type Target = Vec<Complex<f64>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::SubAssign for Polynomial {
    fn sub_assign(&mut self, rhs: Self) {
        assert_eq!(self.len(), rhs.len());
        *self = Self(
            self
            .iter()
            .zip(rhs.iter())
            .map(|(a, b)| a - b)
            .collect()
        )
    }
}

impl Polynomial {
    fn from<T>(coefficients: Vec<T>) -> Self
    where T: Into<f64> {
        Polynomial(
            coefficients
            .into_iter()
            .map(|a| Complex::from(a.into()))
            .collect()
        )
    }

    fn derivative(&self) -> Self {
        Polynomial(
            self
            .iter()
            .enumerate()
            .map(|(n, a)| (n as f64) * a)
            .skip(1)
            .collect()
        )
    }

    fn start_points(&self) -> Self {
        let pi = std::f64::consts::PI;
        let degree = (self.len() - 1) as f64;

        let radius =
        (self.last().unwrap() /
            self
            .iter()
            .find(|&a| a.norm() != 0.0)
            .unwrap()
        )
        .norm()
        .powf(1.0 / degree);

        Polynomial(
            (0..degree as usize)
            .map(|n| 2.0 * pi * n as f64 / degree + pi / 2.0 * degree)
            .map(|arg| Complex::new(arg.cos(), arg.sin()) * radius)
            .collect()
        )
    }

    fn eval(&self, x: Complex<f64>) -> Complex<f64> {
        self
        .iter()
        .rev()
        .fold(Complex::default(), |acc, a| a + x * acc)
    }

    fn solve(&self) -> Self {
        let mut current_guesses = self.start_points();
        let mut ratio = Complex::default();
        let derivative = self.derivative();

        loop {
            current_guesses -=
            Polynomial(
                current_guesses
                .iter()
                .map(|&k| {
                    ratio = self.eval(k) / derivative.eval(k);
                    ratio / (1.0 - ratio *
                        current_guesses
                        .iter()
                        .filter(|&j| *j != k)
                        .fold(Complex::default(), |acc, &j| acc + 1.0 / (k - j))
                    )
                }).collect()
            );

            if ratio.norm() <= 0.0001 { break; }
        }

        current_guesses
    }
}

fn main() {
    let polynomial = Polynomial::from(vec![/*coefficients go here*/]);
    println!("{:#.5?}", polynomial.solve());
}