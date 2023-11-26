use num_complex::Complex;
use std::f64::consts::PI;
use std::time::Instant;


fn start_points(coefficients: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let degree = coefficients.len() as f64 - 1.0;
    let last_non_zero = coefficients
        .iter()
        .rev()
        .find(|&n| n.norm() != 0.0)
        .unwrap();
    let radius = coefficients[0] / last_non_zero.powf(1.0 / degree);

    (0..degree as i64)
        .map(|n| 2.0 * PI * n as f64 / degree + PI / 2.0 * degree)
        .map(|arg| Complex::new(arg.cos(), arg.sin()) * radius)
        .collect()
}


fn differentiate(coefficients: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let degree = coefficients.len() as f64 - 1.0;

    coefficients
        .iter()
        .enumerate()
        .map(|(i, &coeff)| coeff * (degree - i as f64))
        .collect()
}


fn aberth_ehrlich(coefficients: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let degree = coefficients.len() as f64 - 1.0;
    let derivative = differentiate(&coefficients);
    let mut current_guesses = start_points(&coefficients);
    let mut old_guesses: Vec<Complex<f64>> = vec![Complex::default(); degree as usize - 1];
    let mut iteration_count = 0;

    while old_guesses
        .iter()
        .zip(current_guesses.iter())
        .any(|(&old_guess, &current_guess)| (old_guess - current_guess).norm() > 0.000001)
    {
        iteration_count += 1;
        let mut new_guesses: Vec<Complex<f64>> = Default::default();

        for &k in &current_guesses {
            let dividend: Complex<f64> = coefficients
                .iter()
                .enumerate()
                .fold(Complex::default(), |accumulator, (i, &coeff)| {
                    accumulator + coeff * k.powc(Complex::new(degree - i as f64, 0.0))
                });

            let divisor: Complex<f64> = derivative
                .iter()
                .enumerate()
                .fold(Complex::default(), |accumulator, (i, &coeff)| {
                    accumulator + coeff * k.powc(Complex::new(degree - 1.0 - i as f64, 0.0))
                });

            let quotient = dividend / divisor;

            let summation = current_guesses
                .iter()
                .filter(|&j| *j != k)
                .fold(Complex::default(), |accumulator, &j| accumulator + 1.0 / (k - j));

            let kth_guess = quotient / (1.0 - quotient * summation);
            new_guesses.push(kth_guess);
        }

        old_guesses = current_guesses.clone();
        current_guesses = current_guesses
            .iter()
            .zip(new_guesses.iter())
            .map(|(&i, &new_guess)| i - new_guess)
            .collect();
    }

    println!("DEBUG {} iter", iteration_count);
    current_guesses
}


fn main() {
    let coefficients = vec![1, 2, 3, 4];

    let coefficients_complex =
        coefficients.iter().map(|&a| Complex::new(a as f64, 0.0)).collect();

    let start_time = Instant::now();
    let result = aberth_ehrlich(coefficients_complex);
    let end_time = Instant::now();

    let elapsed_time = end_time - start_time;

    println!("DEBUG {:.2?}\n{:.5?}", elapsed_time, result);
}