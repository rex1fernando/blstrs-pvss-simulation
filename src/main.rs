
use blstrs::{pairing, Bls12, G1Affine, G1Projective, G2Affine, G2Prepared, G2Projective, Gt, Scalar};
use group::{ff::Field as _, Group as _, Curve as _};
use pairing::{MultiMillerLoop, MillerLoopResult};
use rand::thread_rng;
use std::time::{Instant, Duration};

fn multi_pairing<'a>(lhs: impl Iterator<Item = &'a G1Projective>, rhs: impl Iterator<Item = &'a G2Projective>) -> Gt {
    <Bls12 as MultiMillerLoop>::multi_miller_loop(
        &lhs.zip(rhs)
           .map(|(a,b)| (a.to_affine(), G2Prepared::from(b.to_affine())))
           .collect::<Vec<(G1Affine, G2Prepared)>>()
           .iter()
           .map(|(g1, g2)| (g1, g2))
           .collect::<Vec<(&G1Affine, &G2Prepared)>>()
           ).final_exponentiation()
}

fn simulate_group_ops(
    num_exps_in_G1: usize,
    mexps_in_G1: &[(usize, usize)], 
    num_exps_in_G2: usize,
    mexps_in_G2: &[(usize, usize)],
    num_pairings: usize,
    multi_pairings: &[(usize, usize)]) 
{
    let mut rng = thread_rng();


    let bases_G1 : Vec<G1Projective> = (0..num_exps_in_G1).map(|_| G1Projective::random(&mut rng)).collect();
    let scalars_G1 : Vec<Scalar> = (0..num_exps_in_G1).map(|_| Scalar::random(&mut rng)).collect();
    let exp_G1_args = bases_G1.iter().zip(scalars_G1.iter());

    let mut mexp_G1_args = Vec::new();
    for (num_mexps, mexp_size) in mexps_in_G1 {
        let bases : Vec<G1Projective> = (0..*mexp_size).map(|_| G1Projective::random(&mut rng)).collect();
        let scalars : Vec<Scalar> = (0..*mexp_size).map(|_| Scalar::random(&mut rng)).collect();
        mexp_G1_args.push((*num_mexps, bases, scalars));
    }

    let bases_G2 : Vec<G2Projective> = (0..num_exps_in_G2).map(|_| G2Projective::random(&mut rng)).collect();
    let scalars_G2 : Vec<Scalar> = (0..num_exps_in_G2).map(|_| Scalar::random(&mut rng)).collect();
    let exp_G2_args = bases_G2.iter().zip(scalars_G2.iter());

    let mut mexp_G2_args = Vec::new();
    for (num_mexps, mexp_size) in mexps_in_G2 {
        let bases : Vec<G2Projective> = (0..*mexp_size).map(|_| G2Projective::random(&mut rng)).collect();
        let scalars : Vec<Scalar> = (0..*mexp_size).map(|_| Scalar::random(&mut rng)).collect();
        mexp_G2_args.push((*num_mexps, bases, scalars));
    }

    let pairing_args_G1 : Vec<G1Affine> = (0..num_pairings).map(|_| G1Projective::random(&mut rng).to_affine()).collect();
    let pairing_args_G2 : Vec<G2Affine> = (0..num_pairings).map(|_| G2Projective::random(&mut rng).to_affine()).collect();
    let pairing_args = pairing_args_G1.iter().zip(pairing_args_G2);

    let mut multi_pairing_args = Vec::new();
    for (num_multi_pairings, multi_pairing_size) in multi_pairings {
        let multi_pairing_args_G1 : Vec<G1Projective> = (0..*multi_pairing_size).map(|_| G1Projective::random(&mut rng)).collect();
        let multi_pairing_args_G2 : Vec<G2Projective> = (0..*multi_pairing_size).map(|_| G2Projective::random(&mut rng)).collect();
        multi_pairing_args.push((*num_multi_pairings, multi_pairing_args_G1, multi_pairing_args_G2));
    }

    let start_time = Instant::now();
    let exp_G1_result : Vec<G1Projective> = exp_G1_args.map(|(base, scalar)| base * scalar).collect();
    mexp_G1_args.iter().for_each( |(num_mexps, bases, scalars)| for i in 0..*num_mexps { G1Projective::multi_exp(bases, scalars); });
    let exp_G2_result : Vec<G2Projective> = exp_G2_args.map(|(base, scalar)| base * scalar).collect();
    mexp_G2_args.iter().for_each( |(num_mexps, bases, scalars)| for i in 0..*num_mexps { G2Projective::multi_exp(bases, scalars); });
    let pairings_result : Vec<Gt> = pairing_args.map(|(a,b)| pairing(&a,&b)).collect();
    multi_pairing_args.iter().for_each( |(num_multi_pairings, args_G1, args_G2)| for i in 0..*num_multi_pairings { multi_pairing(args_G1.iter(), args_G2.iter()); });
    let duration = start_time.elapsed();
    println!("{:?}", duration);
}

fn simulate_groth(n: usize, k: usize, t: usize, l: usize) {
    println!("Groth16, n={}, k={}, t={}, l={}", n, k, t, l);
    println!("Prover:");
    simulate_group_ops(
        n + 2*k + l + 2, // G1 exps
        &[(2, n),
          (n*k + n + k + l + 1, 2)], // G1 mexps
        t + 1, // G2 exps
        &[], // no G2 mexps
        0,     // no pairings,
        &[] // no multi-pairings
        );
    println!("Verifier:");
    simulate_group_ops(
        n + 2, // G1 exps
        &[
            (1, 2),
            (2, n+1),
            (n, k+1),
            (1, l+1),
            (1, k*n*l + l + 1),
            (1, n+2),
            (2, k)
        ], // G1 mexps
        1, // G2 exps
        &[
            (1, t+1),
            (1, k)
        ], // G2 mexps
        0,     // no pairings,
        &[(1,3)] // multi-pairings
        );
    // we need to do one inversion to convert the three pairings into multi-pairings, which I'm not
    // simulating.
}

fn main() {
    simulate_groth(1000, 16, 666, 50);
}
