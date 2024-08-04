
use blstrs::{pairing, Bls12, G1Affine, G1Projective, G2Affine, G2Prepared, G2Projective, Gt, Scalar};
use group::{ff::Field as _, Curve as _, Group};
use pairing::{MultiMillerLoop, MillerLoopResult};
use rand::{thread_rng, RngCore};
use std::{ops::Mul, time::{Duration, Instant}};
use std::hint::black_box;


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

trait MultiExp : Sized {
    fn multi_exp(bases: &[Self], scalars: &[Scalar]) -> Self;
}

impl MultiExp for G1Projective {
    fn multi_exp(bases: &[Self], scalars: &[Scalar]) -> G1Projective {
        G1Projective::multi_exp(bases, scalars)
    }
}

impl MultiExp for G2Projective {
    fn multi_exp(bases: &[Self], scalars: &[Scalar]) -> G2Projective {
        G2Projective::multi_exp(bases, scalars)
    }
}

trait GroupOpsSimulationItem {
    fn simulate(&self);
}


struct Exps<T>
where T : Mul<Scalar> + Group {
    bases : Vec<T>,
    scalars : Vec<Scalar>
}


impl<T> Exps<T>
where T : Mul<Scalar> + Group {
    fn new(mut rng: &mut impl RngCore, num : usize) -> Self {
        Self {
            bases: (0..num).map(|_| T::random(&mut rng)).collect(),
            scalars: (0..num).map(|_| Scalar::random(&mut rng)).collect(),
        }
    }
}

impl<T> GroupOpsSimulationItem for Exps<T>
where T : Mul<Scalar> + Group {
    fn simulate(&self) {
        // use black_box here so the rust compiler doesn't optimize away the "dead" code
        // https://doc.rust-lang.org/stable/std/hint/fn.black_box.html
        let result : Vec<<T as Mul<Scalar>>::Output> = black_box(self.bases
            .iter()
            .zip(&self.scalars)
            .map(|(base, scalar)| *base * *scalar)
            .collect());
    }
}

type G1Exps = Exps<G1Projective>;
type G2Exps = Exps<G2Projective>;

struct MultiExps<T>
where T : Mul<Scalar> + Group {
    num: usize,
    bases : Vec<T>,
    scalars : Vec<Scalar>
}

impl<T> MultiExps<T>
where T : MultiExp + Group + Mul<Scalar> {
    fn new(mut rng: &mut impl RngCore, num : usize, size: usize) -> Self {
        Self {
            num,
            bases: (0..size).map(|_| T::random(&mut rng)).collect(),
            scalars: (0..size).map(|_| Scalar::random(&mut rng)).collect(),
        }
    }
}

impl<T> GroupOpsSimulationItem for MultiExps<T>
where T : MultiExp + Group + Mul<Scalar> {
    fn simulate(&self) {
        // use black_box here so the rust compiler doesn't optimize away the "dead" code
        // https://doc.rust-lang.org/stable/std/hint/fn.black_box.html
        let result : Vec<T> = black_box(
            (0..self.num).map(|_| T::multi_exp(&self.bases, &self.scalars)).collect()
            );
    }
}

type G1MultiExps = MultiExps<G1Projective>;
type G2MultiExps = MultiExps<G2Projective>;

struct Pairings {
    args_G1: Vec<G1Affine>,
    args_G2: Vec<G2Affine>,
}


impl Pairings {
    fn new(mut rng: &mut impl RngCore, num : usize) -> Self {
        Self {
            args_G1: (0..num).map(|_| G1Projective::random(&mut rng).to_affine()).collect(),
            args_G2: (0..num).map(|_| G2Projective::random(&mut rng).to_affine()).collect(),
        }
    }
}

impl GroupOpsSimulationItem for Pairings {
    fn simulate(&self) {
        // use black_box here so the rust compiler doesn't optimize away the "dead" code
        // https://doc.rust-lang.org/stable/std/hint/fn.black_box.html
        let result : Vec<Gt> = black_box(
            self.args_G1.iter()
            .zip(&self.args_G2)
            .map(|(a, b)| pairing(a,b))
            .collect()
            );
    }
}

struct MultiPairings {
    num: usize,
    args_G1: Vec<G1Projective>,
    args_G2: Vec<G2Projective>,
}

impl MultiPairings {
    fn new(mut rng: &mut impl RngCore, num : usize, size: usize) -> Self {
        Self {
            num,
            args_G1: (0..size).map(|_| G1Projective::random(&mut rng)).collect(),
            args_G2: (0..size).map(|_| G2Projective::random(&mut rng)).collect(),
        }
    }
}

impl GroupOpsSimulationItem for MultiPairings {
    fn simulate(&self) {
        // use black_box here so the rust compiler doesn't optimize away the "dead" code
        // https://doc.rust-lang.org/stable/std/hint/fn.black_box.html
        for i in 0..self.num {
            let result : Gt = black_box(
                multi_pairing(self.args_G1.iter(), self.args_G2.iter())
                );
        }
    }
}


pub struct GroupOpsSimulation<'a, R>
where R : RngCore {
    items: Vec<Box<dyn GroupOpsSimulationItem>>,
    rng: &'a mut R 

}


impl<'a, R> GroupOpsSimulation<'a, R>
where R : RngCore
{
    pub fn new(rng: &'a mut R) -> Self {
        Self {
            items: Vec::new(),
            rng
        }
    }

    pub fn simulate(&self) {
        for item in &self.items {
            item.simulate();
        }
    }

    // convenience methods
    pub fn g1_exps(&mut self, num: usize) -> &mut Self {
        self.items.push(Box::new(G1Exps::new(self.rng, num)));
        self
    }
    pub fn g2_exps(&mut self, num: usize) -> &mut Self {
        self.items.push(Box::new(G2Exps::new(self.rng, num)));
        self
    }
    pub fn g1_multi_exps(&mut self, num: usize, size: usize) -> &mut Self {
        self.items.push(Box::new(G1MultiExps::new(self.rng, num, size)));
        self
    }
    pub fn g2_multi_exps(&mut self, num: usize, size: usize) -> &mut Self {
        self.items.push(Box::new(G2MultiExps::new(self.rng, num, size)));
        self
    }
    pub fn pairings(&mut self, num: usize) -> &mut Self {
        self.items.push(Box::new(Pairings::new(self.rng, num)));
        self
    }
    pub fn multi_pairings(&mut self, num: usize, size: usize) -> &mut Self {
        self.items.push(Box::new(MultiPairings::new(self.rng, num, size)));
        self
    }
}
