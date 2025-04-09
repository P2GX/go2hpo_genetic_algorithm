use crate::genetic_algorithm::Solution;
use crate::logical_formula::ConjunctionGenerator;
use num::Num;
use rand::distr::uniform::SampleUniform;
use rand::{prelude::*, rng};
use std::iter::Sum;
use std::ops::SubAssign;

use ontolius::term::simple::SimpleMinimalTerm;
use ontolius::{
    ontology::{HierarchyWalks, OntologyTerms},
    TermId,
};

use crate::logical_formula::DNFVec;
use crate::{
    logical_formula::{Conjunction, TermObservation},
    logical_formula::{DNFBitmask, DNF},
};

//SELECTION OPERATOR
// TRAIT: Selection
// Implementations:
//      TournamentSelection
//      RouletteWheelSelection
//      RankSelection

pub trait Selection<T> {
    fn select<'a>(&mut self, population: &'a Vec<Solution<T>>) -> &'a Solution<T>;
}

pub struct TournamentSelection<'a, R: Rng> {
    tournament_size: usize,
    rng: &'a mut R,
}

impl<'a, R> TournamentSelection<'a, R>
where
    R: Rng,
{
    pub fn new(tournament_size: usize, rng: &'a mut R) -> Self {
        if tournament_size < 1 {
            panic!(
                "tournament_size should be at least 1. Current value: {}",
                tournament_size
            );
        }
        Self {
            tournament_size,
            rng,
        }
    }
}

impl<'a, T, R> Selection<T> for TournamentSelection<'a, R>
where
    T: Clone,
    R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {
        let rand_index = self.rng.random_range(0..population.len());
        let mut best: &Solution<T> = population
            .get(rand_index)
            .expect("The index should not be out of bounds");

        for i in 1..self.tournament_size {
            let rand_index = self.rng.random_range(0..population.len());
            let opponent: &Solution<T> = population
                .get(rand_index)
                .expect("The index should not be out of bounds");

            if best < opponent {
                best = opponent;
            }
        }

        best
    }
}

pub struct RouletteWheelSelection<'a, R: Rng> {
    rng: &'a mut R,
    transform: Box<dyn Fn(f64) -> f64>,
    roulette_wheel: &'a mut RouletteWheel<f64>,
}

impl<'a, T, R> Selection<T> for RouletteWheelSelection<'a, R>
where
    T: Clone,
    R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {

        // check that the Roulette Wheel is initialized
        match self.roulette_wheel.get_tot() {
            Some(_) => {},
            None => {
                let scores: Vec<f64> = population
                    .iter()
                    .map(|sol| sol.get_score())
                    .map(|score| (self.transform)(score))
                    .collect();
                self.roulette_wheel.initialize(scores);
            }
        };

        let index = self.roulette_wheel.select_random_index(self.rng);
        &population[index]
    }
}

impl<'a, R> RouletteWheelSelection<'a, R>
where
    R: Rng,
{
    pub fn new(rng: &'a mut R, transform: Box<dyn Fn(f64) -> f64>, roulette_wheel: &'a mut RouletteWheel<f64>) -> Self {
        Self { rng, transform, roulette_wheel }
    }

    pub fn default(rng: &'a mut R, roulette_wheel: &'a mut RouletteWheel<f64>) -> Self {
        let identity = Box::new(|x: f64| x);
        RouletteWheelSelection::new(rng, identity, roulette_wheel)
    }
}

pub struct RankSelection<'a, R: Rng> {
    rng: &'a mut R,
}

impl<'a, R, T> Selection<T> for RankSelection<'a, R>
where
    T: Clone,
    R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {
        let mut ranked_population: Vec<_> = population.iter().collect();
        ranked_population.sort_unstable_by(|a, b| {
            a.get_score()
                .partial_cmp(&b.get_score())
                .expect("It should be possible to compare the values")
        });

        // Sum of the first (ranked_population.len() - 1) natural numbers, which is the sum of all the indexes / ranks
        let tot: usize = ((ranked_population.len() - 1) * ranked_population.len()) / 2;

        let mut threshold = self.rng.random_range(0..tot);

        for (rank, &sol) in ranked_population.iter().enumerate() {
            threshold -= rank;
            if threshold <= 0 {
                return sol;
            }
        }
        &ranked_population[ranked_population.len() - 1]
    }
}

pub struct RouletteWheel<T> {
    scores: Vec<T>,
    tot: Option<T>,
}

impl<T> RouletteWheel<T>
where
    T: Num + Copy + Sum<T> + PartialOrd + SampleUniform + SubAssign,
{
    pub fn new() -> Self {
        RouletteWheel {
            scores: Vec::new(),
            tot: None,
        }
    }

    pub fn get_tot(&self) -> Option<T> {
        return self.tot;
    }

    pub fn initialize(&mut self, scores: Vec<T>) {
        self.scores = scores;
        self.tot = Some(self.scores.iter().copied().sum());
    }

    /// get the index of the randomly selected element
    pub fn select_random_index<R: Rng>(&self, rng: &mut R) -> usize {
        let total = self
            .tot
            .expect("Tot should be already initialized before calling select");

        let mut threshold = rng.random_range(T::zero()..total);

        for (index, &score) in self.scores.iter().enumerate() {
            threshold -= score;
            if threshold <= T::zero() {
                return index;
            }
        }

        return (self.scores.len() - 1);
    }

    pub fn get_scores(&self) -> &Vec<T> {
        return &self.scores;
    }
}
