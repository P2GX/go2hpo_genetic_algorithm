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
    fn select(&mut self, population: &Vec<Solution<T>>) -> Solution<T>;
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
    fn select(&mut self, population: & Vec<Solution<T>>) -> Solution<T>{
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
        best.clone()
    }
}

pub struct RouletteWheelSelection<'a, T, R: Rng> {
    rng: &'a mut R,
    transform: Box<dyn Fn(f64) -> f64>,
    
    ref_population: Option<*const Vec<Solution<T>>>,
    scores: Vec<f64>,
    tot: Option<f64>,
}

impl<'a, T, R> Selection<T> for RouletteWheelSelection<'a,T,  R>
where
    T: Clone,
    R: Rng,
{
    fn select(&mut self, population: &Vec<Solution<T>>) -> Solution<T> 
    {
        if !self.is_population_same(population) {
            self.initialize(population);
        }

        let index = self.select_random_index();
        population[index].clone()
    }
}

impl<'a, T, R> RouletteWheelSelection<'a, T, R>
where
    R: Rng,
{
    pub fn new(rng: &'a mut R, transform: Box<dyn Fn(f64) -> f64>) -> Self {
        Self { rng, transform, ref_population: None, scores: Vec::new(), tot: None}
    }

    pub fn default(rng: &'a mut R) -> Self {
        let identity = Box::new(|x: f64| x);
        RouletteWheelSelection::new(rng, identity)
    }

    pub fn initialize(&mut self, population: &Vec<Solution<T>>) {
        self.ref_population = Some(population);
        self.scores = population
        .iter()
        .map(|sol| sol.get_score())
        .map(|score| (self.transform)(score))
        .collect();
        self.tot = Some(self.scores.iter().copied().sum());
    }

    // check if ref_population and population point to the same vector
    pub fn is_population_same(&self, population: & Vec<Solution<T>>) -> bool {
        match self.ref_population {
            Some(ptr) => std::ptr::eq(ptr, population),
            None => false,
        }
    }

    fn select_random_index(&mut self) -> usize{
        let mut threshold = self.rng.random_range(0.0..self.tot.expect("It should be Some"));

        for (index, &score) in self.scores.iter().enumerate() {
            threshold -= score;
            if threshold <= 0.0 {
                return index;
            }
        }
        return (self.scores.len() - 1);
    }
    
}

pub struct RankSelection<'a, T, R: Rng> {
    rng: &'a mut R,

    ref_population: Option<*const Vec<Solution<T>>>,
    ranked_population: Option<Vec<Solution<T>>>,
    tot: Option<usize>,
}

impl<'a, R, T> Selection<T> for RankSelection<'a, T, R>
where
    T: Clone,
    R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> Solution<T> {
        if !self.is_population_same(population) {
            self.initialize(population);
        }

        let mut threshold = self.rng.random_range(0..self.tot.expect("Tot should be Some"));
        let ranked = self.ranked_population.as_ref().expect("ranked_population should be Some");
        for (rank, sol) in ranked.iter().enumerate() {
            threshold -= rank;
            if threshold <= 0 {
                return sol.clone();
            }
        }
        ranked[ranked.len() - 1].clone()
    }
}

impl<'a, T, R> RankSelection<'a, T, R>
where 
    R: Rng,
    T: Clone,
{
    pub fn new(rng: &'a mut R) -> Self {
        Self { rng: rng, ref_population: None, ranked_population: None, tot: None }
    }

    pub fn initialize(&mut self, population: &Vec<Solution<T>>) {
        self.ref_population = Some(population);
        
        let mut ranked: Vec<Solution<T>> = population.iter().cloned().collect();
        ranked.sort_unstable_by(|a, b| {
            a.get_score()
                .partial_cmp(&b.get_score())
                .expect("It should be possible to compare the values")
        });
        self.ranked_population = Some(ranked);

        self.tot = Some(((population.len() - 1) * population.len()) / 2);
    }

    // check if ref_population and population point to the same vector
    pub fn is_population_same(&self, population: &Vec<Solution<T>>) -> bool {
        match self.ref_population {
            Some(ptr) => std::ptr::eq(ptr, population),
            None => false,
        }
    }
}