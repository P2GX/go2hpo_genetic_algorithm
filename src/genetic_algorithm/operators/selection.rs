use crate::genetic_algorithm::Solution;
use rand::{prelude::*, rng};
use crate::logical_formula::ConjunctionGenerator;

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



pub trait Selection<T>
{
    fn select<'a>(&mut self, population: &'a Vec<Solution<T>>) -> &'a Solution<T>;
}


pub struct TournamentSelection <'a, R: Rng>{
    tournament_size: usize,
    rng: &'a mut R,
}

impl<'a, R> TournamentSelection <'a, R> 
where 
R: Rng,
{
    pub fn new(tournament_size: usize, rng:  &'a mut R) -> Self{
        if tournament_size < 1{
            panic!("tournament_size should be at least 1. Current value: {}", tournament_size);
        }
        Self{tournament_size, rng}
    }
}

impl<'a, T, R> Selection<T> for TournamentSelection<'a, R>
where 
T: Clone,
R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {
        let rand_index = self.rng.random_range(0..population.len());
        let mut best: &Solution<T> = population.get(rand_index).expect("The index should not be out of bounds");

        for i in 1..self.tournament_size{
            let rand_index = self.rng.random_range(0..population.len());
            let opponent: &Solution<T> = population.get(rand_index).expect("The index should not be out of bounds");
            
            if best < opponent{
                best = opponent;
            }
        }

        best

    }
}

pub struct RouletteWheelSelection<'a, R: Rng>
{
    rng: &'a mut R,
    transform: Box<dyn Fn(f64) -> f64>,
}

impl<'a, T, R> Selection<T> for RouletteWheelSelection<'a, R>
where 
T: Clone,
R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {

        let tot: f64 = population
        .iter()
        .map(|sol| sol.get_score())
        .map(|score| (self.transform)(score))
        .sum();
        
        let mut threshold = self.rng.random_range(0.0..tot);

        for sol in population{
            threshold -= (self.transform)(sol.get_score());
            
            if threshold <= 0.0 {
                return sol;
            }
        }
        &population[population.len() - 1]
    }
}

impl<'a, R> RouletteWheelSelection<'a, R>
where 
R: Rng,{
    pub fn new(rng:  &'a mut R, transform: Box<dyn Fn(f64) -> f64>) -> Self {
        Self { rng, transform }
    }

    pub fn default(rng:  &'a mut R) -> Self {
        let identity = Box::new(|x: f64| x);
        RouletteWheelSelection::new(rng, identity)
    }
}


pub struct RankSelection<'a, R: Rng>{
    rng: &'a mut R,
}

impl<'a, R, T> Selection<T> for RankSelection<'a, R>
where 
T: Clone,
R: Rng,
{
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> &'b Solution<T> {
        let mut ranked_population: Vec<_> = population.iter().collect();
        ranked_population.sort_unstable_by(|a, b| a.get_score().partial_cmp(&b.get_score()).expect("It should be possible to compare the values"));

        // Sum of the first (ranked_population.len() - 1) natural numbers, which is the sum of all the indexes / ranks
        let tot: usize = ((ranked_population.len() - 1) * ranked_population.len()) / 2;

        let mut threshold = self.rng.random_range(0..tot);

        for (rank, &sol) in ranked_population.iter().enumerate(){
            threshold -= rank;
            if threshold <= 0{
                return sol;
            }
        }
        &ranked_population[ranked_population.len() - 1]
    }
}
