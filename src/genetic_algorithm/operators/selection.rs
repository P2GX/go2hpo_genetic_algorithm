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
    fn select<'a>(&self, population: &'a Vec<Solution<T>>) -> &'a Solution<T>;
}


pub struct TournamentSelection{
    tournament_size: usize,
}

impl TournamentSelection {
    pub fn new(tournament_size: usize) -> Self{
        if tournament_size < 1{
            panic!("tournament_size should be at least 1. Current value: {}", tournament_size);
        }
        Self{tournament_size}
    }
}

impl<T> Selection<T> for TournamentSelection
where 
T: Clone,
{
    fn select<'a>(&self, population: &'a Vec<Solution<T>>) -> &'a Solution<T> {
        let mut rng = rand::rng();
        let rand_index = rng.random_range(0..population.len());
        let mut best: &Solution<T> = population.get(rand_index).expect("The index should not be out of bounds");

        for i in 1..self.tournament_size{
            let rand_index = rng.random_range(0..population.len());
            let opponent: &Solution<T> = population.get(rand_index).expect("The index should not be out of bounds");
            
            if best < opponent{
                best = opponent;
            }
        }

        best

    }
}

pub struct RouletteWheelSelection
{
    transform: Box<dyn Fn(f64) -> f64>,
}

impl<T> Selection<T> for RouletteWheelSelection
where 
T: Clone,
{
    fn select<'a>(&self, population: &'a Vec<Solution<T>>) -> &'a Solution<T> {

        let tot: f64 = population
        .iter()
        .map(|sol| sol.get_score())
        .map(|score| (self.transform)(score))
        .sum();
        
        let mut rng = rand::rng();
        let mut threshold = rng.random_range(0.0..tot);

        for sol in population{
            threshold -= (self.transform)(sol.get_score());
            
            if threshold <= 0.0 {
                return sol;
            }
        }
        &population[population.len() - 1]
    }
}

impl RouletteWheelSelection{
    pub fn new(transform: Box<dyn Fn(f64) -> f64>) -> Self {
        Self { transform }
    }

    pub fn default() -> Self {
        let identity = Box::new(|x: f64| x);
        RouletteWheelSelection::new(identity)
    }
}


pub struct RankSelection;

impl<T> Selection<T> for RankSelection
where 
T: Clone,
{
    fn select<'a>(&self, population: &'a Vec<Solution<T>>) -> &'a Solution<T> {
        let mut ranked_population: Vec<_> = population.iter().collect();
        ranked_population.sort_unstable_by(|a, b| a.get_score().partial_cmp(&b.get_score()).expect("It should be possible to compare the values"));

        // Sum of the first (ranked_population.len() - 1) natural numbers, which is the sum of all the indexes / ranks
        let tot: usize = ((ranked_population.len() - 1) * ranked_population.len()) / 2;

        let mut rng = rand::rng();
        let mut threshold = rng.random_range(0..tot);

        for (rank, &sol) in ranked_population.iter().enumerate(){
            threshold -= rank;
            if threshold <= 0{
                return sol;
            }
        }
        &ranked_population[ranked_population.len() - 1]
    }
}
