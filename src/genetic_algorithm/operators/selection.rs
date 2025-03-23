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

pub trait Selection<T> 
where 
T: Clone,
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


