mod base;
mod operators;
mod algorithm;

pub use base::{Solution, FitnessScorer};

pub use operators::{Crossover, Mutation, ElitesSelector};
// pub use operators::{ConjunctionCrossover};

pub use algorithm::GeneticAlgorithm;