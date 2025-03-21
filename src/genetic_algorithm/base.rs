use crate::{
    logical_formula::{Conjunction, TermObservation},
    logical_formula::{DNFBitmask, DNF},
    logical_formula::SatisfactionChecker,
};

//Solution, Individual, Chromosome, GeneticIndividual, FormulaWrapper
// It's just a wrapper over the formula (DNF or Conjunction) to avoid computing the same fitness function more than once 
// and to keep track of the fitness score withoutthe need of an additional vector
//make it implement clone
#[derive(Debug, Clone)]
pub struct Solution<T: Clone> {
    pub formula: T,
    score: Option<f64>,
}

impl<T: Clone> Solution<T> {
    pub fn get_or_score(&mut self, scorer: &dyn FitnessScorer<T>) -> f64 {
        match self.score {
            Some(val) => return val,
            None => {
                self.score = Some(scorer.fitness(&self.formula));
                return self.score.unwrap();
            }
        }
    }

    pub fn get(&self) -> Option<f64> {
        return self.score.clone();
    }
}

// FITNESS SCORER

pub trait FitnessScorer<T> {
    fn fitness(&self, formula: &T) -> f64;
}

pub struct DNFScorer<C: SatisfactionChecker> {
    checker: C,
}

impl<C: SatisfactionChecker, T: DNF> FitnessScorer<T> for DNFScorer<C> {
    fn fitness(&self, formula: &T) -> f64 {
        1.0
    }
}

pub struct ConjunctionScorer<C: SatisfactionChecker> {
    checker: C,
}

impl<C: SatisfactionChecker> FitnessScorer<Conjunction> for ConjunctionScorer<C> {
    fn fitness(&self, formula: &Conjunction) -> f64 {
        1.0
    }
}