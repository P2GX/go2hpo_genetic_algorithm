use crate::{
    logical_formula::{Conjunction, TermObservation},
    logical_formula::{DNFBitmask, DNF},
    logical_formula::SatisfactionChecker,
};

use std::sync::atomic::Ordering;

//Solution, Individual, Chromosome, GeneticIndividual, FormulaWrapper
// It's just a wrapper over the formula (DNF or Conjunction) to avoid computing the same fitness function more than once 
// and to keep track of the fitness score without the need of an additional vector
//make it implement clone
#[derive(Debug, Clone)]
pub struct Solution<T: Clone> {
    pub formula: T,
    score: f64,
}


impl<T> Solution<T>
where
T: Clone{
    pub fn new(formula: T, score: f64) -> Solution<T>{
        Solution {formula, score}
    }
}

impl<T: Clone> Solution<T> {
    // Deprecated, from now on the Solution will always contain a computed score, so the promise change
    // Now I'm promising always a solution and I don't have to carry the FitnessScorer around to compute
    // the fitness in case in which it wasn't computed
    // pub fn get_or_score(&mut self, scorer: &dyn FitnessScorer<T>) -> f64 {
    //     match self.score {
    //         Some(val) => return val,
    //         None => {
    //             self.score = Some(scorer.fitness(&self.formula));
    //             return self.score.unwrap();
    //         }
    //     }
    // }

    pub fn get_score(&self) -> f64 {
        return self.score;
    }

    pub fn get_formula(&self) -> &T{
        return &self.formula;
    }
}

impl<T: Clone> PartialEq for Solution<T> {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}

impl<T: Clone> PartialOrd for Solution<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.score.partial_cmp(&other.score)
    }
}


// FormulaEvaluator and the method is "evaluate" which returns a Solution
// or 
// SolutionGenerator and the method is "generate" which returns a Solution
pub struct FormulaEvaluator<T>{
    scorer: Box<dyn FitnessScorer<T>>,
}

impl<T> FormulaEvaluator<T>
where
T: Clone{
    pub fn evaluate(&self, formula: T) -> Solution<T>{
        let score = self.scorer.fitness(&formula);
        Solution::new(formula, score)
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