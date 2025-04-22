use crate::{
    logical_formula::Conjunction,
    logical_formula::DNF,
    logical_formula::SatisfactionChecker,
};


//Solution, Individual, Chromosome, GeneticIndividual, FormulaWrapper
// It's just a wrapper over the formula (DNF or Conjunction) to avoid computing the same fitness function more than once 
// and to keep track of the fitness score without the need of an additional vector
//make it implement clone
#[derive(Debug, Clone)]
pub struct Solution<T> {
    pub formula: T,
    score: f64,
}


impl<T> Solution<T>{
    pub fn new(formula: T, score: f64) -> Solution<T>{
        Solution {formula, score}
    }
}

impl<T> Solution<T> {

    pub fn get_score(&self) -> f64 {
        return self.score;
    }

    pub fn get_formula(&self) -> &T{
        return &self.formula;
    }
}

impl<T> PartialEq for Solution<T> {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}

impl<T> PartialOrd for Solution<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.score.partial_cmp(&other.score)
    }
}

// Possible names:
// FormulaEvaluator and the method is "evaluate" which returns a Solution
// or 
// SolutionGenerator and the method is "generate" which returns a Solution
pub struct FormulaEvaluator<T>{
    scorer: Box<dyn FitnessScorer<T>>,
}

impl<T> FormulaEvaluator<T> {
    pub fn evaluate(&self, formula: T) -> Solution<T>{
        let score = self.scorer.fitness(&formula);
        Solution::new(formula, score)
    } 
}


// FITNESS SCORER

pub trait FitnessScorer<T> {
    fn fitness(&self, formula: &T) -> f64;
}

pub struct DNFScorer<C> {
    checker: C,
}

impl<C: SatisfactionChecker, T: DNF> FitnessScorer<T> for DNFScorer<C> {
    fn fitness(&self, formula: &T) -> f64 {
        1.0
    }
}

pub struct ConjunctionScorer<C> {
    checker: C,
}

impl<C: SatisfactionChecker> FitnessScorer<Conjunction> for ConjunctionScorer<C> {
    fn fitness(&self, formula: &Conjunction) -> f64 {
        1.0
    }
}

pub struct AccuracyConjunctionScorer<C> { //<C, GO, EXPR, HPO>
    checker: C,
    

}