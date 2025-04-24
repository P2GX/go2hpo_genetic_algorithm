use hpo2gene_mapper::GenePhenotypeMapping;
use ontolius::TermId;

use crate::{
    annotations::GeneSetAnnotations, logical_formula::{Conjunction, SatisfactionChecker, DNF}
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
pub struct FormulaEvaluator<T, P>{
    scorer: Box<dyn FitnessScorer<T, P>>,
}

impl<T, P> FormulaEvaluator<T, P> {
    pub fn evaluate(&self, formula: T, phenotype: P) -> Solution<T>{
        let score = self.scorer.fitness(&formula, &phenotype);
        Solution::new(formula, score)
    } 
}


// FITNESS SCORER

pub trait FitnessScorer<T, P> {
    fn fitness(&self, formula: &T, phenotype: &P) -> f64; //I pass the phenotype by argument, that very likely will always be a HPO TermId term. In this way sperimenting with different phenotypes can be easier
}

// pub trait ConjunctionScorer: FitnessScorer<Conjunction, TermId> {}


pub struct DNFScorer<C> {
    checker: C,
    gene_set : &'static GeneSetAnnotations,
}

impl<C: SatisfactionChecker, T: DNF> FitnessScorer<T, TermId> for DNFScorer<C> {
    fn fitness(&self, formula: &T, phenotype: &TermId) -> f64 {
        todo!()
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