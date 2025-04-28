use std::collections::HashMap;

use hpo2gene_mapper::GenePhenotypeMapping;
use ontolius::TermId;

use crate::{
    annotations::{GeneId, GeneSetAnnotations}, logical_formula::{Conjunction, SatisfactionChecker, DNF}
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

impl<T, P> FormulaEvaluator<T, P> 
where 
T: Clone{
    pub fn evaluate(&self, formula: &T, phenotype: &P) -> Solution<T>{
        let score = self.scorer.fitness(formula, phenotype);
        Solution::new(formula.clone(), score)
    } 
}


// FITNESS SCORER

pub trait FitnessScorer<T, P> {
    fn fitness(&self, formula: &T, phenotype: &P) -> f64; //I pass the phenotype by argument, that very likely will always be a HPO TermId term. In this way sperimenting with different phenotypes can be easier
}

pub enum ScoreMetric {
    Accuracy,
    Precision,
    Recall,
    FScore(f64),
}


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
    gene_set : &'static GeneSetAnnotations,
    score_metric: ScoreMetric,
}

impl<C: SatisfactionChecker> FitnessScorer<Conjunction, TermId> for ConjunctionScorer<C> {
    fn fitness(&self, formula: &Conjunction, phenotype: &TermId) -> f64 {
        let genes_satisfaction: HashMap<GeneId, bool> = self.checker.all_satisfactions(formula);
        match self.score_metric {
            ScoreMetric::Accuracy => self.accuracy(&genes_satisfaction, phenotype),
            ScoreMetric::Precision => self.precision(&genes_satisfaction, phenotype),
            ScoreMetric::Recall => self.recall(&genes_satisfaction, phenotype),
            ScoreMetric::FScore(beta) => self.FScore(&genes_satisfaction, phenotype),
        }
    }
}

impl<C>  ConjunctionScorer<C>{
    pub fn new(checker: C, gene_set : &'static GeneSetAnnotations, score_metric: ScoreMetric) -> Self{
        Self { checker, gene_set, score_metric }
    }

    pub fn accuracy(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        let n_tot = self.gene_set.len();
        if n_tot == 0 { panic!("Gene Set is empty") }

        let correct_predictions = self.gene_set.get_gene_annotations_map().iter()
            .filter(|(gene, gene_annotation)| gene_annotation.contains_phenotype(phenotype) == *genes_satisfaction.get(*gene).expect("There should be an entry for the gene"))
            .count();

        return correct_predictions as f64 / n_tot as f64;
    }

    pub fn precision(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        let mut tp = 0; //true positives
        let mut fp = 0; //false positives

        for (gene, gene_annotation) in self.gene_set.get_gene_annotations_map() {
            let predicted = *genes_satisfaction.get(gene).expect("Missing prediction for gene");
            let actual = gene_annotation.contains_phenotype(phenotype);

            if predicted {
                if actual {
                    tp += 1;
                } else {
                    fp += 1;
                }
            }
        }

        if tp + fp == 0 { // denominator is 0
            return 0.0;  //No positive predictions, return 0.0
        }

        return tp as f64 / (tp + fp) as f64;
    }

    pub fn recall(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        let mut tp = 0;  //true positives
        let mut fn_ = 0; //false negatives
    
        for (gene, gene_annotation) in self.gene_set.get_gene_annotations_map() {
            let predicted = *genes_satisfaction.get(gene).expect("Missing prediction for gene");
            let actual = gene_annotation.contains_phenotype(phenotype);
    
            if actual {
                if predicted {
                    tp += 1;
                } else {
                    fn_ += 1;
                }
            }
        }
    
        if tp + fn_ == 0 {
            return 0.0; 
        }
    
        tp as f64 / (tp + fn_) as f64
    }
    
    pub fn FScore(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        if let ScoreMetric::FScore(beta) = self.score_metric {
            let precision = self.precision(genes_satisfaction, phenotype);
            let recall = self.recall(genes_satisfaction, phenotype);
    
            if precision + recall == 0.0 {
                return 0.0;
            }
    
            let beta_sq = beta * beta;
            (1.0 + beta_sq) * (precision * recall) / (beta_sq * precision + recall)
        } else {
            panic!("FScore was called with wrong score metric");
        }
    }

}



#[cfg(test)]
mod tests {
    use crate::logical_formula::{FormulaGenerator, RedundantRandomConjunctionGenerator};

    use super::*;
    use lazy_static::lazy_static;
    use rand::{rng, rngs::SmallRng, SeedableRng};

    lazy_static! {
        static ref small_test_terms: Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"]
            .into_iter()
            .map(|s| s.parse().unwrap())
            .collect();

        static ref small_test_tissues: Vec<String> = vec!["Colon_Transverse_Muscularis".to_string(),
                                                         "Colon_Transverse_Mixed_Cell".to_string(),
                                                          "Colon_Transverse_Muscularis".to_string(),
                                                          "Testis".to_string(),
                                                          "Small_Intestine_Terminal_Ileum_Mixed_Cell".to_string()];
    }

    #[test]
    fn test_fitness_score() {
        let seed = 42;
        // let mut rng = SmallRng::seed_from_u64(seed);
        // let mut generator = RedundantRandomConjunctionGenerator::new( 
        //     2,
        //     &small_test_terms,
        //     2,
        //     &small_test_tissues,
        //     rng);
        //     let conjunction: Conjunction = generator.generate();

    }



}