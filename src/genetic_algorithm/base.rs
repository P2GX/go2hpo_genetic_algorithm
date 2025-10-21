use std::collections::HashMap;
use std::fmt;
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

impl<T: fmt::Display> fmt::Display for Solution<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<Formula: {}; Score: {}>", self.formula, self.score)
    }
}

// Possible names:
// FormulaEvaluator and the method is "evaluate" which returns a Solution
// or 
// SolutionGenerator and the method is "generate" which returns a Solution
pub struct FormulaEvaluator<'a, T, P> {
    scorer: Box<dyn FitnessScorer<T, P> + 'a>,
}

impl<'a, T, P> FormulaEvaluator<'a, T, P>
where 
T: Clone{
    /// Constructor for FormulaEvaluator
    pub fn new(scorer: Box<dyn FitnessScorer<T, P> + 'a>) -> Self {
        Self { scorer }
    }

    /// Evaluate a formula with the given phenotype
    pub fn evaluate(&self, formula: &T, phenotype: &P) -> Solution<T>{
        let score = self.scorer.fitness(formula, phenotype);
        Solution::new(formula.clone(), score)
    }

    /// Getter for the scorer
    pub fn get_scorer(&self) -> &dyn FitnessScorer<T, P> {
        &*self.scorer
    }
   
}


// FITNESS SCORER

pub trait FitnessScorer<T, P> {
    fn fitness(&self, formula: &T, phenotype: &P) -> f64; //I pass the phenotype by argument, that very likely will always be a HPO TermId term. In this way sperimenting with different phenotypes can be easier
    fn custom_score_fitness(&self, formula: &T, phenotype: &TermId, custom_score_metric: &ScoreMetric) -> f64;
}

pub enum ScoreMetric {
    Accuracy,
    Precision,
    Recall,
    FScore(f64),
}

pub struct DNFScorer<C> {
    conjunction_scorer: ConjunctionScorer<C>,
    penalty_lambda: f64,
}

impl<C> DNFScorer<C> {
    /// Creates a new [`DNFScorer`].
    ///
    /// # Arguments
    /// * `conjunction_scorer` - Scorer that evaluates how well a conjunction
    ///   explains the gene-to-phenotype associations.
    /// * `penalty_lambda` - Weight of the complexity penalty.  
    ///   Each additional conjunction reduces the final score by
    ///   `penalty_lambda`.
    pub fn new(conjunction_scorer: ConjunctionScorer<C>, penalty_lambda: f64) -> Self {
        Self { conjunction_scorer, penalty_lambda }
    }
}


impl<C: SatisfactionChecker, T: DNF> FitnessScorer<T, TermId> for DNFScorer<C> {
    fn fitness(&self, formula: &T, phenotype: &TermId) -> f64 {
            //self.conjunction_scorer.score_metric
            //self.penalty_lambda 
            self._fitness(formula, phenotype, self.conjunction_scorer.get_metric(), self.penalty_lambda)
    }

    fn custom_score_fitness(&self, formula: &T, phenotype: &TermId, custom_score_metric: &ScoreMetric) -> f64 {
        self._fitness(formula, phenotype, custom_score_metric, self.penalty_lambda)
    }

}


impl<'a, C: SatisfactionChecker> DNFScorer<C> {

    fn _fitness<T: DNF>(&self, formula: &T, phenotype: &TermId, custom_score_metric: &ScoreMetric, penalty_lambda: f64) -> f64{
        let genes_satisfaction: HashMap<GeneId, bool> =
            match self.get_genes_satisfaction_for_dnf(formula) {
                Ok(satisfactions) => satisfactions,
                _ => return 0.0,
            };

        // Clauses = number of conjunctions
        let clauses = formula.get_active_conjunctions().len();
        if clauses == 0 {
            return 0.0; // invalid formula gets zero fitness
        }

        let base_score = match custom_score_metric {
            ScoreMetric::Accuracy => self.conjunction_scorer.accuracy(&genes_satisfaction, phenotype),
            ScoreMetric::Precision => self.conjunction_scorer.precision(&genes_satisfaction, phenotype),
            ScoreMetric::Recall => self.conjunction_scorer.recall(&genes_satisfaction, phenotype),
            ScoreMetric::FScore(beta) => self.conjunction_scorer.fs_score(&genes_satisfaction, phenotype),
        };

        // Simple penalty subtraction
        // base_score - penalty_lambda * (clauses as f64)

        let penalized_score = base_score / (1.0 + penalty_lambda * (clauses as f64 - 1.0));
        penalized_score
    }
}

impl<'a, C: SatisfactionChecker> DNFScorer<C> {
    pub fn get_genes_satisfaction_for_dnf<T: DNF>(&self, formula: &T) -> anyhow::Result<HashMap<GeneId, bool>>{
        let conjunctions = formula.get_active_conjunctions();
        if conjunctions.len() == 0 {anyhow::bail!("The conjunction length shouldn't be 0")}

        let mut conjunction_iter = conjunctions.iter();
        let first_conjunction = conjunction_iter.next().expect("It should have at least the first element");
        let mut satisfaction_mask = self.conjunction_scorer.checker.all_satisfactions(first_conjunction);

        for conjunction in conjunction_iter{
            // Check: IF a Conjunction is made of only NOT terms we get a score of 0
            if !conjunction.tissue_expressions.is_empty() || 
                conjunction.term_observations.iter()
                .map(|obs| obs.is_excluded)
                .any(|excluded| excluded == false){

                for (gene, satisfied) in satisfaction_mask.iter_mut(){
                    
                    if *satisfied {continue;}
                    
                    if self.conjunction_scorer.checker.is_satisfied(gene, conjunction){
                        *satisfied = true;
                    }   
                }
            }
        }

        Ok(satisfaction_mask)
    }

    pub fn change_metric(&mut self, new_score_metric: ScoreMetric) {
        self.conjunction_scorer.score_metric = new_score_metric;
    }

    pub fn get_metric(&self) -> &ScoreMetric {
        &self.conjunction_scorer.score_metric
    }

}


pub struct ConjunctionScorer<C> {
    checker: C,
    score_metric: ScoreMetric,
}


impl<C: SatisfactionChecker> FitnessScorer<Conjunction, TermId> for ConjunctionScorer< C>  {
    //Creare un metodo in privato comune _fitness in comune che viene chiamata dai due metodi
    fn fitness(&self, formula: &Conjunction, phenotype: &TermId) -> f64 {
        let genes_satisfaction: HashMap<GeneId, bool> = self.checker.all_satisfactions(formula);
        match self.score_metric {
            ScoreMetric::Accuracy => self.accuracy(&genes_satisfaction, phenotype),
            ScoreMetric::Precision => self.precision(&genes_satisfaction, phenotype),
            ScoreMetric::Recall => self.recall(&genes_satisfaction, phenotype),
            ScoreMetric::FScore(beta) => self.fs_score(&genes_satisfaction, phenotype),
        }
    }

    fn custom_score_fitness(&self, formula: &Conjunction, phenotype: &TermId, custom_score_metric: &ScoreMetric) -> f64{
        let genes_satisfaction: HashMap<GeneId, bool> = self.checker.all_satisfactions(formula);
        match custom_score_metric {
            ScoreMetric::Accuracy => self.accuracy(&genes_satisfaction, phenotype),
            ScoreMetric::Precision => self.precision(&genes_satisfaction, phenotype),
            ScoreMetric::Recall => self.recall(&genes_satisfaction, phenotype),
            ScoreMetric::FScore(beta) => self.fs_score(&genes_satisfaction, phenotype),
        }
    }
}

impl<C: SatisfactionChecker> ConjunctionScorer<C> {
    pub fn new(checker: C, score_metric: ScoreMetric) -> Self {
        Self { checker, score_metric }
    }

    pub fn change_metric(&mut self, new_score_metric: ScoreMetric) {
        self.score_metric = new_score_metric;
    }

    pub fn get_metric(&self) -> &ScoreMetric {
        &self.score_metric
    }

    pub fn accuracy(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        let n_tot = self.checker.get_gene_set().len();
        if n_tot == 0 { panic!("Gene Set is empty") }

        let correct_predictions = self.checker.get_gene_set().get_gene_annotations_map().iter()
            .filter(|(gene, gene_annotation)| gene_annotation.contains_phenotype(phenotype) == *genes_satisfaction.get(*gene).expect("There should be an entry for the gene"))
            .count();

        return correct_predictions as f64 / n_tot as f64;
    }

    pub fn precision(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
        let mut tp = 0; //true positives
        let mut fp = 0; //false positives

        for (gene, gene_annotation) in self.checker.get_gene_set().get_gene_annotations_map() {
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
    
        for (gene, gene_annotation) in self.checker.get_gene_set().get_gene_annotations_map() {
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
    
    pub fn fs_score(&self, genes_satisfaction: &HashMap<GeneId, bool>, phenotype: &TermId) -> f64{
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
    use super::*;
    use std::collections::{HashMap, HashSet};
    use crate::{annotations::GeneAnnotations, logical_formula::{DNFVec, DgeState, TermObservation, TissueExpression}};
    use ontolius::ontology::csr::MinimalCsrOntology;


        fn term(term_str: &str) -> TermId {
        term_str.parse().unwrap()
    }

    fn make_gene_annotations(
        id: GeneId,
        phenotype_terms: &[&str],
        term_annotations: &[&str],
    ) -> GeneAnnotations {
        GeneAnnotations::new(
            id,
            term_annotations.iter().map(|t| term(t)).collect(),
            HashSet::new(),
            phenotype_terms.iter().map(|p| term(p)).collect(),
        )
    }

    struct DummyChecker {
        predictions: HashMap<String, bool>,
        gene_set_option: usize
    }

    impl SatisfactionChecker for DummyChecker {
        
        fn is_satisfied(&self, symbol: &GeneId, _: &Conjunction) -> bool {
            *self.predictions.get(symbol).unwrap()
        }

        fn all_satisfactions(&self, _: &Conjunction) -> HashMap<String, bool> {
            self.predictions.clone()
        }

        fn get_gene_set(&self) -> &GeneSetAnnotations { 
            if self.gene_set_option == 1{
                return Box::leak(Box::new(GeneSetAnnotations::new({
                        let mut map: HashMap<GeneId, GeneAnnotations> = HashMap::new();

                        // gene1: TP
                        map.insert("gene1".into(), make_gene_annotations("gene1".into(), &["HPO:0000001"], &[]));

                        // gene2: FP
                        map.insert("gene2".into(), make_gene_annotations("gene2".into(), &[], &[]));

                        // gene3: FN
                        map.insert("gene3".into(), make_gene_annotations("gene3".into(), &["HPO:0000001"], &[]));

                        // gene4: TN
                        map.insert("gene4".into(), make_gene_annotations("gene4".into(), &[], &[]));

                        map
                    })));
            }else{
                return Box::leak(Box::new(GeneSetAnnotations::new({
                    let mut map = HashMap::new();
                    map.insert("gene1".to_string(), make_gene_annotations("gene1".into(), &["HPO:0000001", "HPO:0000002", "HPO:0000003"], &[]));
                    map.insert("gene2".to_string(), make_gene_annotations("gene2".into(), &["HPO:0000002", "HPO:0000003"], &[]));
                    map.insert("gene3".to_string(), make_gene_annotations("gene3".into(), &["HPO:0000001", "HPO:0000003"], &[]));
                    map.insert("gene4".to_string(), make_gene_annotations("gene4".into(), &[],&[] ));
                    map
                },
                )));

            }
        }
    }





    #[test]
    fn test_precision_recall_fscore_accuracy() {
        let phenotype = term("HPO:0000001");

        let gene_set = Box::leak(Box::new(GeneSetAnnotations::new({
                let mut map: HashMap<GeneId, GeneAnnotations> = HashMap::new();

                // gene1: TP
                map.insert("gene1".into(), make_gene_annotations("gene1".into(), &["HPO:0000001"], &[]));

                // gene2: FP
                map.insert("gene2".into(), make_gene_annotations("gene2".into(), &[], &[]));

                // gene3: FN
                map.insert("gene3".into(), make_gene_annotations("gene3".into(), &["HPO:0000001"], &[]));

                // gene4: TN
                map.insert("gene4".into(), make_gene_annotations("gene4".into(), &[], &[]));

                map
            })));

        let checker = DummyChecker {
            predictions: HashMap::from([
                ("gene1".into(), true),
                ("gene2".into(), true),
                ("gene3".into(), false),
                ("gene4".into(), false),
            ]),
            gene_set_option: 1
        };

        
        let scorer = ConjunctionScorer::new(checker, ScoreMetric::FScore(1.0));
        let conj = Conjunction::new();

        let all_satisfactions = &scorer.checker.all_satisfactions(&conj);
        
        
        let fscore = scorer.fs_score(all_satisfactions, &phenotype);
        assert!((fscore - 0.5).abs() < 1e-6);

        let precision = scorer.precision(all_satisfactions, &phenotype);
        assert!((precision - 0.5).abs() < 1e-6);

        let recall = scorer.recall(all_satisfactions, &phenotype);
        assert!((recall - 0.5).abs() < 1e-6);

        let accuracy = scorer.accuracy(all_satisfactions, &phenotype);
        assert!((accuracy - 0.5).abs() < 1e-6);

    }

    


    // TESTS PER DNF

    #[test]
    fn test_dnfscorer_with_fixed_predictions() {
        let phenotype1 = term("HPO:0000001");
        let phenotype2 = term("HPO:0000002");
        let phenotype3 = term("HPO:0000003");

        let gene_set = Box::leak(Box::new(GeneSetAnnotations::new({
                let mut map = HashMap::new();
                map.insert("gene1".to_string(), make_gene_annotations("gene1".into(), &["HPO:0000001", "HPO:0000002", "HPO:0000003"], &[]));
                map.insert("gene2".to_string(), make_gene_annotations("gene2".into(), &["HPO:0000002", "HPO:0000003"], &[]));
                map.insert("gene3".to_string(), make_gene_annotations("gene3".into(), &["HPO:0000001", "HPO:0000003"], &[]));
                map.insert("gene4".to_string(), make_gene_annotations("gene4".into(), &[],&[] ));
                map
            },
        )));

        let checker = DummyChecker {
            predictions: HashMap::from([
                ("gene1".into(), true),
                ("gene2".into(), true),
                ("gene3".into(), false),
                ("gene4".into(), false),
            ]),
            gene_set_option: 2
        };

        let conjunction = Conjunction::new(); // dummy, content doesn't matter
        let dnf = DNFVec::from_conjunctions(vec![conjunction]);

        let conjunction_scorer = ConjunctionScorer::new(checker, ScoreMetric::Precision);
        let mut scorer = DNFScorer { conjunction_scorer, penalty_lambda: 0.0 };
        
        // PRECISION
        let precision_1 = scorer.fitness(&dnf, &phenotype1);
        assert!((precision_1 - 0.5).abs() < 1e-6);
        let precision_2 = scorer.fitness(&dnf, &phenotype2);
        assert!((precision_2 - 1.0).abs() < 1e-6);
        let precision_3 = scorer.fitness(&dnf, &phenotype3);
        assert!((precision_3 - 1.0).abs() < 1e-6);

        // RECALL
        scorer.change_metric(ScoreMetric::Recall);
        let recall_1 = scorer.fitness(&dnf, &phenotype1);
        assert!((recall_1 - 0.5).abs() < 1e-6);
        let recall_2 = scorer.fitness(&dnf, &phenotype2);
        assert!((recall_2 - 1.0).abs() < 1e-6);
        let recall_3 = scorer.fitness(&dnf, &phenotype3);
        assert!((recall_3 - (2.0 / 3.0)).abs() < 1e-6);

        // F1
        scorer.change_metric(ScoreMetric::FScore(1.0));
        let f1_1 = scorer.fitness(&dnf, &phenotype1);
        assert!((f1_1 - 0.5).abs() < 1e-6);
        let f1_2 = scorer.fitness(&dnf, &phenotype2);
        assert!((f1_2 - 1.0).abs() < 1e-6);
        let f1_3 = scorer.fitness(&dnf, &phenotype3);
        assert!((f1_3 - (2.0*precision_3*recall_3/(precision_3+recall_3))).abs() < 1e-6);

        // ACCURACY
        scorer.change_metric(ScoreMetric::Accuracy);
        let accuracy_1 = scorer.fitness(&dnf, &phenotype1);
        assert!((accuracy_1 - 0.5).abs() < 1e-6);
        let accuracy_2 = scorer.fitness(&dnf, &phenotype2);
        assert!((accuracy_2 - 1.0).abs() < 1e-6);
        let accuracy_3 = scorer.fitness(&dnf, &phenotype3);
        assert!((accuracy_3 - (3.0/4.0)).abs() < 1e-6);

    }

    #[test]
    fn test_dnfscorer_with_real_satisfaction_checker(){

    }
    
    #[test]
    fn test_display_solution() {
        let t1: TermId = "GO:0051146".parse().unwrap();

        let conj = Conjunction {
            term_observations: vec![
                TermObservation::new(t1, false),
            ],
            tissue_expressions: vec![
                TissueExpression::new("Liver".to_string(), DgeState::Up),
            ],
        };

        let sol = Solution::new(conj, 0.95);

        assert_eq!(
            format!("{}", sol),
            "<Formula: (GO:0051146 AND UP(Liver)); Score: 0.95>"
        );
    }



}


