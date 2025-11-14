use crate::annotations::GeneSetAnnotations;

use super::{Conjunction, DNFBitmask, DNFVec, TermObservation, TissueExpression, DNF};
use super::DgeState;
use ontolius::TermId;
use rand::{seq::SliceRandom, Rng};
use std::collections::{HashMap, HashSet};
use bitvec::prelude::BitVec;
use crate::annotations::GeneAnnotations;

pub trait FormulaGenerator {
    type Output;

    fn generate(&mut self) -> Self::Output;
}


pub trait ConjunctionGenerator: FormulaGenerator<Output = Conjunction> {}

pub trait DNFGenerator<T>: FormulaGenerator<Output = T>
where
    T: DNF,
{}


/// Creates Conjunctions randomly that might have redundant terms
pub struct RedundantRandomConjunctionGenerator<'a, R> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
    rng: R,
}

impl<'a, R> FormulaGenerator  for RedundantRandomConjunctionGenerator<'a, R> 
where
    R: Rng,
{
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {
        let mut go_terms: Vec<TermObservation> = Vec::new();
        // Select go terms randomly
        for _ in 0..self.n_go_terms {
            go_terms.push(self.select_random_observation());
        }

        let mut tissue_annots: Vec<TissueExpression> = Vec::new();
        // Select tissues annotations randomly
        for _ in 0..self.n_tissue_terms {
            tissue_annots.push(self.select_random_tissue_annot());
        }
        Conjunction {
            term_observations: go_terms,
            tissue_expressions: tissue_annots,
        }
    }
}

impl<'a, R> ConjunctionGenerator for RedundantRandomConjunctionGenerator<'a, R> where R: Rng {}

impl<'a, R> RedundantRandomConjunctionGenerator<'a, R> 
where 
    R: Rng,
{
    //randomly generate a TermObeservation
    fn select_random_observation(&mut self) -> TermObservation {
        //randomly pick a term
        let term_id = self.go_terms[self.rng.random_range(0..self.go_terms.len())].clone();
        //randomly choose whether true or false
        let is_excluded: bool = self.rng.random();
        return TermObservation::new(term_id, is_excluded);
    }

    //randomly generate a TissueExpression
    fn select_random_tissue_annot(&mut self) -> TissueExpression {
        let term_id = self.tissue_terms[self.rng.random_range(0..self.tissue_terms.len())].clone();

        let state = DgeState::get_random_up_down_only(&mut self.rng);

        return TissueExpression::new(term_id, state);
    }
    
    pub fn new(n_go_terms: usize,
            go_terms: &'a Vec<TermId>,
            n_tissue_terms: usize,
            tissue_terms: &'a Vec<String>,
            rng: R,) -> Self{
        Self { n_go_terms, go_terms, n_tissue_terms, tissue_terms, rng }
    }
}

/// Creates Conjunctions randomly whith non-repeating terms
pub struct RandomConjunctionGenerator<'a, R> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
    rng: R,
}

impl<'a, R> FormulaGenerator for RandomConjunctionGenerator<'a, R> 
where 
    R: Rng,
{
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {

        // Shuffle the annotation vectors
        let mut shuffled_go_terms = self.go_terms.clone();
        shuffled_go_terms.shuffle(&mut self.rng);

        let mut shuffled_tissue_terms = self.tissue_terms.clone();
        shuffled_tissue_terms.shuffle(&mut self.rng);


        // Take the first n according to the value specified in the relative field
        let chosen_go_terms: Vec<TermObservation> = shuffled_go_terms
            .iter()
            .take(self.n_go_terms)
            .map(|term_id| TermObservation::new(term_id.clone(), self.rng.random()))
            .collect();

        let chosen_tissues: Vec<TissueExpression> = shuffled_tissue_terms
            .iter()
            .take(self.n_tissue_terms)
            .map(|term_id| TissueExpression::new(term_id.clone(), DgeState::get_random_up_down_only(&mut self.rng)))
            .collect();


        Conjunction {
            term_observations: chosen_go_terms,
            tissue_expressions: chosen_tissues,
        }
    }
}

impl<'a, R> ConjunctionGenerator for RandomConjunctionGenerator<'a, R> where R: Rng {}

impl<'a, R> RandomConjunctionGenerator<'a, R> 
where 
    R: Rng,
{
    pub fn new(n_go_terms: usize, go_terms: &'a Vec<TermId>, n_tissue_terms: usize, tissue_terms: &'a Vec<String>, rng: R) -> Self{
        Self {n_go_terms, go_terms, n_tissue_terms, tissue_terms, rng}
    }
}


/// Picks some terms from a list of genes
pub struct GenePickerConjunctionGenerator<'a, R> {
    rng: &'a mut R,
    prob_terms: f64,
    prob_tissues: f64,
    candidate_genes: Vec<(&'a String, &'a GeneAnnotations)>,

    // Optional phenotype filter
    target_phenotype: Option<TermId>,

    // Optional limits
    max_terms: Option<usize>,
    max_tissues: Option<usize>,
}

impl<'a, R> GenePickerConjunctionGenerator<'a, R>
where
    R: Rng,
{
    pub fn new(
        rng: &'a mut R,
        prob_terms: f64,
        prob_tissues: f64,
        gene_set: &'a GeneSetAnnotations,
        target_phenotype: Option<TermId>,
        max_terms: Option<usize>,
        max_tissues: Option<usize>,
    ) -> Self {

        let annotations_map = gene_set.get_gene_annotations_map();
        let candidate_genes: Vec<_> = match &target_phenotype {
            Some(phenotype) => annotations_map
                .iter()
                .filter(|(_, ann)| ann.contains_phenotype(phenotype))
                .collect(),
            None => annotations_map.iter().collect(),
        };

        Self {
            rng,
            prob_terms,
            prob_tissues,
            candidate_genes,
            target_phenotype,
            max_terms,
            max_tissues,
        }
    }
}



impl<'a, R> FormulaGenerator for GenePickerConjunctionGenerator<'a, R>
where
    R: Rng,
{
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {
        // let annotations_map = self.gene_set.get_gene_annotations_map();

        // // Filter candidate genes depending on target_phenotype
        // let candidate_genes: Vec<_> = match &self.target_phenotype {
        //     Some(phenotype) => annotations_map
        //         .iter()
        //         .filter(|(_, ann)| ann.contains_phenotype(phenotype))
        //         .collect(),
        //     None => annotations_map.iter().collect(),
        // };

        assert!(
            !self.candidate_genes.is_empty(),
            "No candidate genes available for the given phenotype filter"
        );

        // Pick one at random
        let random_index = self.rng.random_range(0..self.candidate_genes.len());
        let (_gene_id, gene_annotations) = self.candidate_genes[random_index];

        // Sample conjunction from gene annotations
        gene_annotations.into_randomly_sampled_conjunction(
            self.prob_terms,
            self.prob_tissues,
            self.rng,
            self.max_terms,
            self.max_tissues,
        )
    }
}



impl<'a, R> ConjunctionGenerator for GenePickerConjunctionGenerator<'a, R>
where 
R: Rng{}


//DNFBitmask Generators

pub struct RandomDNFBistmaskGenerator<'a, R> {
    conjunction_list: &'a [Conjunction],
    rng: R,
}

impl<'a, R> FormulaGenerator for RandomDNFBistmaskGenerator<'a, R>
where
    R: Rng,
{
    type Output = DNFBitmask<'a>;

    fn generate(&mut self) -> Self::Output {
        let mut mask = BitVec::repeat(false, self.conjunction_list.len());
        for i in 0..mask.len() {
            if self.rng.random_bool(0.5) {
                mask.set(i, true);
            }
        }

        DNFBitmask::new_with_bitmask(self.conjunction_list, mask)
    }
}

impl<'a, R> RandomDNFBistmaskGenerator<'a, R>
where
    R: Rng,
{
    pub fn new(conjunction_list: &'a [Conjunction], rng: R) -> Self {
        Self { conjunction_list, rng }
    }

}

impl<'a, R> DNFGenerator<DNFBitmask<'a>> for RandomDNFBistmaskGenerator<'a, R> where R: Rng {}

//DNFVec Generators

pub struct RandomDNFVecGenerator<'a, CG, R>
where 
    CG: ConjunctionGenerator{
    conjunction_generator: &'a mut  CG,
    num_conjunctions: usize,
    rng: R,
}

impl<'a, CG, R> FormulaGenerator for RandomDNFVecGenerator<'a, CG, R>
where
    CG: ConjunctionGenerator, 
    R: Rng,
{
    type Output = DNFVec;
    
    fn generate(&mut self) -> Self::Output {
        let conjunctions: Vec<Conjunction> = (0..self.num_conjunctions).map(|_| self.conjunction_generator.generate()).collect();
        DNFVec::from_conjunctions(conjunctions)
    }
}

impl<'a, CG, R> RandomDNFVecGenerator<'a, CG, R>
where
    CG: ConjunctionGenerator, 
    R: Rng,
{
    pub fn new(conjunction_generator: &'a mut  CG, num_conjunctions: usize, rng: R) -> Self {
        Self { conjunction_generator, num_conjunctions, rng }
    }

}




//TO DO: fare dei test per testare anche RandomDNFBistmaskGenerator e Generators

#[cfg(test)]
mod tests {
    use crate::logical_formula::Formula;

    use super::*;
    use lazy_static::lazy_static;
    use rand::rng;

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
    fn test_randomly_pick() {
        // let small_test_terms:Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"].into_iter().map(|s| s.parse().unwrap()).collect();
        let mut generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(),
        };

        let observation = generator.select_random_observation();

        // println!("{:?}", observation);
        assert!(
            small_test_terms.contains(&observation.term_id),
            "Term ID is not in the list"
        );

        let tissue = generator.select_random_tissue_annot();
        
        assert!(
            small_test_tissues.contains(&tissue.term_id),
            "Tissue ID is not in the list"
        );

    }

    fn _test_generate_random_conjunction<T: ConjunctionGenerator>(conjunction_generator: &mut T){
        let conjunction: Conjunction = conjunction_generator.generate();
        println!("{:?}", conjunction);

        assert_eq!(conjunction.len(), 4, "The real number of terms differs from the expected one");
        
        for term_obs in conjunction.term_observations{
            assert!(
                small_test_terms.contains(&term_obs.term_id),
                "Term ID is not in the list"
            );
        }

        for tissue_term in conjunction.tissue_expressions{
            assert!(
                small_test_tissues.contains(&tissue_term.term_id),
                "Term ID is not in the list"
            );
        }
    }

    #[test]
    fn test_generate_random_conjunction_redundant() {
        let mut generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(), 
        };
        _test_generate_random_conjunction(&mut generator);
    }

    #[test]
    fn test_generate_random_conjunction() {
        let mut generator = RandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(), 
        };
        _test_generate_random_conjunction(&mut generator);
    }

}
