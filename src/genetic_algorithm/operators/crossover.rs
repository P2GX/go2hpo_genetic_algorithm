use bitvec::vec::BitVec;
use rand::prelude::*;
use std::panic;

use crate::logical_formula::{DNFBitmask, DNFVec, DNF};

use ontolius::term::simple::SimpleMinimalTerm;
use ontolius::{
    ontology::{HierarchyWalks, OntologyTerms},
    TermId,
};


use crate::{
    logical_formula::{Conjunction, TermObservation, TermExpression, DgeState},
};




// CROSSOVER OPERATOR

pub trait Crossover<T> {
    fn crossover(&self, parent1: &T, parent2: &T) -> T;
}


pub struct ConjunctionCrossover{
    parent1_fraction: f64, //0.5 for balanced fraction of the parents contribution on the offspring 
}

impl ConjunctionCrossover{
    pub fn new() -> Self{
        Self {parent1_fraction: 0.5}
    }
    
    pub fn new_with_parent_bias(parent1_fraction: f64) -> Self{
        if parent1_fraction < 0.0 || parent1_fraction > 1.0{
            panic!("parent1_fraction should be a probability (from 0 to 1). The current value is {}", parent1_fraction);
        }
        Self{parent1_fraction}
    }
}

impl ConjunctionCrossover{
    fn field_crossover<T>(&self, parent1_annots: &Vec<T>, parent2_annots: &Vec<T>) -> Vec<T>
    where 
    T: Clone{
        let mut rng = rand::rng();

        let mut offspring_terms = Vec::new();
        
        let amount_parent1: usize = ((parent1_annots.len() as f64) * self.parent1_fraction).round() as usize;
        let amount_parent2: usize = ((parent2_annots.len() as f64) * (1.0 - self.parent1_fraction)).round() as usize;


        let terms_from_parent1 = parent1_annots.choose_multiple(&mut rng, amount_parent1);
        let terms_from_parent2 = parent2_annots.choose_multiple(&mut rng, amount_parent2);

        offspring_terms.extend(terms_from_parent1.cloned());
        offspring_terms.extend(terms_from_parent2.cloned());
        return offspring_terms;
    }
}


impl Crossover<Conjunction> for ConjunctionCrossover{

    fn crossover(&self, parent1: &Conjunction, parent2: &Conjunction) -> Conjunction {
        
        let new_term_observation = self.field_crossover(&parent1.term_observations, &parent2.term_observations);
        let new_tissue_expression = self.field_crossover(&parent1.tissue_expressions, &parent2.tissue_expressions);
        Conjunction {
            term_observations: new_term_observation,
            tissue_expressions: new_tissue_expression,
        }
    }

}

pub struct DNFVecCrossover{
    //this a field that makes the choice between the parents unbalanced towards one of the two
    parent1_fraction: f64, //0.5 for balanced feraction of the parents contribution on the offspring 
}

impl DNFVecCrossover{
    pub fn new() -> Self{
        Self{parent1_fraction: 0.5}
    }

    pub fn new_with_parent_bias(parent1_fraction: f64) -> Self{
        if parent1_fraction < 0.0 || parent1_fraction > 1.0{
            panic!("parent1_fraction should be a probability (from 0 to 1). The current value is {}", parent1_fraction);
        }
        Self{parent1_fraction}
    }
}

impl Crossover<DNFVec> for DNFVecCrossover{
    fn crossover(&self, parent1: &DNFVec, parent2: &DNFVec) -> DNFVec {
        let mut rng = rand::rng();

        let mut offspring_conjunctions: Vec<Conjunction> = Vec::new();

        let active_conjunctions1 = parent1.get_active_conjunctions();
        let active_conjunctions2 = parent2.get_active_conjunctions();
        
        let parent1_amount: usize = ((active_conjunctions1.len() as f64) * self.parent1_fraction).round() as usize;
        let parent2_amount: usize = ((active_conjunctions2.len() as f64) * (1.0 - self.parent1_fraction)).round() as usize;

        // Select random conjunctions from each parent
        let from_parent1 = active_conjunctions1
            .choose_multiple(&mut rng, parent1_amount)
            .copied();

        let from_parent2 = active_conjunctions2
            .choose_multiple(&mut rng, parent2_amount)
            .copied();

        // Merge conjunctions while avoiding duplicates
        offspring_conjunctions.extend(from_parent1.cloned());
        offspring_conjunctions.extend(from_parent2.cloned());


        DNFVec::from_conjunctions(offspring_conjunctions)
    }
}


pub struct DNFBitmaskCrossover{
    parent1_prob: f64,        
}

impl DNFBitmaskCrossover{
    pub fn new() -> Self{
        Self{parent1_prob: 0.5}
    }

    pub fn new_with_parent_bias(parent1_prob: f64) -> Self{
        if parent1_prob < 0.0 || parent1_prob > 1.0{
            panic!("parent1_prob should be a probability (from 0 to 1). The current value is {}", parent1_prob);
        }
        Self{parent1_prob}
    }
}

// I can do a map on the offspring in which for each bit I decide weather to takethe paren1 or 2 given the prob
impl<'a> Crossover<DNFBitmask<'a>> for DNFBitmaskCrossover{
    fn crossover(&self, parent1: &DNFBitmask<'a>, parent2: &DNFBitmask<'a>) -> DNFBitmask<'a> {
        let mut rng = rand::rng();

        if parent1.get_possible_conjunctions() != parent2.get_possible_conjunctions(){
            panic!("The vector of the precomputed references should be the same for parent1 and parent2");
        }
        
        let new_mask: BitVec  = parent1
        .get_conjunction_mask()
        .iter()
        .zip(parent2.get_conjunction_mask().iter())
        .map(|(bit1, bit2)| if rng.random_bool(self.parent1_prob) {*bit1} else {*bit2})
        .collect();

        DNFBitmask::new_with_bitmask(parent1.get_possible_conjunctions(), new_mask)
    }
}





#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_crossover_field_correct_range(){
        let wrong_freq_conj_1 = panic::catch_unwind(|| ConjunctionCrossover::new_with_parent_bias(-0.1));
        assert!(wrong_freq_conj_1.is_err());
        let wrong_freq_conj_2 = panic::catch_unwind(|| ConjunctionCrossover::new_with_parent_bias(1.1));
        assert!(wrong_freq_conj_2.is_err());

        let wrong_freq_dnf_vec_1 = panic::catch_unwind(|| DNFVecCrossover::new_with_parent_bias(-0.1));
        assert!(wrong_freq_dnf_vec_1.is_err());
        let wrong_freq_dnf_vec_2 = panic::catch_unwind(|| DNFVecCrossover::new_with_parent_bias(1.1));
        assert!(wrong_freq_dnf_vec_2.is_err());

        let wrong_prob_dnf_bitmask_1 = panic::catch_unwind(|| DNFBitmaskCrossover::new_with_parent_bias(-0.1));
        assert!(wrong_prob_dnf_bitmask_1.is_err());
        let wrong_prob_dnf_bitmask_2 = panic::catch_unwind(|| DNFBitmaskCrossover::new_with_parent_bias(1.1));
        assert!(wrong_prob_dnf_bitmask_2.is_err());
    }

    #[test]
    fn test_crossover_conjunction_valid_offspring() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let t4: TermId = "GO:0042692".parse().unwrap();

        let parent1 = Conjunction {
            term_observations: vec![
                TermObservation { term_id: t1, is_excluded: false },
                TermObservation { term_id: t2, is_excluded: true },
            ],
            tissue_expressions: vec![
                TermExpression{ term_id: "Heart".to_string(), state: DgeState::Up},
                TermExpression{ term_id: "Brain".to_string(), state: DgeState::Down}
            ],
        };

        let parent2 = Conjunction {
            term_observations: vec![
                TermObservation { term_id: t3, is_excluded: false },
                TermObservation { term_id: t4, is_excluded: true },
            ],
            tissue_expressions: vec![
                TermExpression{ term_id: "Liver".to_string(), state: DgeState::Up},
                TermExpression{ term_id: "Brain".to_string(), state: DgeState::Up}
            ],
        };

        let crossover = ConjunctionCrossover::new();
        let offspring = crossover.crossover(&parent1, &parent2);

        assert!(!offspring.term_observations.is_empty(), "Offspring should not be empty.");

        let n_parent1_terms = offspring.term_observations
            .iter()
            .filter(|term| parent1.term_observations.contains(term) )
            .count();
        let n_parent2_terms = offspring.term_observations
            .iter()
            .filter(|term| parent2.term_observations.contains(term) )
            .count();

        println!("Terms from the father:{} Terms from the mother:{} Terms in total:{}", n_parent1_terms, n_parent2_terms, offspring.term_observations.len());
        assert!(n_parent1_terms > 0 , "Offspring should contain terms from Parent 1.");
        assert!(n_parent2_terms > 0, "Offspring should contain terms from Parent 2.");
        assert!(n_parent1_terms + n_parent2_terms ==  offspring.term_observations.len(), "Offspring terms should all be from the two parents");
    }

    #[test]
    fn test_crossover_conjunction_offspring_not_identical_to_parents() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let t4: TermId = "GO:0042692".parse().unwrap();

        let parent1 = Conjunction {
            term_observations: vec![
                TermObservation { term_id: t1, is_excluded: false },
                TermObservation { term_id: t2, is_excluded: true },
            ],
            tissue_expressions: vec![
                TermExpression{ term_id: "Heart".to_string(), state: DgeState::Up},
                TermExpression{ term_id: "Brain".to_string(), state: DgeState::Down}
            ],
        };

        let parent2 = Conjunction {
            term_observations: vec![
                TermObservation { term_id: t3, is_excluded: false },
                TermObservation { term_id: t4, is_excluded: true },
            ],
            tissue_expressions: vec![
                TermExpression{ term_id: "Colon".to_string(), state: DgeState::Up}]
        };

        let crossover = ConjunctionCrossover::new();
        let offspring = crossover.crossover(&parent1, &parent2);

        assert!(
            offspring.term_observations != parent1.term_observations,
            "Offspring should not be identical to Parent 1."
        );
        assert!(
            offspring.term_observations != parent2.term_observations,
            "Offspring should not be identical to Parent 2."
        );
    }

    //TO DO: Write tests for DNF crossover

}
