//! DNF representations: vector-backed and bitmask-backed forms.
use crate::logical_formula::base::Formula;

use super::conjunctions::{self, Conjunction};
use std::fmt;
use bitvec::prelude::*;

/// Disjunctive normal form interface: a set of active conjunctions (OR of ANDs).
pub trait DNF: Formula{
    type SelectionType;

    /// Get active conjunctions that make up this DNF.
    fn get_active_conjunctions(&self) -> Vec<&Conjunction>;

    /// Activate a conjunction (by value or index, depending on impl).
    fn activate_conjunction(&mut self, conjunction: Self::SelectionType) -> Result<(), String>;
}


/// DNF represented as an owned vector of conjunctions.
/// Example display: `(GO:0008150 AND UP(Liver)) OR (NOT(GO:0003674) AND DOWN(Heart))`.
#[derive(Clone, Debug)]
pub struct DNFVec {
    conjunctions: Vec<Conjunction>,
}

impl PartialEq for DNFVec{
    fn eq(&self, other: &Self) -> bool {
        self.conjunctions == other.conjunctions
    }
}

impl DNF for DNFVec {
    type SelectionType = Conjunction; //the istance itself

    /// Push a conjunction onto the active list.
    fn activate_conjunction(&mut self, conjunction: Self::SelectionType) -> Result<(), String> {
        self.conjunctions.push(conjunction);
        Ok(())
    }

    fn get_active_conjunctions(&self) -> Vec<&Conjunction> {
        self.conjunctions.iter().collect()
    }
}

impl Formula for DNFVec{
    fn len(&self) -> usize {
        self.conjunctions.len()
    }
}

impl DNFVec {
    /// Empty DNF (no active conjunctions).
    pub fn new() -> Self {
        let conjs: Vec<Conjunction> = Vec::new();
        Self {
            conjunctions: conjs,
        }
    }

    /// Build from a vector of conjunctions (all considered active).
    pub fn from_conjunctions(conjunctions: Vec<Conjunction>) -> Self{
        Self { conjunctions }
    }

    /// Mutable access to the underlying conjunction list.
    pub fn get_mut_active_conjunctions(&mut self) -> &mut Vec<Conjunction>{
        return &mut self.conjunctions;
    }
}

impl fmt::Display for DNFVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let conjs = self.get_active_conjunctions();
        if conjs.is_empty() {
            return write!(f, "FALSE"); 
        }

        let parts: Vec<String> = conjs.iter().map(|c| format!("{}", c)).collect();
        write!(f, "{}", parts.join(" OR "))
    }
}

/// DNF represented by a precomputed slice of conjunctions plus a bitmask of active ones.
pub struct DNFBitmask<'a> {
    conjunctions: &'a [Conjunction], //Precomputed conjunctions
    conjunction_mask: BitVec,
}

impl PartialEq for DNFBitmask<'_>{
    fn eq(&self, other: &Self) -> bool {
        self.conjunction_mask == other.conjunction_mask && self.conjunctions == other.conjunctions
    }
}

impl<'a> DNF for DNFBitmask<'a> {
    type SelectionType = usize; //the index

    /// Set the bit at `index` to activate a conjunction.
    fn activate_conjunction(&mut self, index: usize) -> Result<(), String> {
        if index < self.conjunction_mask.len() {
            self.conjunction_mask.set(index, true);
            Ok(())
        } else {
            Err(format!(
                "Index {} is out of bounds. Max index: {}",
                index,
                self.conjunction_mask.len() - 1
            ))
        }
    }

    fn get_active_conjunctions(&self) -> Vec<&Conjunction> {
        self.conjunction_mask
            .iter()
            .enumerate()
            .filter(|(_, bit)| **bit)
            .map(|(i, _)| &self.conjunctions[i])
            .collect()
    }
}

impl<'a> Formula for DNFBitmask<'a>{
    fn len(&self) -> usize {
        self.conjunction_mask.count_ones()
    }
}

impl<'a> DNFBitmask<'a> {
    /// Create an empty (all false) bitmask over provided conjunctions.
    pub fn new(conjunctions: &'a [Conjunction]) -> Self {
        Self {
            conjunction_mask: bitvec![0; conjunctions.len()], // Initialize bitmask (all `false`)
            conjunctions,
        }
    }

    /// Create with a custom bitmask (must match conjunctions length).
    pub fn new_with_bitmask(conjunctions: &'a[Conjunction], conjunction_mask: BitVec) -> Self{
        Self { conjunctions, conjunction_mask }
    }

    /// Access the full list of possible conjunctions (active + inactive).
    pub fn get_possible_conjunctions(&self) -> &'a [Conjunction]{
        return self.conjunctions;
    }

    /// Access the raw bitmask of active conjunctions.
    pub fn get_conjunction_mask(&self) -> &BitVec{
        &self.conjunction_mask
    }

}

impl DNFBitmask<'_>{
    /// Returns the total number of conjunctions (both active and inactive).
    /// This is different from len(), which only counts active ones.
    pub fn total_conjunctions_count(&self) -> usize {
        self.conjunctions.len()  // Total number of conjunctions (active + inactive)
    }

    /// Flip activation state of a conjunction at `index`.
    pub fn toggle_conjunction(&mut self, index: usize) -> Result<(), String>{
        if index < self.conjunction_mask.len() {
            let current_value = self.conjunction_mask[index];
            self.conjunction_mask.set(index, !current_value);
            Ok(())
        } else {
            Err(format!(
                "Index {} is out of bounds. Max index: {}",
                index,
                self.conjunction_mask.len() - 1
            ))
        }
    }
}

impl<'a> fmt::Display for DNFBitmask<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let conjs = self.get_active_conjunctions();
        if conjs.is_empty() {
            return write!(f, "FALSE"); // it is empty so it is always false
        }

        let parts: Vec<String> = conjs.iter().map(|c| format!("{}", c)).collect();
        write!(f, "{}", parts.join(" OR "))
    }
}

// pub struct DNF{
//     conjunctions : Vec<Conjunction>,
// }

// impl DNF{
//     pub fn len(&self) -> usize{
//         return self.conjunctions.len();
//     }
// }

// pub trait DNFGenerator {
//     fn generate(&self) -> DNF;
// }

// pub struct RandomPrecomputedDNFGenerator<'a> {
//     n_conjunctions: usize,
//     terms: &'a Vec<Conjunction>,
// }

// impl<'a> DNFGenerator for RandomPrecomputedDNFGenerator<'a> {
//     fn generate(&self) -> DNF {
//         todo!()
//     }
// }

#[cfg(test)]
mod tests {
    use ontolius::TermId;

    use crate::logical_formula::{DgeState, TermObservation, TissueExpression};

    use super::*;

    #[test]
    fn test_dnf_vec() {
        let conjunction1 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };

        let mut dnf_vec = DNFVec::new();
        dnf_vec.activate_conjunction(conjunction1.clone()).unwrap();
        dnf_vec.activate_conjunction(conjunction2.clone()).unwrap();

        let active = dnf_vec.get_active_conjunctions();
        assert_eq!(active.len(), 2, "DNFVec should contain 2 active conjunctions");
        assert_eq!(dnf_vec.len(), active.len(), "DNFVec len method should be equal to the number of active conjunctions");
    }

    #[test]
    fn test_dnf_bitmask_valid_index() {
        let conjunction1 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };

        let precomputed_conjunctions = [conjunction1, conjunction2];
        let mut dnf_bitmask = DNFBitmask::new(&precomputed_conjunctions);

        let result = dnf_bitmask.activate_conjunction(0);
        assert!(result.is_ok(), "Activation should succeed for index 0");

        let active = dnf_bitmask.get_active_conjunctions();
        assert_eq!(active.len(), 1, "Only one conjunction should be active");
    }

    #[test]
    fn test_dnf_bitmask_invalid_index() {
        let conjunction1 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![], tissue_expressions: vec![] };

        let precomputed_conjunctions = [conjunction1, conjunction2];
        let mut dnf_bitmask = DNFBitmask::new(&precomputed_conjunctions);

        let result = dnf_bitmask.activate_conjunction(5);  // Out of bounds
        assert!(result.is_err(), "Activation should fail because the index is not valid");
        assert_eq!(
            result.unwrap_err(),
            "Index 5 is out of bounds. Max index: 1"
        );

        let active = dnf_bitmask.get_active_conjunctions();
        assert_eq!(active.len(), 0, "No conjunctions should be active");
    }

    // #[test]
    // fn test_display_dnfvec() {
    //     let t1: TermId = "GO:0051146".parse().unwrap();
    //     let t2: TermId = "GO:0051216".parse().unwrap();

    //     let conj1 = Conjunction {
    //         term_observations: vec![
    //             TermObservation::new(t1.clone(), false),
    //         ],
    //         tissue_expressions: vec![
    //             TissueExpression::new("Liver".to_string(), DgeState::Up),
    //         ],
    //     };

    //     let conj2 = Conjunction {
    //         term_observations: vec![
    //             TermObservation::new(t2.clone(), true),
    //         ],
    //         tissue_expressions: vec![
    //             TissueExpression::new("Heart".to_string(), DgeState::Down),
    //         ],
    //     };

    //     let dnf = DNFVec::from_conjunctions(vec![conj1, conj2]);

    //     assert_eq!(
    //         format!("{}", dnf),
    //         "(GO:0051146 AND UP(Liver)) OR (NOT(GO:0051216) AND DOWN(Heart))"
    //     );
    // }


    #[test]
    fn test_display_dnfvec_and_bitmask() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0051216".parse().unwrap();

        let conj1 = Conjunction {
            term_observations: vec![
                TermObservation::new(t1.clone(), false),
            ],
            tissue_expressions: vec![
                TissueExpression::new("Liver".to_string(), DgeState::Up),
            ],
        };

        let conj2 = Conjunction {
            term_observations: vec![
                TermObservation::new(t2.clone(), true),
            ],
            tissue_expressions: vec![
                TissueExpression::new("Heart".to_string(), DgeState::Down),
            ],
        };

        // Test DNFVec
        let dnf_vec = DNFVec::from_conjunctions(vec![conj1.clone(), conj2.clone()]);
        println!("DNF Vec: {}",dnf_vec);
        assert_eq!(
            format!("{}", dnf_vec),
            "(GO:0051146 AND UP(Liver)) OR (NOT(GO:0051216) AND DOWN(Heart))"
        );

        // Test DNFBitmask
        let binding = [conj1.clone(), conj2.clone()];
        let mut dnf_bitmask = DNFBitmask::new(&binding);
        dnf_bitmask.activate_conjunction(0).unwrap();
        dnf_bitmask.activate_conjunction(1).unwrap();
        println!("DNF Bitmask: {}",dnf_bitmask);
        assert_eq!(
            format!("{}", dnf_bitmask),
            "(GO:0051146 AND UP(Liver)) OR (NOT(GO:0051216) AND DOWN(Heart))"
        );
    }


}
