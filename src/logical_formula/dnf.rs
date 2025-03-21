use super::conjunctions::{self, Conjunction};
use bitvec::prelude::*;

pub trait DNF {
    type SelectionType;

    fn get_active_conjunctions(&self) -> Vec<&Conjunction>;

    fn activate_conjunction(&mut self, conjunction: Self::SelectionType) -> Result<(), String>;

    fn len(&self) -> usize;
}

pub struct DNFVec {
    conjunctions: Vec<Conjunction>,
}

impl DNF for DNFVec {
    type SelectionType = Conjunction; //the istance itself

    fn activate_conjunction(&mut self, conjunction: Self::SelectionType) -> Result<(), String> {
        self.conjunctions.push(conjunction);
        Ok(())
    }

    fn get_active_conjunctions(&self) -> Vec<&Conjunction> {
        self.conjunctions.iter().collect()
    }

    fn len(&self) -> usize {
        self.conjunctions.len()
    }
}

impl DNFVec {
    pub fn new() -> Self {
        let conjs: Vec<Conjunction> = Vec::new();
        Self {
            conjunctions: conjs,
        }
    }
    pub fn get_mut_active_conjunctions(&mut self) -> &mut Vec<Conjunction>{
        return &mut self.conjunctions;
    }
}

pub struct DNFBitmask<'a> {
    conjunction_mask: BitVec,
    conjunctions: &'a [Conjunction], //Precomputed conjunctions
}

impl<'a> DNF for DNFBitmask<'a> {
    type SelectionType = usize; //the index

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

    fn len(&self) -> usize {
        self.conjunction_mask.count_ones()
    }
}

impl<'a> DNFBitmask<'a> {
    pub fn new(conjunctions: &'a [Conjunction]) -> Self {
        Self {
            conjunction_mask: bitvec![0; conjunctions.len()], // Initialize bitmask (all `false`)
            conjunctions,
        }
    }
}

impl DNFBitmask<'_>{
    /// Returns the total number of conjunctions (both active and inactive).
    /// This is different from len(), which only counts active ones.
    pub fn total_conjunctions_count(&self) -> usize {
        self.conjunctions.len()  // Total number of conjunctions (active + inactive)
    }

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
    use super::*;

    #[test]
    fn test_dnf_vec() {
        let conjunction1 = Conjunction { term_observations: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![] };

        let mut dnf_vec = DNFVec::new();
        dnf_vec.activate_conjunction(conjunction1.clone()).unwrap();
        dnf_vec.activate_conjunction(conjunction2.clone()).unwrap();

        let active = dnf_vec.get_active_conjunctions();
        assert_eq!(active.len(), 2, "DNFVec should contain 2 active conjunctions");
        assert_eq!(dnf_vec.len(), active.len(), "DNFVec len method should be equal to the number of active conjunctions");
    }

    #[test]
    fn test_dnf_bitmask_valid_index() {
        let conjunction1 = Conjunction { term_observations: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![] };

        let precomputed_conjunctions = [conjunction1, conjunction2];
        let mut dnf_bitmask = DNFBitmask::new(&precomputed_conjunctions);

        let result = dnf_bitmask.activate_conjunction(0);
        assert!(result.is_ok(), "Activation should succeed for index 0");

        let active = dnf_bitmask.get_active_conjunctions();
        assert_eq!(active.len(), 1, "Only one conjunction should be active");
    }

    #[test]
    fn test_dnf_bitmask_invalid_index() {
        let conjunction1 = Conjunction { term_observations: vec![] };
        let conjunction2 = Conjunction { term_observations: vec![] };

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
}
