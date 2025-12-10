//! Elite-preservation strategies for GA populations.
use crate::genetic_algorithm::Solution;

/// SELECTION OF ELITE POPULATION (SURVIVAL OF BEST SOLUTIONS)

pub trait ElitesSelector<T: Clone> {
    /// Copy elite solutions from the previous population into `next_population`.
    /// Returns the new length of `next_population`.
    fn pass_elites(
        &self,
        next_population: &mut Vec<Solution<T>>,
        previous_population: &Vec<Solution<T>>,
        already_sorted: bool,
    ) -> usize;
}

pub struct ElitesByNumberSelector {
    number_of_elites: usize,
}

impl ElitesByNumberSelector {
    /// Keep the top `number_of_elites` solutions (requires population > elites).
    pub fn new(number_of_elites: usize) -> Self {
        if number_of_elites == 0 {
            panic!("Number of elites must be greater than 0");
        }
        Self { number_of_elites }
    }
}


impl<T: Clone> ElitesSelector<T> for ElitesByNumberSelector {
    /// If `already_sorted` is false, sorts by descending score, then takes the top-N.
    fn pass_elites(
        &self,
        next_population: &mut Vec<Solution<T>>,
        previous_population: &Vec<Solution<T>>,
        already_sorted: bool,
    ) -> usize {
        if previous_population.len() <= self.number_of_elites {
            panic!("Population size should be bigger than elites number. Population size: {}, elites number: {}", previous_population.len(),  self.number_of_elites);
        }

        let mut sorted_population = previous_population.clone();

        if !already_sorted {
            sorted_population.sort_by(|a, b| {
                b.get_score()
                    .partial_cmp(&a.get_score())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        let elites = sorted_population
            .iter()
            .take(self.number_of_elites)
            .cloned();

        next_population.extend(elites);
        next_population.len()
    }
}

pub struct ElitesByThresholdSelector {
    elite_cutoff: f64,
    maximum_number: Option<usize>,
}

impl ElitesByThresholdSelector {
    /// Keep all solutions with score >= cutoff; optionally cap at `maximum_number`.
    pub fn new(elite_cutoff: f64, maximum_number: Option<usize>) -> Self {
        Self {
            elite_cutoff,
            maximum_number,
        }
    }
}

impl<T: Clone> ElitesSelector<T> for ElitesByThresholdSelector {
    /// Selects all with score ≥ cutoff; respects optional max cap.
    fn pass_elites(
        &self,
        next_population: &mut Vec<Solution<T>>,
        previous_population: &Vec<Solution<T>>,
        _already_sorted: bool,
    ) -> usize {
        let mut elites = previous_population.iter().cloned().filter(|sol| {
            sol.get_score() >= self.elite_cutoff as f64
        });

        if let Some(max_n) = self.maximum_number {
            next_population.extend(elites.take(max_n));
        } else {
            next_population.extend(elites);
        }
        next_population.len()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_elites_by_number_selector() {
        // Create a mock population of Solutions
        let population: Vec<Solution<&str>> = vec![
            Solution::new("A", 1.0),
            Solution::new("B", 2.0),
            Solution::new("C", 3.0),
            Solution::new("D", 4.0),
        ];

        // Create selector with top 2 elites
        let selector = ElitesByNumberSelector::new(2);

        let mut next_population: Vec<Solution<&str>> = Vec::new();

        let final_len = selector.pass_elites(&mut next_population, &population, false);

        // Check that 2 elites were passed
        assert_eq!(final_len, 2);

        // Check that the elites are the ones with the highest scores
        let scores: Vec<f64> = next_population.iter().map(|s| s.get_score()).collect();
        assert_eq!(scores, vec![4.0, 3.0]);
    }

    #[test]
    fn test_elites_by_threshold_selector() {
        let population: Vec<Solution<&str>> = vec![
            Solution::new("A", 1.0),
            Solution::new("B", 2.5),
            Solution::new("C", 3.0),
            Solution::new("D", 4.5),
        ];

        // Case 1: threshold without maximum_number
        let selector_no_limit = ElitesByThresholdSelector::new(3.0, None);
        let mut next_population: Vec<Solution<&str>> = Vec::new();
        let final_len = selector_no_limit.pass_elites(&mut next_population, &population, false);

        assert_eq!(final_len, 2); // "C" and "D" survive
        let scores: Vec<f64> = next_population.iter().map(|s| s.get_score()).collect();
        assert_eq!(scores, vec![3.0, 4.5]);

        // Case 2: threshold with maximum_number
        let selector_with_limit = ElitesByThresholdSelector::new(1.5, Some(2));
        let mut next_population: Vec<Solution<&str>> = Vec::new();
        let final_len = selector_with_limit.pass_elites(&mut next_population, &population, false);

        assert_eq!(final_len, 2); // but limited by max 2
        let scores: Vec<f64> = next_population.iter().map(|s| s.get_score()).collect();
        assert!(scores.contains(&3.0) || scores.contains(&4.5)); 
    }

}


