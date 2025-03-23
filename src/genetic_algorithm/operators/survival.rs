use crate::genetic_algorithm::Solution;

// SELECTION OF ELITE POPULATION (SURVIVAL OF BEST SOLUTIONS)

pub trait ElitesSelector<T: Clone> {
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

impl<T: Clone> ElitesSelector<T> for ElitesByNumberSelector {
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

impl<T: Clone> ElitesSelector<T> for ElitesByThresholdSelector {
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
