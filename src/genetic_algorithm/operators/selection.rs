//! Selection strategies (tournament, roulette, rank) for choosing parents.
use crate::genetic_algorithm::Solution;
use crate::logical_formula::ConjunctionGenerator;
use num::Num;
use rand::distr::uniform::SampleUniform;
use rand::{prelude::*, rng};
use std::iter::Sum;
use std::ops::SubAssign;

use ontolius::term::simple::SimpleMinimalTerm;
use ontolius::{
    ontology::{HierarchyWalks, OntologyTerms},
    TermId,
};

use crate::logical_formula::DNFVec;
use crate::{
    logical_formula::{Conjunction, TermObservation},
    logical_formula::{DNFBitmask, DNF},
};

//SELECTION OPERATOR
// TRAIT: Selection
// Implementations:
//      TournamentSelection
//      RouletteWheelSelection
//      RankSelection

pub trait Selection<T> {
    /// Pick one solution from a population according to the strategy.
    fn select(&mut self, population: &Vec<Solution<T>>) -> Solution<T>;
}

pub struct TournamentSelection<'a, R: Rng> {
    tournament_size: usize,
    rng: &'a mut R,
}

impl<'a, R> TournamentSelection<'a, R>
where
    R: Rng,
{
    /// Create a tournament selector with a given tournament size (>=1).
    pub fn new(tournament_size: usize, rng: &'a mut R) -> Self {
        if tournament_size < 1 {
            panic!(
                "tournament_size should be at least 1. Current value: {}",
                tournament_size
            );
        }
        Self {
            tournament_size,
            rng,
        }
    }
}

impl<'a, T, R> Selection<T> for TournamentSelection<'a, R>
where
    T: Clone,
    R: Rng,
{
    /// Sample `tournament_size` individuals and return the best by score.
    fn select(&mut self, population: &Vec<Solution<T>>) -> Solution<T> {
        let rand_index = self.rng.random_range(0..population.len());
        let mut best: &Solution<T> = population
            .get(rand_index)
            .expect("The index should not be out of bounds");

        for i in 1..self.tournament_size {
            let rand_index = self.rng.random_range(0..population.len());
            let opponent: &Solution<T> = population
                .get(rand_index)
                .expect("The index should not be out of bounds");

            if best < opponent {
                best = opponent;
            }
        }
        best.clone()
    }
}

pub struct RouletteWheelSelection<'a, T, R: Rng> {
    rng: &'a mut R,
    transform: Box<dyn Fn(f64) -> f64>,

    ref_population: Option<*const Vec<Solution<T>>>,
    scores: Vec<f64>,
    tot: Option<f64>,
}

impl<'a, T, R> Selection<T> for RouletteWheelSelection<'a, T, R>
where
    T: Clone,
    R: Rng,
{
    /// Fitness-proportional selection; rebuilds weights if population changed.
    fn select(&mut self, population: &Vec<Solution<T>>) -> Solution<T> {
        if !self.is_population_same(population) {
            self.initialize(population);
        }

        let index = self.select_random_index();
        population[index].clone()
    }
}

impl<'a, T, R> RouletteWheelSelection<'a, T, R>
where
    R: Rng,
{
    /// Build a roulette selector; `transform` can rescale scores (e.g., square).
    pub fn new(rng: &'a mut R, transform: Box<dyn Fn(f64) -> f64>) -> Self {
        Self {
            rng,
            transform,
            ref_population: None,
            scores: Vec::new(),
            tot: None,
        }
    }

    /// Identity-transform default.
    pub fn default(rng: &'a mut R) -> Self {
        let identity = Box::new(|x: f64| x);
        RouletteWheelSelection::new(rng, identity)
    }

    pub fn initialize(&mut self, population: &Vec<Solution<T>>) {
        self.ref_population = Some(population);
        self.scores = population
            .iter()
            .map(|sol| sol.get_score())
            .map(|score| (self.transform)(score))
            .collect();
        self.tot = Some(self.scores.iter().copied().sum());
    }

    // check if ref_population and population point to the same vector
    pub fn is_population_same(&self, population: &Vec<Solution<T>>) -> bool {
        match self.ref_population {
            Some(ptr) => std::ptr::eq(ptr, population),
            None => false,
        }
    }

    fn select_random_index(&mut self) -> usize {
        let mut threshold = self
            .rng
            .random_range(0.0..self.tot.expect("It should be Some"));

        for (index, &score) in self.scores.iter().enumerate() {
            threshold -= score;
            if threshold <= 0.0 {
                return index;
            }
        }
        return self.scores.len() - 1;
    }
}

pub struct RankSelection<'a, T, R: Rng> {
    rng: &'a mut R,

    ref_population: Option<*const Vec<Solution<T>>>,
    ranked_population: Option<Vec<Solution<T>>>,
    tot: Option<usize>,
}

impl<'a, R, T> Selection<T> for RankSelection<'a, T, R>
where
    T: Clone,
    R: Rng,
{
    /// Rank-based selection: higher ranks have higher draw probability.
    fn select<'b>(&mut self, population: &'b Vec<Solution<T>>) -> Solution<T> {
        if !self.is_population_same(population) {
            self.initialize(population);
        }

        let mut threshold = self
            .rng
            .random_range(0..self.tot.expect("Tot should be Some"));
        let ranked = self
            .ranked_population
            .as_ref()
            .expect("ranked_population should be Some");
        for (rank, sol) in ranked.iter().enumerate() {
            threshold = threshold.saturating_sub(rank);
            if threshold <= 0 {
                return sol.clone();
            }
        }
        ranked[ranked.len() - 1].clone()
    }
}

impl<'a, T, R> RankSelection<'a, T, R>
where
    R: Rng,
    T: Clone,
{
    pub fn new(rng: &'a mut R) -> Self {
        Self {
            rng: rng,
            ref_population: None,
            ranked_population: None,
            tot: None,
        }
    }

    pub fn initialize(&mut self, population: &Vec<Solution<T>>) {
        self.ref_population = Some(population);

        let mut ranked: Vec<Solution<T>> = population.iter().cloned().collect();
        ranked.sort_unstable_by(|a, b| {
            a.get_score()
                .partial_cmp(&b.get_score())
                .expect("It should be possible to compare the values")
        });
        self.ranked_population = Some(ranked);

        self.tot = Some(((population.len() - 1) * population.len()) / 2);
    }

    // check if ref_population and population point to the same vector
    pub fn is_population_same(&self, population: &Vec<Solution<T>>) -> bool {
        match self.ref_population {
            Some(ptr) => std::ptr::eq(ptr, population),
            None => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::mock::StepRng;
    use rand::rngs::SmallRng;

    use rand::distr::uniform;
    use rand::Rng;
    use std::cell::RefCell;

    use lazy_static::lazy_static;

    lazy_static! {
        static ref SMALL_POPULATION: Vec<Solution<DNFVec>> = (0..10)
            .map(|i| i as f64 / 10.0)
            .map(|fake_score| Solution::new(DNFVec::new(), fake_score))
            .collect();
        static ref SEED_LIST: Vec<usize> = vec![2, 9, 12, 15, 17, 18, 21, 30, 42];
    }

    // #[test]
    fn test_selection_solution_from_population_tournament(seed: u64) {
        // let mut rng = StepRng::new(1, 4); // deterministic
        let mut rng = SmallRng::seed_from_u64(seed);
        // let mut rng = rand::rng(); // non deterministic
        let mut selector = TournamentSelection::new(5, &mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert!(SMALL_POPULATION.contains(&selected));
    }

    // #[test]
    fn test_selection_solution_from_population_roulette(seed: u64) {
        // let mut rng = StepRng::new(1, 4); // deterministic
        let mut rng = SmallRng::seed_from_u64(seed);
        // let mut rng = rand::rng(); // non deterministic
        let mut selector = RouletteWheelSelection::default(&mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert!(SMALL_POPULATION.contains(&selected));
    }

    // #[test]
    fn test_selection_solution_from_population_rank(seed: u64) {
        // let mut rng = StepRng::new(1, 4); // deterministic
        let mut rng = SmallRng::seed_from_u64(seed);
        // let mut rng = rand::rng(); // non deterministic
        let mut selector = RankSelection::new(&mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert!(SMALL_POPULATION.contains(&selected));
    }

    // #[test]
    fn test_selection_best_tournament(seed: u64) {
        let predefined_seed = seed;
        let tournament_size = 5;

        let mut rng_pre = SmallRng::seed_from_u64(predefined_seed);
        let mut max = 0;
        for _ in (0..tournament_size) {
            let next_rand = rng_pre.random_range(0..SMALL_POPULATION.len());
            if next_rand > max {
                max = next_rand;
            }
        }

        let mut rng = SmallRng::seed_from_u64(predefined_seed);
        let mut selector = TournamentSelection::new(tournament_size, &mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert_eq!(selected.get_score(), max as f64 / 10.0)
    }

    // #[test]
    fn test_selection_best_roulette(seed: u64) {
        let predefined_seed = seed;
        let tournament_size = 5;

        let mut rng_pre = SmallRng::seed_from_u64(predefined_seed);
        let mut index = SMALL_POPULATION.len() - 1;
        let tot = SMALL_POPULATION.iter().map(|sol| sol.get_score()).sum();
        let mut threshold = rng_pre.random_range(0.0..tot);

        let scores = SMALL_POPULATION.iter().map(|sol| sol.get_score());
        for (i, score) in scores.enumerate() {
            threshold -= score;
            if threshold <= 0.0 {
                index = i;
                break;
            }
        }

        let mut rng = SmallRng::seed_from_u64(predefined_seed);
        let mut selector = RouletteWheelSelection::default(&mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert_eq!(selected, SMALL_POPULATION[index])
    }

    // #[test]
    fn test_selection_best_rank(seed: u64) {
        let predefined_seed = seed;
        let tournament_size = 5;

        let mut ranked: Vec<Solution<DNFVec>> = SMALL_POPULATION.iter().cloned().collect();
        ranked.sort_unstable_by(|a, b| {
            a.get_score()
                .partial_cmp(&b.get_score())
                .expect("It should be possible to compare the values")
        });

        let mut rng_pre = SmallRng::seed_from_u64(predefined_seed);
        let mut choice = ranked[ranked.len() - 1].clone();
        let tot = ((ranked.len() - 1) * ranked.len()) / 2;
        let mut threshold = rng_pre.random_range(0..tot);

        for (rank, sol) in ranked.iter().enumerate() {
            threshold = threshold.saturating_sub(rank);
            if threshold <= 0 {
                choice = sol.clone();
                break;
            }
        }

        let mut rng = SmallRng::seed_from_u64(predefined_seed);
        let mut selector = RankSelection::new(&mut rng);

        let selected = selector.select(&SMALL_POPULATION);
        assert_eq!(selected, choice);
    }

    // #[test]
    fn test_panics_if_tournament_size_zero(seed: u64) {
        let result = std::panic::catch_unwind(|| {
            // let mut rng = StepRng::new(0, 1);
            let mut rng = SmallRng::seed_from_u64(seed);
            TournamentSelection::new(0, &mut rng);
        });
        assert!(result.is_err());
    }

    macro_rules! seed_test {
        ($name:ident, $func:ident, $seed:expr) => {
            #[test]
            fn $name() {
                $func($seed);
            }
        };
    }

    //test_selection_solution_from_population_tournament
    seed_test!(
        test_selection_solution_from_population_tournament_9,
        test_selection_solution_from_population_tournament,
        9
    );
    seed_test!(
        test_selection_solution_from_population_tournament_12,
        test_selection_solution_from_population_tournament,
        12
    );
    seed_test!(
        test_selection_solution_from_population_tournament_15,
        test_selection_solution_from_population_tournament,
        15
    );
    seed_test!(
        test_selection_solution_from_population_tournament_18,
        test_selection_solution_from_population_tournament,
        18
    );
    seed_test!(
        test_selection_solution_from_population_tournament_30,
        test_selection_solution_from_population_tournament,
        30
    );

    //test_selection_solution_from_population_roulette
    seed_test!(
        test_selection_solution_from_population_roulette_9,
        test_selection_solution_from_population_roulette,
        9
    );
    seed_test!(
        test_selection_solution_from_population_roulette_12,
        test_selection_solution_from_population_roulette,
        12
    );
    seed_test!(
        test_selection_solution_from_population_roulette_15,
        test_selection_solution_from_population_roulette,
        15
    );
    seed_test!(
        test_selection_solution_from_population_roulette_18,
        test_selection_solution_from_population_roulette,
        18
    );
    seed_test!(
        test_selection_solution_from_population_roulette_30,
        test_selection_solution_from_population_roulette,
        30
    );

    //test_selection_solution_from_population_rank
    seed_test!(
        test_selection_solution_from_population_rank_9,
        test_selection_solution_from_population_rank,
        9
    );
    seed_test!(
        test_selection_solution_from_population_rank_12,
        test_selection_solution_from_population_rank,
        12
    );
    seed_test!(
        test_selection_solution_from_population_rank_15,
        test_selection_solution_from_population_rank,
        15
    );
    seed_test!(
        test_selection_solution_from_population_rank_18,
        test_selection_solution_from_population_rank,
        18
    );
    seed_test!(
        test_selection_solution_from_population_rank_30,
        test_selection_solution_from_population_rank,
        30
    );

    //test_selection_best_tournament
    seed_test!(
        test_selection_best_tournament_9,
        test_selection_best_tournament,
        9
    );
    seed_test!(
        test_selection_best_tournament_12,
        test_selection_best_tournament,
        12
    );
    seed_test!(
        test_selection_best_tournament_15,
        test_selection_best_tournament,
        15
    );
    seed_test!(
        test_selection_best_tournament_18,
        test_selection_best_tournament,
        18
    );
    seed_test!(
        test_selection_best_tournament_30,
        test_selection_best_tournament,
        30
    );

    //test_selection_best_roulette
    seed_test!(
        test_selection_best_roulette_9,
        test_selection_best_roulette,
        9
    );
    seed_test!(
        test_selection_best_roulette_12,
        test_selection_best_roulette,
        12
    );
    seed_test!(
        test_selection_best_roulette_15,
        test_selection_best_roulette,
        15
    );
    seed_test!(
        test_selection_best_roulette_18,
        test_selection_best_roulette,
        18
    );
    seed_test!(
        test_selection_best_roulette_30,
        test_selection_best_roulette,
        30
    );

    //test_selection_best_rank
    seed_test!(test_selection_best_rank_9, test_selection_best_rank, 9);
    seed_test!(test_selection_best_rank_12, test_selection_best_rank, 12);
    seed_test!(test_selection_best_rank_15, test_selection_best_rank, 15);
    seed_test!(test_selection_best_rank_18, test_selection_best_rank, 18);
    seed_test!(test_selection_best_rank_30, test_selection_best_rank, 30);

    //test_panics_if_tournament_size_zero
    seed_test!(
        test_panics_if_tournament_size_zero_9,
        test_panics_if_tournament_size_zero,
        9
    );
    seed_test!(
        test_panics_if_tournament_size_zero_12,
        test_panics_if_tournament_size_zero,
        12
    );
    seed_test!(
        test_panics_if_tournament_size_zero_15,
        test_panics_if_tournament_size_zero,
        15
    );
    seed_test!(
        test_panics_if_tournament_size_zero_18,
        test_panics_if_tournament_size_zero,
        18
    );
    seed_test!(
        test_panics_if_tournament_size_zero_30,
        test_panics_if_tournament_size_zero,
        30
    );
}
