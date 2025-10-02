use std::sync::Arc;

use ontolius::TermId;
use rand::Rng;

use super::base::{Solution, FitnessScorer, FormulaEvaluator}; //to change

use super::operators::{Selection, Crossover, Mutation, ElitesSelector};
use crate::genetic_algorithm::ScoreMetric;
use crate::logical_formula::{Formula, FormulaGenerator};


//GeneticAlgorithm, GAEstimator
pub struct GeneticAlgorithm<'a, T, R> {
    population: Vec<Solution<T>>,
    evaluator: FormulaEvaluator<'a, T, TermId>,
    selection: Box<dyn Selection<T> + 'a>,
    crossover: Box<dyn Crossover<T> + 'a>,
    mutation: Box<dyn Mutation<T> + 'a>,
    elites_selector: Box<dyn ElitesSelector<T> + 'a>,
    formula_generator: Box<dyn FormulaGenerator<Output = T> + 'a>,
    mutation_rate: f64,
    generations: usize,
    rng: R,
    phenotype: TermId,
}


impl<'a, T: Clone, R: Rng> GeneticAlgorithm<'a, T, R> {

    pub fn initialize_population(&mut self, len: usize) -> Result<(), String> {
        self.population = (0..len)
                        .map(|_| self.formula_generator.generate())
                        .map(|formula| self.evaluator.evaluate(&formula, &self.phenotype))
                        .collect();
        Ok(())
    }

    // // Maybe pass the list of GO annotations and tissues as arguments in the fit
    // pub fn fit(&mut self) -> Solution<T> {
        
    //     for _ in 0..self.generations {
    //         // New population for next generation
    //         let mut evolved_population: Vec<Solution<T>> =
    //             Vec::with_capacity(self.population.len());

    //         // elitism, return the size of the population already occupied by the "survivors" of the previous generation
    //         let number_of_elites =
    //             self.elites_selector
    //                 .pass_elites(&mut evolved_population, &self.population, false);

    //         // new generation
    //         for _ in number_of_elites..self.population.len() {
    //             // Selection and crossover
    //             let parent1 = self.selection.select(&self.population);
    //             let parent2 = self.selection.select(&self.population);
    //             let mut new_formula = self.crossover.crossover(parent1.get_formula(), parent2.get_formula());

    //             // Mutate with a certain probability
    //             if self.rng.random::<f64>() < self.mutation_rate {
    //                 self.mutation.mutate(&mut new_formula);
    //             }

    //             //Evaluate the the new formula
    //             let offspring = self.evaluator.evaluate(&new_formula, &self.phenotype);

    //             evolved_population.push(offspring);
    //         }

    //         self.population = evolved_population;

    //     }

    //     // Return the solution with the best score, where all the other solutions of the last generation are kept in self.population
    //     let best = self
    //         .population
    //         .iter()
    //         .max_by(|a, b| {
    //             a.get_score()
    //                 .partial_cmp(&b.get_score())
    //                 .unwrap_or(std::cmp::Ordering::Equal)
    //         })
    //         .unwrap()
    //         .clone();

    //     best
    // }

    /// Get a reference to the current population
    pub fn get_population(&self) -> &Vec<Solution<T>> {
        &self.population
    }


}



// impl<'a, T: Clone, R: Rng> GeneticAlgorithm<'a, T, R>{
//     pub fn fit_with_history(&mut self) -> Vec<Vec<Solution<T>>> {
//         let mut history = Vec::with_capacity(self.generations);

//         for _ in 0..self.generations {
//             let mut evolved_population: Vec<Solution<T>> =
//                 Vec::with_capacity(self.population.len());

//             // Elitism
//             let number_of_elites =
//                 self.elites_selector
//                     .pass_elites(&mut evolved_population, &self.population, false);

//             // New generation
//             for _ in number_of_elites..self.population.len() {
//                 let parent1 = self.selection.select(&self.population);
//                 let parent2 = self.selection.select(&self.population);
//                 let mut new_formula =
//                     self.crossover.crossover(parent1.get_formula(), parent2.get_formula());

//                 if self.rng.random::<f64>() < self.mutation_rate {
//                     self.mutation.mutate(&mut new_formula);
//                 }

//                 let offspring = self.evaluator.evaluate(&new_formula, &self.phenotype);
//                 evolved_population.push(offspring);
//             }

//             // Save the population before moving on
//             history.push(evolved_population.clone());

//             // Update population
//             self.population = evolved_population;
//         }

//         history
//     }
// }



impl<'a, T: Clone, R: Rng> GeneticAlgorithm<'a, T, R>{
    /// Helper: evolve one generation and return the new population
    pub fn evolve_one_generation(&mut self) -> Vec<Solution<T>> {
        let mut evolved_population: Vec<Solution<T>> =
            Vec::with_capacity(self.population.len());

        // Elitism
        let number_of_elites =
            self.elites_selector
                .pass_elites(&mut evolved_population, &self.population, false);

        // New generation
        for _ in number_of_elites..self.population.len() {
            let parent1 = self.selection.select(&self.population);
            let parent2 = self.selection.select(&self.population);
            let mut new_formula =
                self.crossover.crossover(parent1.get_formula(), parent2.get_formula());

            if self.rng.random::<f64>() < self.mutation_rate {
                self.mutation.mutate(&mut new_formula);
            }

            let offspring = self.evaluator.evaluate(&new_formula, &self.phenotype);
            evolved_population.push(offspring);
        }

        evolved_population
    }

    /// Original version: return only the best final solution
    pub fn fit(&mut self) -> Solution<T> {
        for _ in 0..self.generations {
            let evolved_population = self.evolve_one_generation();
            self.population = evolved_population;
        }

        // Return best of last generation
        self.population
            .iter()
            .max_by(|a, b| {
                a.get_score()
                    .partial_cmp(&b.get_score())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap()
            .clone()
    }

    /// New version: return all populations across generations
    pub fn fit_with_history(&mut self) -> Vec<Vec<Solution<T>>> {
        let mut history = Vec::with_capacity(self.generations);

        for _ in 0..self.generations {
            let evolved_population = self.evolve_one_generation();

            // Save snapshot of this generation
            history.push(evolved_population.clone());

            // Update current population
            self.population = evolved_population;
        }

        history
    }


}



impl<'a, T: Clone + Formula, R: Rng> GeneticAlgorithm<'a, T, R> {
    /// New version: return statistics (min, avg, max score) for each generation
    pub fn fit_with_stats_history(&mut self) -> Vec<(f64, f64, f64, usize, f64, usize, f64, f64)> {
        let mut stats_history = Vec::with_capacity(self.generations);

        for ith_gen in 0..self.generations {
            println!("Generation number: {}", ith_gen);
            let evolved_population = self.evolve_one_generation();

            // Compute scores
            let scores: Vec<f64> = evolved_population
                .iter()
                .map(|s| s.get_score())
                .collect();

            let min_score = scores
                .iter()
                .cloned()
                .fold(f64::INFINITY, f64::min);
            let max_score = scores
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            let avg_score = scores.iter().sum::<f64>() / scores.len() as f64;

            // Compute Precision Score and Recall Score
            let formulas: Vec<&T> = evolved_population.
                iter()
                .map(|s| s.get_formula())
                .collect();

            // Formula with the highest score (first best)
            let best_formula: &T = evolved_population
                .iter()
                .max_by(|a, b| a.get_score().partial_cmp(&b.get_score()).unwrap())
                .map(|s| s.get_formula()).unwrap();

            let best_one_precision = self.evaluator.get_scorer().custom_score_fitness(best_formula, &self.phenotype, &ScoreMetric::Precision);
            let best_one_recall = self.evaluator.get_scorer().custom_score_fitness(best_formula, &self.phenotype, &ScoreMetric::Recall);

            // Compute avg length
            let (sum, min_len, max_len) = evolved_population
                .iter()
                .map(|s| s.get_formula().len())
                .fold((0usize, usize::MAX, 0usize), |(sum, min, max), len| {
                (sum + len, min.min(len), max.max(len))
            });

            let avg_len = sum as f64 / evolved_population.len() as f64;
                            

            stats_history.push((min_score, avg_score, max_score, min_len, avg_len, max_len, best_one_precision, best_one_recall));

            // Update current population
            self.population = evolved_population;
        }

        stats_history
    }
}



//Two constructors:
//      - One in which the initial population of solution is passed
//      - One in which only the population size is passed (initizialize_population will be called)

impl<'a, T: Clone, R: Rng> GeneticAlgorithm<'a, T, R> {
    /// Constructor with an already initialized population
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_population(
        population: Vec<Solution<T>>,
        evaluator: FormulaEvaluator<'a, T, TermId>,
        selection: Box<dyn Selection<T> + 'a>,
        crossover: Box<dyn Crossover<T> + 'a>,
        mutation: Box<dyn Mutation<T> + 'a>,
        elites_selector: Box<dyn ElitesSelector<T> + 'a>,
        formula_generator: Box<dyn FormulaGenerator<Output = T> + 'a>,
        mutation_rate: f64,
        generations: usize,
        rng: R,
        phenotype: TermId,
    ) -> Self {
        Self {
            population,
            evaluator,
            selection,
            crossover,
            mutation,
            elites_selector,
            formula_generator,
            mutation_rate,
            generations,
            rng,
            phenotype,
        }
    }
    /// Constructor with only population size (population will be generated)
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_size(
        population_size: usize,
        evaluator: FormulaEvaluator<'a, T, TermId>,
        selection: Box<dyn Selection<T> + 'a>,
        crossover: Box<dyn Crossover<T> + 'a>,
        mutation: Box<dyn Mutation<T> + 'a>,
        elites_selector: Box<dyn ElitesSelector<T> + 'a>,
        mut formula_generator: Box<dyn FormulaGenerator<Output = T> + 'a>,
        mutation_rate: f64,
        generations: usize,
        rng: R,
        phenotype: TermId,
    ) -> Self {
        let population: Vec<Solution<T>> = (0..population_size)
            .map(|_| formula_generator.generate())
            .map(|formula| evaluator.evaluate(&formula, &phenotype))
            .collect();

        Self {
            population,
            evaluator,
            selection,
            crossover,
            mutation,
            elites_selector,
            formula_generator,
            mutation_rate,
            generations,
            rng,
            phenotype,
        }
    }
}
