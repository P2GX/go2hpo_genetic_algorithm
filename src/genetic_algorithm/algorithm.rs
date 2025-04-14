use rand::Rng;

use super::base::{Solution, FitnessScorer, FormulaEvaluator}; //to change

use super::operators::{Selection, Crossover, Mutation, ElitesSelector};
use crate::logical_formula::FormulaGenerator;


//GeneticAlgorithm, GAEstimator
pub struct GeneticAlgorithm<T, R> {
    
    population: Vec<Solution<T>>,
    evaluator: FormulaEvaluator<T>,
    scorer: Box<dyn FitnessScorer<T>>, // probably this can be removed since there is FormulaEvaluator
    selection: Box<dyn Selection<T>>,
    crossover: Box<dyn Crossover<T>>,
    mutation: Box<dyn Mutation<T>>,
    elites_selector: Box<dyn ElitesSelector<T>>,
    formula_generator: Box<dyn FormulaGenerator<Output = T>>,
    mutation_rate: f64,
    generations: usize,
    rng: R,
}

impl<T: Clone, R: Rng> GeneticAlgorithm<T, R> {
    //TO DO: Two constructors:
    //      - One in which the initial population of solution is passed
    //      - One in which only the population size is passed (initizialize_population will be called)

    pub fn initialize_population(&mut self, len: usize) -> Result<(), String> {
        self.population = (0..len)
                        .map(|_| self.formula_generator.generate())
                        .map(|formula| self.evaluator.evaluate(formula))
                        .collect();
        Ok(())
    }

    pub fn fit(&mut self) -> Solution<T> {
        //Initialization: Score initial population, , there is no need for this anymore because now 
            // the solution is evaluated in the moment in which it is created
        // for solution in self.population.iter_mut() {
        //     solution.get_or_score(&*self.scorer);
        // }

        for _ in 0..self.generations {
            // New population for next generation
            let mut evolved_population: Vec<Solution<T>> =
                Vec::with_capacity(self.population.len());

            // elitism, return the size of the population already occupied by the "survivors" pof the previous generation
            let number_of_elites =
                self.elites_selector
                    .pass_elites(&mut evolved_population, &self.population, false);

            // new generation
            for _ in number_of_elites..self.population.len() {
                // Selection and crossover
                let parent1 = self.selection.select(&self.population);
                let parent2 = self.selection.select(&self.population);
                let mut new_formula = self.crossover.crossover(parent1.get_formula(), parent2.get_formula());

                // Mutate with a certain probability
                if self.rng.random::<f64>() < self.mutation_rate {
                    self.mutation.mutate(&mut new_formula);
                }

                //Evaluate the the new formula
                let offspring = self.evaluator.evaluate(new_formula);

                evolved_population.push(offspring);
            }

            self.population = evolved_population;


            // Score population of current generation, there is no need for this anymore because now 
            // the solution is evaluated in the moment in which it is created
            //
            // for solution in self.population.iter_mut() {
            //     solution.get_or_score(&*self.scorer);
            // }
        }

        // Return the solution with the best score, where all the other solutions of the last generation are kept in self.population
        let best = self
            .population
            .iter()
            .max_by(|a, b| {
                a.get_score()
                    .partial_cmp(&b.get_score())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap()
            .clone();

        // OR maybe I could order them and return the best n

        best
    }
}
