use ontolius::TermId;
use rand::Rng;

use super::base::{Solution, FitnessScorer, FormulaEvaluator}; //to change

use super::operators::{Selection, Crossover, Mutation, ElitesSelector};
use crate::logical_formula::FormulaGenerator;


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

    // Maybe pass the list of GO annotations and tissues as arguments in the fit
    pub fn fit(&mut self) -> Solution<T> {
        
        //TO DO: initialize population if not already initialized
        //Initialization: Score initial population
            // the solution is evaluated in the moment in which it is created


        for _ in 0..self.generations {
            // New population for next generation
            let mut evolved_population: Vec<Solution<T>> =
                Vec::with_capacity(self.population.len());

            // elitism, return the size of the population already occupied by the "survivors" of the previous generation
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
                let offspring = self.evaluator.evaluate(&new_formula, &self.phenotype);

                evolved_population.push(offspring);
            }

            self.population = evolved_population;

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
