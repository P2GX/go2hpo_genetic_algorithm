use super::base::{Solution, FitnessScorer}; //to change

use super::ga_operators::{Selection, Crossover, Mutation, ElitesSelector};


//GeneticAlgorithm, GAEstimator
pub struct GeneticAlgorithm<T: Clone> {
    population: Vec<Solution<T>>,
    scorer: Box<dyn FitnessScorer<T>>,
    selection: Box<dyn Selection<Solution<T>>>,
    crossover: Box<dyn Crossover<Solution<T>>>,
    mutation: Box<dyn Mutation<Solution<T>>>,
    elites_selector: Box<dyn ElitesSelector<T>>,
    mutation_rate: f64,
    generations: usize,
}

impl<T: Clone> GeneticAlgorithm<T> {
    //TO DO: Two constructors:
    //      - One in which the initial population of solution is passed
    //      - One in which only the cardinality of the set of solutions per generation is passed (initizialize_population will be called)

    pub fn initialize_population(&mut self, len: usize) -> Result<(), String> {
        todo!()
    }

    pub fn fit(&mut self) -> Solution<T> {
        //Initialization: Score initial population
        for solution in self.population.iter_mut() {
            solution.get_or_score(&*self.scorer);
        }

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
                let mut offspring = self.crossover.crossover(&parent1, &parent2);

                // Mutate with a certain probability
                if rand::random::<f64>() < self.mutation_rate {
                    self.mutation.mutate(&mut offspring);
                }

                evolved_population.push(offspring);
            }

            self.population = evolved_population;

            // Score population of current generation
            for solution in self.population.iter_mut() {
                solution.get_or_score(&*self.scorer);
            }
        }

        // Return the solution with the best score, where all the other solutions of the last generation are kept in self.population
        let best = self
            .population
            .iter()
            .max_by(|a, b| {
                a.get()
                    .unwrap()
                    .partial_cmp(&b.get().unwrap())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap()
            .clone();

        // OR maybe I could order them and return the best n

        best
    }
}
