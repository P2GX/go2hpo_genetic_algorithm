use bitvec::ptr::Mut;
use rand::prelude::*;

use crate::{conjunctions::Conjunction, dnf::{DNFBitmask, DNF}, satisfaction_checker::{self, SatisfactionChecker}};

pub trait Selection<T> {
    fn select(&self, population: &[T]) -> T;
}

pub trait Crossover<T> {
    fn crossover(&self, parent1: &T, parent2: &T) -> T;
}

pub trait Mutation<T> {
    fn mutate(&self, individual: &mut T);
}

pub struct ConjunctionMutation;

impl Mutation<Conjunction> for ConjunctionMutation {
    fn mutate(&self, individual: &mut Conjunction) {
        let mut rng = rand::rng();
        let rnd_num = rng.random_range(0..=7);
        match rnd_num{
            0 => self.mutate_with_parent_term(individual),
            1 => self.mutate_with_child_term(individual),
            2 => self.delete_term(individual),
            3 => self.add_random_term(individual),
            4 => self.toggle_term_status(individual),
            5 => self.delete_gene_expression_term(individual),
            6 => self.add_gene_expression_term(individual),
            7 => self.mutate_with_child_term(individual),
            _ => todo!()
        }
    }
}

impl ConjunctionMutation{
    //Exchange an HPO term with a parent
    pub fn mutate_with_parent_term(&self, individual: &mut Conjunction){
        todo!()
    }
    //Exchange an HPO term with one of its children
    pub fn mutate_with_child_term(&self, individual: &mut Conjunction){
        todo!()
    }
    //Delete an HPO term
    pub fn delete_term(&self, individual: &mut Conjunction){
        todo!()
    }
    //Add a random HPO term
    pub fn add_random_term(&self, individual: &mut Conjunction){
        todo!()
    }
    //Toggle the status of a gene expression term from lo to hi or vice versa (And also of GO terms with is_excluded)
    pub fn toggle_term_status(&self, individual: &mut Conjunction){
        todo!()
    }
    //Delete a gene expression term
    pub fn delete_gene_expression_term(&self, individual: &mut Conjunction){
        todo!()
    }
    //Add a gene expression term from a random tissue
    pub fn add_gene_expression_term(&self, individual: &mut Conjunction){
        todo!()
    }
}

pub struct DNFVecMutation;

impl<T> Mutation<T> for DNFVecMutation {
    fn mutate(&self, individual: &mut T) {
        todo!();
    }
}

pub struct DNFBitmaskMutation;

impl<T> Mutation<T> for DNFBitmaskMutation {
    fn mutate(&self, individual: &mut T) {
        todo!();
    }
}

//Solution, Individual, Chromosome, GeneticIndividual, FormulaWrapper
// It's just a wrapper over the formula (DNF or Conjunction) to avoid computing the same fitness function more than once
//make it implement clone
#[derive(Debug, Clone)]
pub struct Solution<T: Clone> {
    pub formula: T,
    score: Option<f64>,
}

impl<T: Clone> Solution<T> {
    pub fn get_or_score(&mut self, scorer: &dyn FitnessScorer<T>) -> f64 {
        match self.score {
            Some(val) => return val,
            None => {
                self.score = Some(scorer.fitness(&self.formula));
                return self.score.unwrap();
            }
        }
    }

    pub fn get(&self) -> Option<f64>{
        return self.score.clone();
    }
}

pub trait FitnessScorer <T>{
    fn fitness(&self, formula: &T) -> f64;
}

pub struct DNFScorer<C: SatisfactionChecker> {
    checker: C,
}

impl<C: SatisfactionChecker, T: DNF> FitnessScorer<T> for DNFScorer<C> {
    fn fitness(&self, formula: &T) -> f64 {
        1.0
    }
}

pub struct ConjunctionScorer<C: SatisfactionChecker> {
    checker: C,
}

impl<C: SatisfactionChecker> FitnessScorer<Conjunction> for ConjunctionScorer<C> {
    fn fitness(&self, formula: &Conjunction) -> f64 {
        1.0
    }
}
 
pub trait ElitesSelector<T: Clone> {
    fn pass_elites(&self, next_population: &mut Vec<Solution<T>>, previous_population: &Vec<Solution<T>>, already_sorted: bool) -> usize;
}

pub struct ElitesByNumberSelector{
    number_of_elites: usize,
}

impl<T: Clone> ElitesSelector<T> for ElitesByNumberSelector{
    fn pass_elites(&self, next_population: &mut Vec<Solution<T>>, previous_population: &Vec<Solution<T>>, already_sorted: bool) -> usize {
        if previous_population.len() <= self.number_of_elites {
            panic!("Population size should be bigger than elites number. Population size: {}, elites number: {}", previous_population.len(),  self.number_of_elites);
        }

        let mut sorted_population = previous_population.clone();

        if !already_sorted{
            sorted_population.sort_by(|a, b| {
                b.get()
                .expect("The score should already be evaluated")
                .partial_cmp(&a.get().expect("The score should already be evaluated"))
                .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
    
        let elites = sorted_population.iter().take(self.number_of_elites).cloned();

        next_population.extend(elites);
        next_population.len()
    }
}

pub struct ElitesByThresholdSelector{
    elite_cutoff: f64,
    maximum_number: Option<usize>,
}

impl<T: Clone> ElitesSelector<T> for ElitesByThresholdSelector{
    fn pass_elites(&self, next_population: &mut Vec<Solution<T>>, previous_population: &Vec<Solution<T>>, _already_sorted: bool) -> usize{
        let mut elites = previous_population.iter()
        .cloned()
        .filter(|sol| sol.get().expect("The score should already be evaluated") >= self.elite_cutoff as f64);

        if let Some(max_n) = self.maximum_number{
            next_population.extend(elites.take(max_n));
        }else{
            next_population.extend(elites);
        }
        next_population.len()
    }
}



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
            let mut evolved_population: Vec<Solution<T>> = Vec::with_capacity(self.population.len());

            // elitism, return the size of the population already occupied by the "survivors" pof the previous generation
            let number_of_elites = self.elites_selector.pass_elites(&mut evolved_population, &self.population, false);

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
        let best = self.population.iter().max_by(|a, b| {
            a.get()
                .unwrap()
                .partial_cmp(&b.get().unwrap())
                .unwrap_or(std::cmp::Ordering::Equal)
        }).unwrap().clone();
        
        // OR maybe I could order them and return the best n

        best
        
    }
}
