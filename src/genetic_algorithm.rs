//Solution, Individual, Chromosome, GeneticIndividual



pub trait Selection<T>{
    fn select(&self, population: &[T]) -> T;
}

pub trait Crossover<T>{
    fn crossover(&self, parent1: &T, parent2: &T) -> T;
}

pub trait Mutation<T>{
    fn mutate(&self, individual: &mut T);
}




//GeneticAlgorithm, GAEstimator
pub struct GeneticAlgorithm<T>{
    population: Vec<T>,
    selection: Box<dyn Selection<T>>,
    crossover: Box<dyn Crossover<T>>,
    mutation: Box<dyn Mutation<T>>,
    mutation_rate: f64,
    generations: usize,    
}

impl<T> GeneticAlgorithm<T>{
    
    pub fn fit(&self) -> T{
        todo!()
    }
}