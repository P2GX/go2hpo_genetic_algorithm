use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::phenotype2genes;
use ontolius::TermId;

fn main() {
    let phenotype2genes = phenotype2genes();

    // Filter terms with 20-80 genes (promising range based on past results)
    let mut candidates: Vec<(TermId, usize)> = phenotype2genes
        .iter()
        .filter(|(_hpo, genes)| genes.len() >= 20 && genes.len() <= 80)
        .map(|(hpo, genes)| (hpo.clone(), genes.len()))
        .collect();

    // Sort by gene count
    candidates.sort_by_key(|&(_, count)| count);

    println!("Top 20 promising HPO terms (20-80 genes):");
    for (i, (hpo, count)) in candidates.iter().take(20).enumerate() {
        println!("{}. {}: {} genes", i + 1, hpo, count);
    }
}
