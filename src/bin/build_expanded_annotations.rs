use std::fs;
use std::path::Path;

use anyhow::Result;
use go2hpo_genetic_algorithm::{
    utils::fixtures::gene_set_annotations::{
        expand_gene_set_annotations, gene_set_annotations, go_ontology,
    },
};

fn main() -> Result<()> {
    println!("GO annotation pre-expansion builder");
    println!("==================================");

    // Load ontology and base gene set
    let go = go_ontology();
    let base = gene_set_annotations();

    // Expand and collect stats
    println!("Expanding GO annotations with ancestors...");
    let expanded = expand_gene_set_annotations(&base, &go);

    let original_total: usize = base
        .get_gene_annotations_map()
        .values()
        .map(|ann| ann.get_term_annotations().len())
        .sum();
    let expanded_total: usize = expanded
        .get_gene_annotations_map()
        .values()
        .map(|ann| ann.get_term_annotations().len())
        .sum();
    let gene_count = expanded.get_gene_annotations_map().len();
    let expansion_factor = expanded_total as f64 / original_total.max(1) as f64;

    let mut max_expansion = 0usize;
    let mut avg_expansion = 0f64;
    let mut genes_with_expansion = 0usize;

    for (gene, ann) in expanded.get_gene_annotations_map() {
        let base_count = base
            .get_gene_annotations(gene)
            .map(|g| g.get_term_annotations().len())
            .unwrap_or(0);
        let exp_count = ann.get_term_annotations().len();
        if exp_count > base_count {
            genes_with_expansion += 1;
            max_expansion = max_expansion.max(exp_count);
            avg_expansion += exp_count as f64;
        }
    }
    if genes_with_expansion > 0 {
        avg_expansion /= genes_with_expansion as f64;
    }

    // Memory estimate (rough, matches analyzer)
    let term_id_size = std::mem::size_of::<ontolius::TermId>();
    let hashset_overhead = 32usize;
    let bytes_per_annotation = term_id_size + hashset_overhead;
    let original_mb = (original_total * bytes_per_annotation) as f64 / 1_000_000.0;
    let expanded_mb = (expanded_total * bytes_per_annotation) as f64 / 1_000_000.0;

    // Persist expanded set
    let cache_path = Path::new("cache/gene_set_annotations_expanded.bincode");
    if let Some(parent) = cache_path.parent() {
        fs::create_dir_all(parent)?;
    }
    expanded.save_bincode(cache_path)?;

    println!("\nRESULTS:");
    println!("========");
    println!("Genes: {}", gene_count);
    println!("Original annotations: {}", original_total);
    println!("Expanded annotations: {}", expanded_total);
    println!("Expansion factor: {:.2}x", expansion_factor);
    println!("Per-gene: max {}, avg {:.1}", max_expansion, avg_expansion);
    println!(
        "Memory: {:.1} MB -> {:.1} MB (Δ {:.1} MB)",
        original_mb,
        expanded_mb,
        expanded_mb - original_mb
    );
    println!("Saved expanded gene set to {}", cache_path.display());

    Ok(())
}
