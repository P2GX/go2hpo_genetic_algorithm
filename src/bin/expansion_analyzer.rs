use ontolius::{
    ontology::{csr::MinimalCsrOntology, HierarchyWalks, OntologyTerms},
    TermId,
};
use std::collections::{HashMap, HashSet};

// Import the same utilities used by the main application
use go2hpo_genetic_algorithm::{
    annotations::GeneSetAnnotations, logical_formula::TissueExpression,
};

// Copy the data loading functions from the fixtures
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Result};
use flate2::bufread::GzDecoder;
use oboannotation::{
    go::{stats::get_annotation_map, GoAnnotations, GoGafAnnotationLoader},
    io::AnnotationLoader,
};

fn open_for_reading(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(flate2::read::GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn load_go_annotations() -> Result<GoAnnotations> {
    let goa_path_str: &str = "data/gaf/goa_human.gaf.gz";
    println!("Loading GO annotations from {}", goa_path_str);
    let reader: Box<dyn BufRead> = open_for_reading(goa_path_str)?;

    let loader = GoGafAnnotationLoader;
    let annotations = match loader.load_from_buf_read(reader) {
        Ok(annotations) => {
            println!(
                "Loaded {:?} annotations. Skipped {:?} negated annotations",
                annotations.annotations.len(),
                annotations.negated_annotation_count
            );
            annotations
        }
        Err(e) => {
            bail!("Could not load GOA: {e}")
        }
    };
    Ok(annotations)
}

fn load_go_ontology() -> Result<MinimalCsrOntology> {
    // Use the filtered full GO ontology for realistic expansion stats.
    let go_path = "data/go/go-basic.filtered.json.gz";
    let reader = flate2::bufread::GzDecoder::new(BufReader::new(
        File::open(go_path).expect("The GO ontology file should exist"),
    ));

    let parser = ontolius::io::OntologyLoaderBuilder::new()
        .obographs_parser()
        .build();
    let go: MinimalCsrOntology = parser
        .load_from_read(reader)
        .expect("The ontology file should be OK");

    println!(
        "Loaded GO ontology with {} terms",
        go.iter_term_ids().count()
    );
    Ok(go)
}

fn load_gene_set_annotations() -> Result<GeneSetAnnotations> {
    use std::path::Path;

    let cache_path = Path::new("cache/gene_set_annotations.bincode");

    // Try to load from cache first
    if let Ok(gs) = GeneSetAnnotations::load_bincode(cache_path) {
        println!("Loaded gene set from cache with {} genes", gs.len());
        return Ok(gs);
    }

    println!("Cache not found, loading fresh data...");

    // Load GO annotations
    let go_annotations = load_go_annotations()?;
    let gene2go: HashMap<String, HashSet<TermId>> = get_annotation_map(&go_annotations);

    // For analysis purposes, create dummy tissue and phenotype data
    // so genes with GO annotations are included
    let mut gene2tissues: HashMap<String, HashSet<TissueExpression>> = HashMap::new();
    let mut gene2hpo: HashMap<String, HashSet<TermId>> = HashMap::new();

    // Add dummy data for genes that have GO annotations
    for gene_id in gene2go.keys() {
        gene2tissues.insert(gene_id.clone(), HashSet::new());
        gene2hpo.insert(gene_id.clone(), HashSet::new());
    }

    let gene_set = GeneSetAnnotations::from(gene2go, gene2tissues, gene2hpo);
    println!("Loaded fresh gene set with {} genes", gene_set.len());

    Ok(gene_set)
}

#[derive(Debug)]
struct ExpansionStats {
    gene_count: usize,
    original_total_annotations: usize,
    expanded_total_annotations: usize,
    expansion_factor: f64,
    max_expansion_per_gene: usize,
    avg_expansion_per_gene: f64,
    genes_with_expansion: usize,
}

fn analyze_expansion_impact(
    gene_set: &GeneSetAnnotations,
    go_ontology: &MinimalCsrOntology,
) -> Result<ExpansionStats> {
    println!("\nBuilding ancestor cache...");
    // Precompute ancestor relationships for efficiency
    let mut ancestor_cache: HashMap<TermId, HashSet<TermId>> = HashMap::new();

    for term in go_ontology.iter_term_ids() {
        let ancestors: HashSet<TermId> = go_ontology.iter_ancestor_ids(term).cloned().collect();
        ancestor_cache.insert(term.clone(), ancestors);
    }

    println!("Analyzing gene annotation expansion...");

    let mut original_total_annotations = 0;
    let mut expanded_total_annotations = 0;
    let mut max_expansion_per_gene = 0;
    let mut avg_expansion_per_gene = 0.0;
    let mut genes_with_expansion = 0;

    let gene_annotations = gene_set.get_gene_annotations_map();

    for (gene_id, gene_annot) in gene_annotations.iter() {
        let direct_terms = gene_annot.get_term_annotations();
        let original_count = direct_terms.len();
        original_total_annotations += original_count;

        // Calculate expanded annotations
        let mut expanded_terms = HashSet::new();

        for term in direct_terms {
            // Add the direct term
            expanded_terms.insert(term.clone());

            // Add all ancestors
            if let Some(ancestors) = ancestor_cache.get(term) {
                expanded_terms.extend(ancestors.clone());
            }
        }

        let expanded_count = expanded_terms.len();
        expanded_total_annotations += expanded_count;

        if expanded_count > original_count {
            genes_with_expansion += 1;
            max_expansion_per_gene = max_expansion_per_gene.max(expanded_count);
            avg_expansion_per_gene += expanded_count as f64;
        }
    }

    if genes_with_expansion > 0 {
        avg_expansion_per_gene /= genes_with_expansion as f64;
    }

    let expansion_factor = if original_total_annotations > 0 {
        expanded_total_annotations as f64 / original_total_annotations as f64
    } else {
        1.0
    };

    Ok(ExpansionStats {
        gene_count: gene_annotations.len(),
        original_total_annotations,
        expanded_total_annotations,
        expansion_factor,
        max_expansion_per_gene,
        avg_expansion_per_gene,
        genes_with_expansion,
    })
}

fn estimate_memory_usage(stats: &ExpansionStats) -> (f64, f64) {
    // Estimate memory usage per annotation
    // TermId is ~24 bytes (String overhead), HashSet entry has overhead
    let term_id_size = std::mem::size_of::<TermId>(); // ~24 bytes
    let hashset_entry_overhead = 32; // conservative estimate
    let bytes_per_annotation = term_id_size + hashset_entry_overhead;

    let original_memory_mb =
        (stats.original_total_annotations * bytes_per_annotation) as f64 / 1_000_000.0;
    let expanded_memory_mb =
        (stats.expanded_total_annotations * bytes_per_annotation) as f64 / 1_000_000.0;

    (original_memory_mb, expanded_memory_mb)
}

fn main() -> Result<()> {
    println!("GO Annotation Expansion Impact Analysis");
    println!("=======================================");

    // Load data
    let gene_set = load_gene_set_annotations()?;
    let go_ontology = load_go_ontology()?;

    // Run analysis
    let stats = analyze_expansion_impact(&gene_set, &go_ontology)?;

    // Calculate memory estimates
    let (original_mb, expanded_mb) = estimate_memory_usage(&stats);

    // Print results
    println!("\nRESULTS:");
    println!("========");
    println!("Genes analyzed: {}", stats.gene_count);
    println!("Genes with GO annotations: {}", stats.genes_with_expansion);
    println!("");
    println!(
        "Original direct annotations: {}",
        stats.original_total_annotations
    );
    println!(
        "Expanded annotations (with ancestors): {}",
        stats.expanded_total_annotations
    );
    println!("Expansion factor: {:.2}x", stats.expansion_factor);
    println!("");
    println!("Per-gene statistics:");
    println!(
        "  Max expansion: {} annotations",
        stats.max_expansion_per_gene
    );
    println!(
        "  Avg expansion: {:.1} annotations",
        stats.avg_expansion_per_gene
    );
    println!("");
    println!("Memory estimates:");
    println!("  Original: {:.1} MB", original_mb);
    println!("  Expanded: {:.1} MB", expanded_mb);
    println!(
        "  Increase: {:.1} MB (+{:.1}x)",
        expanded_mb - original_mb,
        stats.expansion_factor
    );

    // Analysis conclusion
    println!("\nCONCLUSION:");
    println!("===========");
    if stats.expansion_factor > 10.0 {
        println!(
            "❌ SEVERE EXPANSION: {:.1}x factor suggests pre-expansion is impractical",
            stats.expansion_factor
        );
    } else if stats.expansion_factor > 5.0 {
        println!(
            "⚠️  MODERATE EXPANSION: {:.1}x factor - considerable memory increase",
            stats.expansion_factor
        );
    } else {
        println!(
            "✅ MILD EXPANSION: {:.1}x factor - pre-expansion could be viable",
            stats.expansion_factor
        );
    }

    if expanded_mb > 100.0 {
        println!(
            "❌ Memory usage ({:.1} MB) exceeds reasonable limits",
            expanded_mb
        );
    } else if expanded_mb > 50.0 {
        println!(
            "⚠️  Memory usage ({:.1} MB) is significant but manageable",
            expanded_mb
        );
    } else {
        println!("✅ Memory usage ({:.1} MB) is acceptable", expanded_mb);
    }

    println!("\nRecommendation: Runtime caching is likely preferable to pre-expansion.");

    Ok(())
}
