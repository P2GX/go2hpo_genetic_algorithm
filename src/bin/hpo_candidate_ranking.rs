use anyhow::{Context, Result};
use clap::Parser;
use csv::Writer;
use flate2::read::GzDecoder;
use oboannotation::{
    go::{stats::get_annotation_map, GoAnnotations, GoGafAnnotationLoader},
    io::AnnotationLoader,
};
use ontolius::TermId;
use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    gene_set_annotations_expanded as load_gene_set_annotations, go_ontology,
};
use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
};

#[derive(Parser, Debug)]
#[command(name = "hpo_candidate_ranking")]
#[command(about = "Rank HPO terms for GA candidate screening", long_about = None)]
struct Args {
    #[arg(long, default_value = "data/hpo2gene/phenotype_to_genes.txt")]
    hpo_path: PathBuf,

    #[arg(long, default_value = "data/gaf/goa_human.gaf.gz")]
    goa_path: PathBuf,

    #[arg(long, default_value = "data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")]
    gtex_path: PathBuf,

    #[arg(long, default_value = "data/map/hgnc_complete_set.txt")]
    hgnc_map_path: PathBuf,

    #[arg(long, default_value = "research/research-analysis/hpo_candidate_ranking.csv")]
    output: PathBuf,

    #[arg(long, default_value_t = 20)]
    min_genes: usize,

    #[arg(long, default_value_t = 200)]
    max_genes: usize,

    #[arg(long, default_value_t = 2)]
    min_go_support: usize,
}

#[derive(Debug, Clone)]
struct CandidateRow {
    hpo_id: String,
    hpo_name: String,
    gene_count: usize,
    go_cov: usize,
    gtex_cov: usize,
    both_cov: usize,
    cov_frac: f64,
    best_go_p: f64,
    best_go_score: f64,
    enriched_terms: usize,
    tissue_median_ratio: f64,
    tissue_score: f64,
    combined_score: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let (hpo_to_genes, hpo_names) = load_hpo_mappings(&args.hpo_path)?;
    let ensembl_to_symbol = load_hgnc_map(&args.hgnc_map_path)?;
    let gtex_specificity = load_gtex_specificity(&args.gtex_path, &ensembl_to_symbol)?;
    let go_annotations = load_go_annotations(&args.goa_path)?;
    let go_map = get_annotation_map(&go_annotations);
    let go = go_ontology();
    let gene_set_annotations = load_gene_set_annotations(&go);

    let mut rows = Vec::new();
    for (hpo_id, genes) in hpo_to_genes {
        let gene_count = genes.len();
        if gene_count < args.min_genes || gene_count > args.max_genes {
            continue;
        }

        let go_cov = genes.iter().filter(|g| go_map.contains_key(*g)).count();
        let gtex_cov = genes.iter().filter(|g| gtex_specificity.contains_key(*g)).count();
        let both_cov = genes
            .iter()
            .filter(|g| go_map.contains_key(*g) && gtex_specificity.contains_key(*g))
            .count();
        let cov_frac = both_cov as f64 / gene_count as f64;

        let (best_go_p, best_go_score, enriched_terms) =
            best_go_enrichment(&gene_set_annotations, &hpo_id, args.min_go_support);

        let (tissue_median_ratio, tissue_score) = tissue_specificity_score(&genes, &gtex_specificity);

        let combined_score = best_go_score + tissue_score;

        rows.push(CandidateRow {
            hpo_id: hpo_id.clone(),
            hpo_name: hpo_names
                .get(&hpo_id)
                .cloned()
                .unwrap_or_else(|| "Unknown".to_string()),
            gene_count,
            go_cov,
            gtex_cov,
            both_cov,
            cov_frac,
            best_go_p,
            best_go_score,
            enriched_terms,
            tissue_median_ratio,
            tissue_score,
            combined_score,
        });
    }

    rows.sort_by(|a, b| {
        b.combined_score
            .partial_cmp(&a.combined_score)
            .unwrap_or(Ordering::Equal)
    });

    write_csv(&args.output, &rows)?;

    println!(
        "Wrote {} candidate rows to {}",
        rows.len(),
        args.output.display()
    );
    println!("Top 10 candidates:");
    for row in rows.iter().take(10) {
        println!(
            "{}\t{}\tgenes={}\tgo_score={:.2}\ttissue_score={:.2}\tcombined={:.2}",
            row.hpo_id, row.hpo_name, row.gene_count, row.best_go_score, row.tissue_score, row.combined_score
        );
    }

    Ok(())
}

fn load_hpo_mappings(
    path: &PathBuf,
) -> Result<(HashMap<String, HashSet<String>>, HashMap<String, String>)> {
    let file = File::open(path).with_context(|| format!("open HPO file {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut hpo_to_genes: HashMap<String, HashSet<String>> = HashMap::new();
    let mut hpo_names: HashMap<String, String> = HashMap::new();

    let mut lines = reader.lines();
    lines.next(); // header
    for line in lines {
        let line = line?;
        let mut parts = line.split('\t');
        let hpo_id = match parts.next() {
            Some(v) => v,
            None => continue,
        };
        let hpo_name = match parts.next() {
            Some(v) => v,
            None => continue,
        };
        parts.next();
        let gene_symbol = match parts.next() {
            Some(v) => v,
            None => continue,
        };
        if gene_symbol.is_empty() {
            continue;
        }
        hpo_to_genes
            .entry(hpo_id.to_string())
            .or_default()
            .insert(gene_symbol.to_string());
        hpo_names
            .entry(hpo_id.to_string())
            .or_insert_with(|| hpo_name.to_string());
    }

    Ok((hpo_to_genes, hpo_names))
}

fn load_hgnc_map(path: &PathBuf) -> Result<HashMap<String, String>> {
    let file = File::open(path).with_context(|| format!("open HGNC map {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("HGNC map missing header"))??;
    let columns: Vec<&str> = header.split('\t').collect();
    let ensembl_index = columns
        .iter()
        .position(|&col| col == "ensembl_gene_id")
        .ok_or_else(|| anyhow::anyhow!("HGNC header missing ensembl_gene_id"))?;
    let symbol_index = columns
        .iter()
        .position(|&col| col == "symbol")
        .ok_or_else(|| anyhow::anyhow!("HGNC header missing symbol"))?;

    let mut map = HashMap::new();
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= ensembl_index || fields.len() <= symbol_index {
            continue;
        }
        let ensembl = fields[ensembl_index].trim();
        let symbol = fields[symbol_index].trim();
        if !ensembl.is_empty() && !symbol.is_empty() {
            map.insert(ensembl.to_string(), symbol.to_string());
        }
    }
    Ok(map)
}

fn load_gtex_specificity(
    gtex_path: &PathBuf,
    ensembl_to_symbol: &HashMap<String, String>,
) -> Result<HashMap<String, f64>> {
    let file =
        File::open(gtex_path).with_context(|| format!("open GTEx file {}", gtex_path.display()))?;
    let gz = GzDecoder::new(file);
    let reader = BufReader::new(gz);
    let mut lines = reader.lines();

    // GCT format: first two lines are metadata.
    lines.next();
    lines.next();

    let mut specificity = HashMap::new();
    for line in lines {
        let line = line?;
        let mut parts = line.split('\t');
        let ensembl_id = match parts.next() {
            Some(v) => v,
            None => continue,
        };
        let _description = parts.next();

        if ensembl_id == "Name" {
            continue;
        }

        let base_id = ensembl_id.split('.').next().unwrap_or(ensembl_id);
        let symbol = match ensembl_to_symbol.get(base_id) {
            Some(v) => v,
            None => continue,
        };

        let mut values = Vec::new();
        for v in parts {
            if let Ok(val) = v.parse::<f64>() {
                values.push(val);
            }
        }
        if values.is_empty() {
            continue;
        }
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let median = values[values.len() / 2];
        let max = *values.last().unwrap_or(&0.0);
        let ratio = if median > 0.0 {
            max / median
        } else {
            max / 1e-6
        };

        specificity.insert(symbol.to_string(), ratio);
    }

    Ok(specificity)
}

fn load_go_annotations(goa_path: &PathBuf) -> Result<GoAnnotations> {
    let file =
        File::open(goa_path).with_context(|| format!("open GOA {}", goa_path.display()))?;
    let reader: Box<dyn BufRead> = Box::new(BufReader::new(GzDecoder::new(file)));
    let loader = GoGafAnnotationLoader;
    loader.load_from_buf_read(reader).context("load GO annotations")
}

fn best_go_enrichment(
    gene_set_annotations: &GeneSetAnnotations,
    hpo_id: &str,
    min_support: usize,
) -> (f64, f64, usize) {
    let hpo_term: TermId = match hpo_id.parse() {
        Ok(term) => term,
        Err(_) => return (1.0, 0.0, 0),
    };

    let gene_map = gene_set_annotations.get_gene_annotations_map();
    let total_genes = gene_map.len();
    if total_genes == 0 {
        return (1.0, 0.0, 0);
    }

    let mut pos_total = 0usize;
    let mut total_counts: HashMap<TermId, usize> = HashMap::new();
    let mut pos_counts: HashMap<TermId, usize> = HashMap::new();

    for ann in gene_map.values() {
        let is_pos = ann.contains_phenotype(&hpo_term);
        if is_pos {
            pos_total += 1;
        }

        for term in ann.get_term_annotations() {
            *total_counts.entry(term.clone()).or_default() += 1;
            if is_pos {
                *pos_counts.entry(term.clone()).or_default() += 1;
            }
        }
    }

    if pos_total == 0 {
        return (1.0, 0.0, 0);
    }

    let mut best_p = 1.0;
    let mut best_fold_score = 0.0;
    let mut enriched_terms = 0;
    for (term, total_with_term) in total_counts {
        let pos_with_term = *pos_counts.get(&term).unwrap_or(&0);
        if pos_with_term < min_support {
            continue;
        }

        let neg_total = total_genes.saturating_sub(pos_total);
        let neg_with_term = total_with_term.saturating_sub(pos_with_term);
        let pos_rate = (pos_with_term + 1) as f64 / (pos_total + 2) as f64;
        let neg_rate = (neg_with_term + 1) as f64 / (neg_total + 2) as f64;
        let fold_enrichment = pos_rate / neg_rate;
        let fold_score = fold_enrichment.log2();
        if fold_score > best_fold_score {
            best_fold_score = fold_score;
        }

        let successes = total_with_term as u64;
        let failures = (total_genes.saturating_sub(total_with_term)) as u64;
        let draws = pos_total as u64;

        let mut p_value = if successes == 0 || draws == 0 {
            1.0
        } else {
            Hypergeometric::new(successes, failures, draws)
                .ok()
                .map(|hg| 1.0 - hg.cdf((pos_with_term as u64).saturating_sub(1)))
                .unwrap_or(1.0)
        };
        if p_value < 0.0 {
            p_value = 0.0;
        } else if p_value > 1.0 {
            p_value = 1.0;
        }
        if p_value == 0.0 {
            p_value = 1e-300;
        }

        if p_value <= 0.05 {
            enriched_terms += 1;
        }
        if p_value < best_p {
            best_p = p_value;
        }
    }

    (best_p, best_fold_score, enriched_terms)
}

fn tissue_specificity_score(
    genes: &HashSet<String>,
    gtex_specificity: &HashMap<String, f64>,
) -> (f64, f64) {
    let mut values: Vec<f64> = genes
        .iter()
        .filter_map(|g| gtex_specificity.get(g).cloned())
        .collect();
    if values.is_empty() {
        return (0.0, 0.0);
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let median = values[values.len() / 2];
    let score = if median > 0.0 { median.log2() } else { 0.0 };
    (median, score)
}

fn write_csv(path: &PathBuf, rows: &[CandidateRow]) -> Result<()> {
    let mut writer = Writer::from_path(path)?;
    writer.write_record([
        "hpo_id",
        "hpo_name",
        "gene_count",
        "go_cov",
        "gtex_cov",
        "both_cov",
        "cov_frac",
        "best_go_p",
        "best_go_score",
        "enriched_terms",
        "tissue_median_ratio",
        "tissue_score",
        "combined_score",
    ])?;
    for row in rows {
        writer.serialize((
            &row.hpo_id,
            &row.hpo_name,
            row.gene_count,
            row.go_cov,
            row.gtex_cov,
            row.both_cov,
            row.cov_frac,
            row.best_go_p,
            row.best_go_score,
            row.enriched_terms,
            row.tissue_median_ratio,
            row.tissue_score,
            row.combined_score,
        ))?;
    }
    writer.flush()?;
    Ok(())
}
