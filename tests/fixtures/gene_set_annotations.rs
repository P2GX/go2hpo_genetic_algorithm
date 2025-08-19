use std::{collections::{HashMap, HashSet}, fs::File, io::{self, BufRead, BufReader}, iter::zip, path::Path};
use anyhow::bail;
use gtex_analyzer::expression_analysis::{GtexSummary, GtexSummaryLoader};
use hpo2gene_mapper::{mapper::GenePhenotypeMemoryMapper, GenePhenotypeMapping};
use ontolius::{io::OntologyLoaderBuilder, ontology::csr::MinimalCsrOntology, TermId};
use rstest::{fixture, rstest};
use go2hpo_genetic_algorithm::{annotations::{GeneId, GeneIdMapper, GeneSetAnnotations, GtexSummaryParser}, logical_formula::TissueExpression};

use oboannotation::{go::{stats::get_annotation_map, GoAnnotations, GoGafAnnotationLoader}, io::AnnotationLoader};

// TISSUE EXPRESSION DATA
fn remove_ensembl_id_version_number<T>(hashmap: HashMap<String, T>) -> HashMap<String, T>{
    let mut new_hashmap: HashMap<String, T> = HashMap::new();
    for (ensembl_id, content) in hashmap{
        let base_id = ensembl_id.split('.').next().unwrap();
        new_hashmap.insert(base_id.to_string(), content);
    }   
    new_hashmap
}

fn get_map() -> io::Result<HashMap<String, GeneId>>{
    let file = File::open("data/map/hgnc_complete_set.txt")?;
    let reader = BufReader::new(file);

    let mut map: HashMap<String, GeneId> = HashMap::new();

    let mut lines = reader.lines();

    let header = lines.next().unwrap()?; 
    let columns: Vec<&str> = header.split('\t').collect();

    let ensembl_index = columns.iter().position(|&col| col == "ensembl_gene_id").unwrap();
    let symbol_index = columns.iter().position(|&col| col == "symbol").unwrap();

    for line in lines {
        if let Ok(line) = line {
            let fields: Vec<&str> = line.split('\t').collect();
            if let (Some(ensembl), Some(symbol)) = (fields.get(ensembl_index), fields.get(symbol_index)) {
                if !ensembl.is_empty() && !symbol.is_empty() {
                    map.insert(ensembl.to_string(), symbol.to_string());
                }
            }
        }
    }
    Ok(map)
}

#[fixture]
pub fn gtex_summary_sample() -> io::Result<GtexSummary>{
    let file_path: &str = "data/gtex/GTEx_RNASeq_gene_median_tpm_HEAD.gct";

    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let summary_loader = GtexSummaryLoader::new(Some(10), None);
    let summary = summary_loader.load_summary(reader)?;
    Ok(summary)
}

#[fixture]
pub fn gene2tissue_expr_sample(gtex_summary_sample: io::Result<GtexSummary>) -> HashMap<String, HashSet<TissueExpression>>{
    let tissue_expressions = GtexSummaryParser::parse(&gtex_summary_sample.expect("It should be Ok"));
    tissue_expressions
}

#[fixture]
pub fn gtex_summary() -> io::Result<GtexSummary>{
    let file_path: &str = "data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz";

    let file = File::open(file_path)?;
    let gz = flate2::read::GzDecoder::new(file);
    let reader = BufReader::new(gz);

    let summary_loader = GtexSummaryLoader::new(None, None);
    let summary = summary_loader.load_summary(reader)?;
    Ok(summary)
}

#[fixture]
pub fn gene2tissue_expr(gtex_summary: io::Result<GtexSummary>) -> HashMap<String, HashSet<TissueExpression>>{
    let tissue_expressions = remove_ensembl_id_version_number(GtexSummaryParser::parse(&gtex_summary.expect("It should be Ok")));
    let gene_id_mapper = GeneIdMapper::new(get_map().unwrap());
    gene_id_mapper.map_keys(tissue_expressions)
}


// GO ONTOLOGY
#[fixture]
pub fn go_sample() -> MinimalCsrOntology{
    let go_path = "data/go/go.toy.json.gz";
    let reader = flate2::bufread::GzDecoder::new(BufReader::new(
        File::open(go_path).expect("The file should be in the repo"),
    ));

    let parser = OntologyLoaderBuilder::new().obographs_parser().build();
    let go: MinimalCsrOntology = parser
        .load_from_read(reader)
        .expect("The ontology file should be OK");
    go
}
 
// GO ANNOTATION DATA

fn open_for_reading<P: AsRef<Path>>(goa_path: P) -> anyhow::Result<Box<dyn BufRead>> {
    Ok(if let Some(extension) = goa_path.as_ref().extension() {
        if extension == "gz" {   
            Box::new(BufReader::new(flate2::read::GzDecoder::new(File::open(goa_path)?)))
        } else {
            Box::new(BufReader::new(File::open(goa_path)?))
        }
    } else {
        Box::new(BufReader::new(File::open(goa_path)?))
    })
}

fn load_and_get_go_annotations() -> anyhow::Result<GoAnnotations>{
    let goa_path_str: &str = "data/gaf/goa_human.gaf.gz";
    println!("processing {}", goa_path_str);
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

#[fixture]
pub fn gene2go_terms() -> HashMap<String, HashSet<TermId>>{
    let annotations = load_and_get_go_annotations().expect("Load GoAnnotations");
    get_annotation_map(&annotations)
}


// GENE2HPO_DATA
#[fixture]
pub fn gene2phenotypes() -> HashMap<String, HashSet<TermId>>{
    let path = Path::new("data/hpo2gene/phenotype_to_genes.txt");
    let mapper = GenePhenotypeMemoryMapper::from_file(path).unwrap();

    let gene_to_hpo = mapper.gene_to_hpo().expect("It should be Ok");
    gene_to_hpo
}


// GeneSetAnnotations
#[fixture]
pub fn gene_set_annotations_sample(gene2go_terms: HashMap<String, HashSet<TermId>>, gene2tissue_expr_sample:  HashMap<String, HashSet<TissueExpression>>,gene2phenotypes: HashMap<String, HashSet<TermId>>) -> GeneSetAnnotations{
    GeneSetAnnotations::from(gene2go_terms , 
                            gene2tissue_expr_sample , 
                            gene2phenotypes)
}

#[fixture]
pub fn gene_set_annotations(gene2go_terms: HashMap<String, HashSet<TermId>>, gene2tissue_expr:  HashMap<String, HashSet<TissueExpression>>,gene2phenotypes: HashMap<String, HashSet<TermId>>) -> GeneSetAnnotations{
    GeneSetAnnotations::from(gene2go_terms , 
                            gene2tissue_expr , 
                            gene2phenotypes)
}

#[rstest]
pub fn check_if_genes_are_the_same(gene2go_terms: HashMap<String, HashSet<TermId>>,
    gene2tissue_expr:  HashMap<String, HashSet<TissueExpression>>,
                                    gene2phenotypes: HashMap<String, HashSet<TermId>>)
{
    // Get key sets
    let keys_go: HashSet<_> = gene2go_terms.keys().cloned().collect();
    let keys_tissue: HashSet<_> = gene2tissue_expr.keys().cloned().collect();
    let keys_phenotypes: HashSet<_> = gene2phenotypes.keys().cloned().collect();

    // Get intersection
    let intersection: HashSet<_> = keys_go
        .intersection(&keys_tissue)
        .cloned()
        .collect::<HashSet<_>>()
        .intersection(&keys_phenotypes)
        .cloned()
        .collect();

    // Debug or assert
    dbg!(&intersection);
    assert!(
        !intersection.is_empty(),
        "No common gene keys across all three datasets."
    );
}


#[rstest]
pub fn check_gene_keys_pairwise_intersection(
    gene2go_terms: HashMap<String, HashSet<TermId>>,
    gene2tissue_expr: HashMap<String, HashSet<TissueExpression>>,
    gene2phenotypes: HashMap<String, HashSet<TermId>>,
) {
    let keys_go: HashSet<_> = gene2go_terms.keys().cloned().collect();
    let keys_tissue: HashSet<_> = gene2tissue_expr.keys().cloned().collect();
    let keys_phenotypes: HashSet<_> = gene2phenotypes.keys().cloned().collect();

    // Pairwise intersections
    let go_tissue: HashSet<_> = keys_go.intersection(&keys_tissue).cloned().collect();
    let go_pheno: HashSet<_> = keys_go.intersection(&keys_phenotypes).cloned().collect();
    let tissue_pheno: HashSet<_> = keys_tissue.intersection(&keys_phenotypes).cloned().collect();

    dbg!(&keys_go.len(), &keys_tissue.len(), &keys_phenotypes.len());
    dbg!(&go_tissue.len(), &go_pheno.len(), &tissue_pheno.len());

    // Print void intersections
    if go_tissue.is_empty() {
        println!("No common genes between GO terms and tissue expression");
    }
    if go_pheno.is_empty() {
        println!("No common genes between GO terms and phenotypes");
    }
    if tissue_pheno.is_empty() {
        println!("No common genes between tissue expression and phenotypes");
    }

    let all_three: HashSet<_> = go_tissue.intersection(&keys_phenotypes).cloned().collect();
    dbg!(&all_three.len());

    if all_three.is_empty() {
        println!("No gene is in common in all three datasets");
    }
}




// 


#[rstest]
pub fn test_gene_set_annotations(gene_set_annotations: GeneSetAnnotations){
    assert!(gene_set_annotations.len() > 0) 
}
