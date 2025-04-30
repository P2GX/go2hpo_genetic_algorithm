use hpo2gene_mapper::mapper::{GenePhenotypeMemoryMapper};
use anyhow::Result;
use ontolius::TermId;
use std::path::Path;
use hpo2gene_mapper::GenePhenotypeMapping;

#[test]
fn test_memory_mapper_from_sample_file() -> Result<()> {
    let hpo1: TermId = "HP:0002188".parse().unwrap();
    
    let omim1: TermId = "OMIM:619031".parse().unwrap();
    let omim2: TermId = "OMIM:619036".parse().unwrap();


    let path = Path::new("data/hpo2gene/phenotype_to_genes.txt");
    let mapper = GenePhenotypeMemoryMapper::from_file(path).unwrap();

    let gene_to_hpo = mapper.gene_to_hpo()?;
    assert!(gene_to_hpo.contains_key("ALG14"));
    assert!(gene_to_hpo["ALG14"].contains(&hpo1));

    let hpo_to_genes = mapper.hpo_to_genes()?;
    assert!(hpo_to_genes[&hpo1].contains("ALG14"));
    assert!(hpo_to_genes[&hpo1].contains("AIFM1"));

    let gene_to_disease = mapper.gene_to_diseases()?;
    let diseases = &gene_to_disease["ALG14"];
    assert!(diseases.contains(&omim1));
    assert!(diseases.contains(&omim2));

    let hpo_name_to_ncbi = mapper.hpo_name_to_ncbi_ids()?;
    assert!(hpo_name_to_ncbi["Delayed CNS myelination"].contains("199857"));

    let ncbi_to_symbol = mapper.ncbi_to_symbol()?;
    assert_eq!(ncbi_to_symbol.get("199857"), Some(&"ALG14".to_string()));
    assert_eq!(ncbi_to_symbol.get("9131"), Some(&"AIFM1".to_string()));

    Ok(())
}

