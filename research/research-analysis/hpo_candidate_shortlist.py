import csv
from pathlib import Path


def main() -> None:
    input_path = Path("research/research-analysis/hpo_candidate_ranking.csv")
    output_path = Path("research/research-analysis/hpo_candidate_shortlist.csv")

    min_cov = 0.95
    min_go_score = 2.0
    min_tissue_score = 2.0

    rows = []
    with input_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_count = int(float(row["gene_count"]))
            cov_frac = float(row["cov_frac"])
            go_score = float(row["best_go_score"])
            tissue_score = float(row["tissue_score"])
            combined = float(row["combined_score"])

            if gene_count < 20 or gene_count > 200:
                continue
            if cov_frac < min_cov:
                continue
            if go_score < min_go_score:
                continue
            if tissue_score < min_tissue_score:
                continue

            row["_combined"] = combined
            row["_gene_count"] = gene_count
            rows.append(row)

    rows.sort(key=lambda r: (-r["_combined"], r["_gene_count"]))

    fieldnames = [
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
    ]

    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    print(f"Wrote {len(rows)} candidates to {output_path}")
    print("Top 10 shortlist:")
    for row in rows[:10]:
        print(
            f"{row['hpo_id']}\t{row['hpo_name']}\t"
            f"genes={row['gene_count']}\t"
            f"go={float(row['best_go_score']):.2f}\t"
            f"tissue={float(row['tissue_score']):.2f}\t"
            f"combined={float(row['combined_score']):.2f}"
        )


if __name__ == "__main__":
    main()
