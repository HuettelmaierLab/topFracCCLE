# topFracCCLE - Identify human cell lines using RNA-seq derived sequence variations

# Cell Line Classification with CCLE SNP Data

This script provides a function to classify a query cell line (from a VCF file) by comparing its SNPs (single nucleotide polymorphisms) to those in the reference Cancer Cell Line Encyclopedia (CCLE) dataset. It returns the best-matching cell lines based on the fraction of shared mutations.

## Function

### `classify_cell_line_CCLE_top_mut_frac`

#### Purpose

Given a query VCF file containing sequence variants for a cell line, this function identifies the most similar CCLE cell lines by comparing their SNP profiles.

#### Arguments

- **vcf_file**: Path to the VCF file with sequence variations for the query cell line.
- **genome**: Genome version, either `"GRCh37"` or `"GRCh38"` (default: `"GRCh38"`).
- **lines**: An optional vector of cell line names to restrict the search. If `NULL`, all cell lines in the dataset are considered.
- **num_out**: Number of top matching cell lines to return (default: `3`).
- **minSNPs**: Minimum number of SNPs required in a cell line for it to be considered (default: `10`).
- **data_dir**: Directory containing the CCLE SNP and metadata files (default: `"data/"`).


#### Output

A data frame with the top matching cell lines, including:
- Cell line name
- RRID (Research Resource Identifier)
- Number of matching SNPs
- Total SNPs reported for the cell line
- Fraction of SNPs found

#### Example Usage

```r
source("classify_cell_line_CCLE_top_mut_frac_sparse_matrix.r")

result <- classify_cell_line_CCLE_top_mut_frac(
    vcf_file = "query.vcf",
    genome = "GRCh38",
    lines = NULL,
    num_out = 5,
    minSNPs = 15,
    data_dir = "data/"
)
print(result)
```

#### Requirements

- R packages: `prodlim`, `Matrix`
- CCLE mutation and metadata files in specified `data_dir`

#### Notes

- Supports both genome builds GRCh37 and GRCh38.

#### Citation

- A detailed description and evaluation can be found in our pre-print article <https://doi.org/10.1101/2025.08.22.671765>
