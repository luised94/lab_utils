# ORC4R Screen - MSA Pipeline

Scripts for generating multiple sequence alignments and conservation
analysis of ORC subunits and control genes for publication.

## Execution Order

1. `orc4r-screen_01_protein-accessions.tsv` - input accession table (94 accessions, 14 organisms, 7 genes)
2. `orc4r-screen_02_download-blosum62.R` - download BLOSUM62 matrix from NCBI FTP (run once)
3. `orc4r-screen_03_fetch-and-align.R` - fetch sequences from UniProt/NCBI, align with DECIPHER, compute per-position conservation
4. `orc4r-screen_04_msa-visualization.R` - generate zoomed alignment plots at positions of interest

## Requirements

- R >= 4.4.3
- Bioconductor: Biostrings (>= 2.74.1), DECIPHER (>= 3.2.0)
- GitHub: YuLab-SMU/ggtree (>= 4.1.2), YuLab-SMU/ggmsa (>= 1.17.0)
- CRAN: ggplot2 (>= 4.0.2)
- vctrs >= 0.7.1 (required by ggmsa; update if ggmsa fails to load with a namespace error)

A pinned `renv.lock` is provided as `orc4r-screen_05_renv.lock` for exact
reproducibility. To restore the environment:

```bash
cp orc4r-screen_05_renv.lock renv.lock
Rscript -e 'renv::restore()'
```

## Output

All output is written to `~/data/protein_files/`:

- `stable_<GENE>.fasta` - raw fetched sequences (7 files)
- `BLOSUM62.rds` - substitution matrix with provenance metadata
- `alignments/stable_<GENE>_aligned.fasta` - aligned sequences
- `alignments/stable_<GENE>_conservation.tsv` - per-position conservation relative to S. cerevisiae
- `alignments/stable_summary.tsv` - alignment summary across all genes
- `alignments/plots/stable_<GENE>_pos<N>.png` and `.svg` - zoomed MSA plots
- `alignments/plots/stable_identity_similarity.tsv` - pairwise identity and similarity

## Notes

- Run scripts from the directory containing this README and the `renv` folder.
- Update `CONTACT_EMAIL` in Script 03 before distributing.
- NCBI rate limiting: Script 03 uses a 0.6s delay with jitter. Increase `REQUEST_DELAY_SEC` if you encounter 429 errors.
- DECIPHER alignment is nondeterministic; minor differences between runs are expected. Sequence ordering in Script 04 normalizes display order.
