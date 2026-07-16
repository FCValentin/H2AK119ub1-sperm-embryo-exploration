# Sperm-derived H2AK119ub1 is required for embryonic development in *Xenopus laevis*

**Valentin FRANCOIS--CAMPION, Florian Berger, Mami Oikawa, Maissa Goumeidane, Romain Guibeaux & Jérôme Jullien**

📄 **Published in *Nature Communications*, April 2025**
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41467--025--58615--7-blue)](https://doi.org/10.1038/s41467-025-58615-7)

---

## About this repository

This repository contains all bioinformatics scripts used for the data analysis published in François-Campion et al., *Nature Communications*, 2025.

Raw data are deposited on ENA under accession number **PRJEB56442**.

### Scientific context

This study investigates the role of the paternal epigenetic mark H2AK119ub1 (a Polycomb Repressive Complex 1 histone modification) in sperm and its transmission to the embryo during fertilisation. We demonstrate that sperm-derived H2AK119ub1 is required for proper early embryonic development in *Xenopus laevis*, establishing a link between paternal epigenetic information and developmental programming.

**Key analyses:**
- ChIP-seq profiling of H2AK119ub1, H3K4me3, H3K27me3 in sperm and embryos at multiple developmental stages
- Drosophila spike-in normalisation strategy with hybrid *Xenopus laevis* / *Drosophila melanogaster* genome alignment
- Bulk RNA-seq transcriptomic analysis across four species (*Xenopus laevis*, mouse, rabbit, human)
- Genomic characterisation, TF motif analysis, permutation tests

---

## Repository structure

| File | Description |
|------|-------------|
| `ChIP_Sequencing.sh` | Full ChIP-seq processing pipeline (alignment, filtering, peak calling) |
| `H2A_enrichment.sh` | H2AK119ub1 particle enrichment detection |
| `ProcessingChIPdata.R` | R scripts: permutations, TF motif heatmaps, genomic characterisation |
| `Chip.yml` | Conda environment for ChIP-seq analysis |
| `rnaseq_align.yml` | Conda environment for RNA-seq alignment |
| `MaskMerging.py` | Python script for N-mask hybrid genome merging |

**External scripts (GitLab Nantes Université):**
- RNA-seq alignment (Snakemake): https://gitlab.univ-nantes.fr/E114424Z/rnaseq_align
- RNA-seq differential expression (R): https://gitlab.univ-nantes.fr/E114424Z/BulkRNAseq

---

## Prerequisites

**Genome references:**
- *Xenopus laevis* v9.2 — [Xenbase](https://download.xenbase.org/xenbase/Genomics/JGI/Xenla9.2/)
- 601DNA spike-in sequence (100% ubiquitylated semi-synthetic nucleosome)
- CpG islands & repeat elements annotation — [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/database/)

**Tools:** see `Chip.yml` and `rnaseq_align.yml` for conda environments

---

## Citation

```bibtex
@article{francoiscampion2025,
  title={Sperm-derived H2AK119ub1 is required for embryonic development in Xenopus laevis},
  author={François-Campion, Valentin and Berger, Florian and Oikawa, Mami and Goumeidane, Maissa and Guibeaux, Romain and Jullien, Jérôme},
  journal={Nature Communications},
  year={2025},
  doi={10.1038/s41467-025-58615-7}
}
```

---

## Author

**Valentin FRANCOIS--CAMPION** — [GitHub](https://github.com/FCValentin) · [ResearchGate](https://www.researchgate.net/profile/Valentin-Francois-Campion)
