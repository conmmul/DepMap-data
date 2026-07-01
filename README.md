## KRAS Allele-Specific Gene Dependencies in Colorectal Cancer

![Python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)
![Jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?style=for-the-badge&logo=Jupyter&logoColor=white)
![pandas](https://img.shields.io/badge/Pandas-2C2D72?style=for-the-badge&logo=pandas&logoColor=white)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

### Motivation
When I first started working on this project during my Summer 2024 research in the Westover lab (UTSW), the idea was to investigate KRAS co-mutations and how they affect lethality using the same general workflow as shown here. I worked for a while figuring out pandas and scipy to produce some volcano plots comparing co-mutant lethality. I stumbled upon the main findings of this project mostly on accident. I was looking up some of the genes that were consistently significant by ANOVA and noticed that many of them were involved in ribosome biogenesis, which has been an interesting field of study in cancer biology recently. This brought me to the idea that this could suggest that G12D cell lines may depend on ribosomal biogenesis regulation more than other mutants. This is consistent with findings that G12D preferentially signals through PI3K/AKT/mTOR, which regulates ribosome biogenesis. I learned so much working on this project from hours on stack overflow and different Python library documentations. I'm hoping these findings can be experimentally validated some day, but the work was rewarding regardless.

### Background

KRAS is the most frequently mutated oncogene in human cancer. In colorectal adenocarcinoma (COAD), approximately 40% of tumors harbor a KRAS mutation, with G12D being the most prevalent allele (~12% of all CRC). Different KRAS alleles have distinct biochemical properties — G12D preferentially activates PI3K signaling, G12V shows the strongest intrinsic GTPase impairment, and A146T drives faster nucleotide exchange with a weaker transformation potential — but whether these differences translate into allele-specific genetic vulnerabilities is not well established.

This project uses genome-wide CRISPR knockout data from the [DepMap](https://depmap.org/portal/) 24Q2 release to identify genes and pathways that are selectively essential in KRAS-mutant colorectal cancer cell lines. By comparing CRISPR gene effect scores (Chronos) between KRAS mutant and wild-type COAD lines, we nominate potential therapeutic vulnerabilities that may differ across alleles.

### Key Findings

- **Ribosome biogenesis is the top enriched pathway** in KRAS G12D COAD lines (KEGG_RIBOSOME: NES = −1.86, FDR = 0), followed by proteasome, spliceosome, RNA polymerase, and DNA replication pathways — all with FDR < 0.001
- **RPS7** is the top individual gene hit (p = 6.98 × 10⁻⁵), consistent with oncogenic KRAS elevating demand on protein synthesis machinery and with RPS7's known role in p53 activation under ribosomal stress
- One-way ANOVA across G12D, G12V, A146T, and WT groups identifies a shared set of core cellular machinery dependencies, with the strongest and most consistent signal in G12D lines
- The pattern is consistent with the oncogenic stress hypothesis: hyperproliferation driven by mutant KRAS creates elevated dependency on ribosomal and proteasomal function relative to non-transformed KRAS WT lines

### Repository Structure

```
kras-depmap/
├── notebooks/
│   ├── DepMap_GSEA.ipynb     # Differential dependency analysis + pre-ranked GSEA (G12D vs WT)
│   ├── anova.ipynb           # One-way ANOVA across G12D, G12V, A146T, and WT groups
│   └── boxplots.ipynb        # Publication-quality box plots for top hits
├── outputs/                  # Generated figures (not tracked by git)
├── environment/
│   └── requirements.txt      # Python dependencies with versions
└── README.md
```

### Data

All data are from the **DepMap Public 24Q2 release** and must be downloaded locally from the [DepMap portal](https://depmap.org/portal/). The specific files used are:

| File | Description |
|------|-------------|
| `CRISPRGeneEffect.csv` | Chronos CRISPR gene effect scores across all cell lines |
| `OmicsSomaticMutations.csv` | Somatic mutation calls per cell line |
| `Model.csv` | Cell line metadata including cancer type (OncotreeCode) |
| `CRISPRGeneDependency.csv` | Probabilistic dependency scores |

The data used in this analysis have also been archived at Figshare:
**DOI:** [10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)

Gene sets (KEGG, PID, BioCarta) are from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) v2024.1.

### Reproducing the Analysis

**1. Clone the repository**
```bash
git clone https://github.com/conmmul/kras-depmap.git
cd kras-depmap
```

**2. Install dependencies**
```bash
pip install -r environment/requirements.txt
```

**3. Download data**

Download the four DepMap files listed above from the [DepMap portal](https://depmap.org/portal/) and place them in a local directory. Download the KEGG, PID, and BioCarta gene sets from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/).

**4. Configure paths**

Each notebook has a config cell near the top. Update `DATA_DIR` and (for the GSEA notebook) `GENESETS_DIR` to point to your local copies of the data.

**5. Run notebooks**

Run the notebooks in this order:
1. `DepMap_GSEA.ipynb` — differential dependency and pathway enrichment
2. `anova.ipynb` — multi-allele ANOVA
3. `boxplots.ipynb` — visualization of top hits

### Dependencies

See `environment/requirements.txt`. Core dependencies:

- Python 3.10+
- pandas, numpy, scipy, matplotlib
- gseapy
- statsmodels

### Methods Summary

COAD cell lines were stratified by KRAS genotype using `OmicsSomaticMutations.csv` and `Model.csv`. The wild-type group excludes all KRAS-mutant lines regardless of allele. For the two-group comparison (G12D vs WT), a Welch's t-test was applied across all genes followed by Benjamini-Hochberg FDR correction. Genes were ranked by mean Chronos score in G12D lines and used as input for pre-ranked GSEA against KEGG gene sets (gseapy, 1000 permutations). For the multi-allele comparison, a one-way ANOVA was applied across G12D, G12V, A146T, and WT groups with BH correction.
