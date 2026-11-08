# Gene Expression Data Analysis – Periodontitis

## Project Overview

This project analyzes small RNA sequencing data to compare miRNA expression levels between normal control samples and patients with periodontitis.

The goal is to identify miRNAs with higher or lower expression in patients and explore differences using statistical analysis and data visualization.

## Dataset Information

- **GEO Accession:** GSE352947
- **Organism:** Human
- **Data Type:** Small RNA-Seq
- **Normalization:** TMM + CPM normalization
- **Total miRNAs:** 2,588
- **Control Samples:** 6 (C1–C6)
- **Patient Samples:** 6 (P1–P6)

Dataset source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352947

## Tools and Technologies

- Python
- Pandas
- NumPy
- Matplotlib
- Seaborn
- SciPy
- Statsmodels
- Excel

## Project Structure

```text
CRG_Project_02_Gene_Expression_Analysis/
│
├── data/
│   └── GSE352947_periodontitis_smallRNAseq_processed.xlsx
│
├── notebooks/
│   └── gene_expression_analysis.py
│
├── results/
│   ├── expression_distribution.png
│   ├── top20_miRNA_heatmap.png
│   └── statistical_analysis.csv
│
└── README.md
```

## Analysis Workflow

1. Loaded and inspected the Excel dataset.
2. Cleaned the miRNA expression data.
3. Selected normalized log2 expression values for 12 samples.
4. Calculated average expression for control and patient groups.
5. Identified the top 10 higher and lower expressed miRNAs.
6. Created an expression distribution boxplot.
7. Generated a heatmap of the 20 most variable miRNAs.
8. Performed Welch's independent-samples t-tests.
9. Applied Benjamini-Hochberg False Discovery Rate (FDR) correction.
10. Exported statistical results to CSV.

## Key Findings

- 2,588 miRNAs were analyzed.
- 1,157 miRNAs showed higher average expression in patients.
- 971 miRNAs showed lower average expression in patients.
- 460 miRNAs had invalid p-values in the Welch t-test analysis.
- 176 miRNAs had unadjusted p-values below 0.05.
- No miRNAs remained statistically significant after FDR correction at the 0.05 threshold.

## Visualizations

### Expression Distribution

The boxplot compares normalized miRNA expression across all control and patient samples.

### miRNA Heatmap

The heatmap displays the 20 most variable miRNAs across the 12 samples, using row-wise standardized expression values.

## Statistical Notes and Limitations

- Each group contains only six samples, so statistical power is limited.
- Welch's t-test was used to compare the two independent groups.
- Multiple testing correction was performed using the Benjamini-Hochberg method.
- The 460 invalid p-values require careful interpretation because the corresponding expression data had zero variance within both groups.
- The analysis is exploratory and does not establish causation or clinical significance.
- No statistically significant miRNAs were identified after FDR correction.

## How to Run the Project

Install the required libraries:

```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels openpyxl
```

Run the analysis from the project root directory:

```bash
python notebooks/gene_expression_analysis.py
```

Generated graphs and statistical results will be saved in the `results` folder.

## Author

Gene Expression Data Analysis Project

Created as part of a bioinformatics and research portfolio.