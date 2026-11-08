import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


# ==========================================
# 1. PROJECT PATHS
# ==========================================

project_folder = Path(__file__).resolve().parent.parent

file_path = (
    project_folder
    / "data"
    / "GSE352947_periodontitis_smallRNAseq_processed.xlsx"
)

results_folder = project_folder / "results"
results_folder.mkdir(parents=True, exist_ok=True)

print("===================================")
print("GENE EXPRESSION DATA ANALYSIS")
print("===================================")

if not file_path.exists():
    raise FileNotFoundError(f"Dataset file not found: {file_path}")

print("Dataset found successfully!")


# ==========================================
# 2. EXCEL WORKBOOK INSPECTION
# ==========================================

excel_file = pd.ExcelFile(file_path)

print("\nAvailable Sheets:")
print(excel_file.sheet_names)

for sheet in ["Info", "Group", "Data"]:
    print("\n" + "=" * 50)
    print("SHEET:", sheet)
    print("=" * 50)

    sheet_df = pd.read_excel(
        file_path,
        sheet_name=sheet,
        header=None
    )

    print("Shape:", sheet_df.shape)
    print(sheet_df.head(10).to_string())


# ==========================================
# 3. READ DATA WITH CORRECT HEADERS
# ==========================================

df_clean = pd.read_excel(
    file_path,
    sheet_name="Data",
    header=1
)

print("\n========== CORRECT DATA HEADERS ==========")
print(df_clean.columns.tolist())

print("\n========== FIRST 3 DATA ROWS ==========")
print(df_clean.head(3).to_string())

print("\n========== DATASET SIZE ==========")
print("Rows:", df_clean.shape[0])
print("Columns:", df_clean.shape[1])


# ==========================================
# 4. CREATE CLEAN DATASET
# ==========================================

sample_cols = [
    "C1", "C2", "C3", "C4", "C5", "C6",
    "P1", "P2", "P3", "P4", "P5", "P6"
]

control_cols = ["C1", "C2", "C3", "C4", "C5", "C6"]
patient_cols = ["P1", "P2", "P3", "P4", "P5", "P6"]

if df_clean.shape[1] < 18:
    raise ValueError("Dataset does not contain expected sample columns.")

analysis_df = df_clean.iloc[
    :, [0, 1] + list(range(6, 18))
].copy()

analysis_df.columns = ["ID", "miRNA"] + sample_cols

for col in sample_cols:
    analysis_df[col] = pd.to_numeric(
        analysis_df[col],
        errors="coerce"
    )

analysis_df["miRNA"] = analysis_df["miRNA"].astype(str)

print("\n========== CLEAN DATASET ==========")
print("Shape:", analysis_df.shape)

print("\nFirst 5 Rows:")
print(analysis_df.head().to_string())

print("\nMissing Values:")
print(analysis_df.isnull().sum())

print("\nData Types:")
print(analysis_df.dtypes)


# ==========================================
# 5. GROUP COMPARISON
# ==========================================

analysis_df["Control_Average"] = analysis_df[control_cols].mean(axis=1)
analysis_df["Patient_Average"] = analysis_df[patient_cols].mean(axis=1)

analysis_df["Expression_Difference"] = (
    analysis_df["Patient_Average"]
    - analysis_df["Control_Average"]
)

print("\n========== GROUP COMPARISON ==========")

print(
    analysis_df[
        [
            "miRNA",
            "Control_Average",
            "Patient_Average",
            "Expression_Difference"
        ]
    ].head(10).to_string(index=False)
)

print("\n========== SUMMARY ==========")

print("Average control expression:", analysis_df["Control_Average"].mean())
print("Average patient expression:", analysis_df["Patient_Average"].mean())
print(
    "Higher in patients:",
    (analysis_df["Expression_Difference"] > 0).sum()
)
print(
    "Lower in patients:",
    (analysis_df["Expression_Difference"] < 0).sum()
)
print(
    "Equal expression:",
    (analysis_df["Expression_Difference"] == 0).sum()
)


# ==========================================
# 6. TOP 10 HIGHER AND LOWER miRNAs
# ==========================================

top_higher = analysis_df.nlargest(10, "Expression_Difference")
top_lower = analysis_df.nsmallest(10, "Expression_Difference")

print("\n========== TOP 10 HIGHER IN PATIENTS ==========")

print(
    top_higher[
        [
            "miRNA",
            "Control_Average",
            "Patient_Average",
            "Expression_Difference"
        ]
    ].to_string(index=False)
)

print("\n========== TOP 10 LOWER IN PATIENTS ==========")

print(
    top_lower[
        [
            "miRNA",
            "Control_Average",
            "Patient_Average",
            "Expression_Difference"
        ]
    ].to_string(index=False)
)


# ==========================================
# 7. EXPRESSION DISTRIBUTION BOXPLOT
# ==========================================

plot_df = analysis_df.melt(
    id_vars=["ID", "miRNA"],
    value_vars=sample_cols,
    var_name="Sample",
    value_name="Expression"
)

plot_df["Group"] = plot_df["Sample"].apply(
    lambda x: "Control" if x.startswith("C") else "Patient"
)

plt.figure(figsize=(12, 6))

sns.boxplot(
    data=plot_df,
    x="Sample",
    y="Expression",
    hue="Group"
)

plt.title("miRNA Expression Distribution Across Samples")
plt.xlabel("Samples")
plt.ylabel("Normalized Expression (log2)")
plt.tight_layout()

boxplot_path = results_folder / "expression_distribution.png"
plt.savefig(boxplot_path, dpi=300)
plt.show()
plt.close()

print("\nBoxplot saved:", boxplot_path)


# ==========================================
# 8. HEATMAP - TOP 20 VARIABLE miRNAs
# ==========================================

expression_data = analysis_df.set_index("miRNA")[sample_cols]

top_variable = expression_data.var(axis=1).nlargest(20)
heatmap_data = expression_data.loc[top_variable.index]

row_std = heatmap_data.std(axis=1).replace(0, np.nan)

heatmap_data = (
    heatmap_data.sub(heatmap_data.mean(axis=1), axis=0)
    .div(row_std, axis=0)
).fillna(0)

plt.figure(figsize=(12, 9))

sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",
    center=0,
    linewidths=0.3
)

plt.title("Top 20 Variable miRNAs Across Samples")
plt.xlabel("Samples")
plt.ylabel("miRNAs")
plt.tight_layout()

heatmap_path = results_folder / "top20_miRNA_heatmap.png"
plt.savefig(heatmap_path, dpi=300)
plt.show()
plt.close()

print("Heatmap saved:", heatmap_path)


# ==========================================
# 9. STATISTICAL ANALYSIS - WELCH T-TEST
# ==========================================

results = []

for _, row in analysis_df.iterrows():
    control_values = row[control_cols].astype(float)
    patient_values = row[patient_cols].astype(float)

    statistic, p_value = ttest_ind(
        patient_values,
        control_values,
        equal_var=False,
        nan_policy="omit"
    )

    results.append({
        "miRNA": row["miRNA"],
        "Control_Average": row["Control_Average"],
        "Patient_Average": row["Patient_Average"],
        "Expression_Difference": row["Expression_Difference"],
        "P_Value": p_value
    })

stats_df = pd.DataFrame(results)
stats_df["Significant"] = stats_df["P_Value"] < 0.05

print("\n========== STATISTICAL ANALYSIS ==========")
print(stats_df.head(10).to_string(index=False))

print(
    "\nSignificant miRNAs (unadjusted p < 0.05):",
    stats_df["Significant"].sum()
)


# ==========================================
# 10. INVESTIGATE INVALID P-VALUES
# ==========================================

invalid_mask = ~np.isfinite(
    pd.to_numeric(stats_df["P_Value"], errors="coerce")
)

invalid_miRNAs = stats_df.loc[invalid_mask, "miRNA"]

invalid_data = analysis_df[
    analysis_df["miRNA"].isin(invalid_miRNAs)
].copy()

invalid_data["Control_Variance"] = invalid_data[control_cols].var(axis=1)
invalid_data["Patient_Variance"] = invalid_data[patient_cols].var(axis=1)

print("\n========== INVALID P-VALUE ANALYSIS ==========")
print("Invalid miRNAs:", len(invalid_data))

if len(invalid_data) > 0:
    print(
        invalid_data[
            [
                "miRNA",
                "Control_Variance",
                "Patient_Variance"
            ]
        ].head(20).to_string(index=False)
    )

    print("\nInvalid p-value variance summary:")
    print(
        invalid_data[
            ["Control_Variance", "Patient_Variance"]
        ].describe()
    )
else:
    print("No invalid p-values found.")


# ==========================================
# 11. FDR CORRECTION - BENJAMINI-HOCHBERG
# ==========================================

p_values = pd.to_numeric(
    stats_df["P_Value"],
    errors="coerce"
).to_numpy(dtype=float)

valid_mask = np.isfinite(p_values)

adjusted = np.full(len(stats_df), np.nan)

if valid_mask.any():
    adjusted[valid_mask] = multipletests(
        p_values[valid_mask],
        method="fdr_bh"
    )[1]

stats_df["FDR_Adjusted_P_Value"] = adjusted

stats_df["FDR_Significant"] = (
    np.isfinite(adjusted) & (adjusted < 0.05)
)

print("\n========== FDR CORRECTION ==========")

print(
    stats_df[
        [
            "miRNA",
            "P_Value",
            "FDR_Adjusted_P_Value",
            "FDR_Significant"
        ]
    ].head(10).to_string(index=False)
)


# ==========================================
# 12. FDR DIAGNOSTIC CHECK
# ==========================================

print("\n========== FDR DIAGNOSTIC ==========")

print("Total p-values:", len(p_values))
print("Valid p-values:", valid_mask.sum())
print("Invalid p-values:", (~valid_mask).sum())

if valid_mask.any():
    print("Minimum valid p-value:", p_values[valid_mask].min())
    print("Maximum valid p-value:", p_values[valid_mask].max())

print("Adjusted values count:", np.isfinite(adjusted).sum())

print(
    "Significant miRNAs after FDR correction:",
    stats_df["FDR_Significant"].sum()
)


# ==========================================
# 13. SAVE STATISTICAL RESULTS
# ==========================================

stats_output = results_folder / "statistical_analysis.csv"
stats_df.to_csv(stats_output, index=False)

print("\nStatistical results saved:", stats_output)


# ==========================================
# 14. FINAL SUMMARY
# ==========================================

print("\n===================================")
print("ANALYSIS COMPLETED SUCCESSFULLY")
print("===================================")

print("Total miRNAs analyzed:", len(analysis_df))
print("Unadjusted significant miRNAs:", stats_df["Significant"].sum())
print("FDR significant miRNAs:", stats_df["FDR_Significant"].sum())

print("\nGenerated Files:")
print("- expression_distribution.png")
print("- top20_miRNA_heatmap.png")
print("- statistical_analysis.csv")

print("\nResults folder:", results_folder)