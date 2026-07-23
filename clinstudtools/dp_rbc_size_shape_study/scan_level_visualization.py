import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from regressions import regression_comp

# ==========================================
# 1. Configuration
# ==========================================
STUDY_MODE = "RBC_SIZE"  # Change to "PLT_SIZE" or "RBC_SIZE"
DATA_FILE = rf"side_results/regression_sweep_{STUDY_MODE}.csv"
OUTPUT_DIR = "./visualizations/regression"

# Set visualization style
sns.set_theme(style="whitegrid", palette="muted")


# ==========================================
# 2. Helper Functions
# ==========================================
def calculate_r_squared(x, y):
    """Calculates R^2 using Pearson correlation, dropping NaNs safely."""
    mask = ~np.isnan(x) & ~np.isnan(y)
    if mask.sum() < 2:
        return 0.0
    r, _ = pearsonr(x[mask], y[mask])
    return r ** 2


def get_threshold_columns(df, prefix):
    """Extracts a sorted list of threshold columns and their float values."""
    cols = [c for c in df.columns if c.startswith(prefix)]
    th_mapping = {c: float(re.search(r"at_(\d+\.\d+)um", c).group(1)) for c in cols}
    sorted_cols = sorted(cols, key=lambda c: th_mapping[c])
    return sorted_cols, [th_mapping[c] for c in sorted_cols]


# ==========================================
# 3. Plotting Functions
# ==========================================
def plot_correlation_curve(thresholds, scan_r2, roi_r2, target_class, best_scan_th, best_roi_th):
    """Plots how R^2 changes as the threshold shifts."""
    plt.figure(figsize=(10, 6))

    plt.plot(thresholds, scan_r2, label=f"Full Scan AI (Peak: {best_scan_th}μm)", linewidth=2.5, marker='o',
             markersize=4)
    plt.plot(thresholds, roi_r2, label=f"ROI AI (Peak: {best_roi_th}μm)", linewidth=2.5, marker='s', markersize=4)

    plt.axvline(x=best_scan_th, color='blue', linestyle='--', alpha=0.5)
    plt.axvline(x=best_roi_th, color='orange', linestyle='--', alpha=0.5)

    plt.title(f"Threshold Optimization ($R^2$ Correlation) - {target_class}", fontsize=14, pad=15)
    plt.xlabel("Cell Size Threshold (μm)", fontsize=12)
    plt.ylabel("$R^2$ Score (Agreement with Reviewers)", fontsize=12)
    plt.legend(title="AI Context Window")

    plt.tight_layout()
    safe_name = target_class.replace(' ', '_').replace('&', 'and')
    plt.savefig(os.path.join(OUTPUT_DIR, f"{safe_name}_R2_Curve.png"), dpi=300)
    plt.close()


def plot_optimal_scatter(df, x_col, y_col, target_class, title_prefix, context_type):
    """Plots the Deming regression scatter at the optimal threshold."""
    plt.figure(figsize=(8, 8))

    # Drop NaNs before passing to Deming regression to avoid matrix math errors
    valid_mask = ~df[x_col].isna() & ~df[y_col].isna()
    x_clean = df.loc[valid_mask, x_col].values
    y_clean = df.loc[valid_mask, y_col].values

    if len(x_clean) < 2:
        print(f"  -> Not enough data to plot {title_prefix}")
        plt.close()
        return

    # Call your custom Deming Regression script
    reg_dict = regression_comp(x_clean, y_clean, reg_method="deming", res_str=True)

    slope = reg_dict["slope"]
    intercept = reg_dict["intercept"]
    r2 = reg_dict["correlation_coefficient"] ** 2

    # Format the strings from your script by removing newlines for a cleaner legend
    slope_str = reg_dict["slope_str"].replace('\n', ' ')
    int_str = reg_dict["intercept_str"].replace('\n', ' ')

    # Plot the raw scatter points
    sns.scatterplot(
        x=x_clean,
        y=y_clean,
        alpha=0.6,
        edgecolor='w',
        s=50,
        color='steelblue'
    )

    # Plot the Deming Regression line mathematically
    x_range = np.array([0, x_clean.max()])
    y_pred = slope * x_range + intercept

    plt.plot(
        x_range,
        y_pred,
        color='red',
        linewidth=2,
        label=f"Deming Fit\nSlope: {slope_str}\nIntercept: {int_str}"
    )

    plt.title(f"{title_prefix}\n{target_class} ($R^2$ = {r2:.3f})", fontsize=14, pad=15)
    plt.xlabel("Reviewer Marked Percentage in ROI (%)", fontsize=12)
    plt.ylabel(f"AI Calculated Percentage in {context_type} (%)", fontsize=12)

    # Ensure axes start near 0 and maintain a square aspect ratio
    max_val = max(x_clean.max(), y_clean.max()) * 1.1
    plt.xlim(-1, max_val)
    plt.ylim(-1, max_val)

    # Draw a faint y=x reference line
    plt.plot([-1, max_val], [-1, max_val], color='gray', linestyle=':', alpha=0.5, label='Perfect 1:1 Agreement')
    plt.legend(loc='best')

    plt.tight_layout()
    safe_name = target_class.replace(' ', '_').replace('&', 'and')
    plt.savefig(os.path.join(OUTPUT_DIR, f"{safe_name}_{context_type}_Optimal_Scatter.png"), dpi=300)
    plt.close()


# ==========================================
# 4. Execution Logic
# ==========================================
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Loading regression data from {DATA_FILE}...")
    try:
        df = pd.read_csv(DATA_FILE)
    except FileNotFoundError:
        print(f"Error: {DATA_FILE} not found. Run scan_level_regression.py first.")
        exit()

    target_classes = df['Target_Class'].unique()

    for target in target_classes:
        print(f"\nAnalyzing: {target}")
        class_df = df[df['Target_Class'] == target]

        scan_cols, thresholds = get_threshold_columns(class_df, 'Scan_Percent_at_')
        roi_cols, _ = get_threshold_columns(class_df, 'ROI_Percent_at_')

        if not scan_cols or not roi_cols:
            print(f"  -> Skipping. Missing threshold columns for {target}.")
            continue

        reviewer_vals = class_df['percent_in_roi_rev'].values

        # 1. Sweep to find R^2 for every threshold to identify the peak
        scan_r2_scores = []
        roi_r2_scores = []

        for s_col, r_col in zip(scan_cols, roi_cols):
            scan_r2_scores.append(calculate_r_squared(reviewer_vals, class_df[s_col].values))
            roi_r2_scores.append(calculate_r_squared(reviewer_vals, class_df[r_col].values))

        # 2. Identify the peak thresholds
        best_scan_idx = np.argmax(scan_r2_scores)
        best_roi_idx = np.argmax(roi_r2_scores)

        best_scan_th = thresholds[best_scan_idx]
        best_roi_th = thresholds[best_roi_idx]

        best_scan_col = scan_cols[best_scan_idx]
        best_roi_col = roi_cols[best_roi_idx]

        print(f"  -> Peak Full-Scan AI threshold: {best_scan_th}μm (R^2 = {scan_r2_scores[best_scan_idx]:.3f})")
        print(f"  -> Peak ROI AI threshold: {best_roi_th}μm (R^2 = {roi_r2_scores[best_roi_idx]:.3f})")

        # 3. Generate the visualizations
        plot_correlation_curve(
            thresholds, scan_r2_scores, roi_r2_scores,
            target_class=target, best_scan_th=best_scan_th, best_roi_th=best_roi_th
        )

        plot_optimal_scatter(
            class_df, x_col='percent_in_roi_rev', y_col=best_scan_col,
            target_class=target, title_prefix=f"Full-Scan AI Regression at {best_scan_th}μm",
            context_type="Full_Scan"
        )

        plot_optimal_scatter(
            class_df, x_col='percent_in_roi_rev', y_col=best_roi_col,
            target_class=target, title_prefix=f"ROI AI Regression at {best_roi_th}μm",
            context_type="ROI"
        )

    print(f"\nComplete! Visualizations saved to '{OUTPUT_DIR}'.")