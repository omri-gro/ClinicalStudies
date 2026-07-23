import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from regressions import deming_regression_pe

# ==========================================
# 1. Configuration
# ==========================================
STUDY_MODE = "PLT_SIZE"  # Change to "PLT_SIZE" or "RBC_SIZE"
SWEEP_DATA = rf"side_results/regression_sweep_{STUDY_MODE}.csv"
MAPPING_FILE = rf"./mapping/{STUDY_MODE}_tasks_mapping.csv"
OUTPUT_DIR = "./side_results"

# Define the grid search parameters
# Step size reduced to 0.5 per your parameters
CONFIGS = {
    'RBC Macrocyte': {'min': 7.0, 'max': 12.0, 'step': 0.5},
    'RBC Microcyte': {'min': 4.0, 'max': 8.0, 'step': 0.5},
    'PLT Large (includes Giant)': {'min': 4.0, 'max': 8.0, 'step': 0.5}
}


# ==========================================
# 2. Helper Functions
# ==========================================
def calculate_r(x, y):
    """Calculates Pearson R correlation."""
    mask = ~np.isnan(x) & ~np.isnan(y)
    if mask.sum() < 2:
        return 0.0
    r, _ = pearsonr(x[mask], y[mask])
    return r


# ==========================================
# 3. Execution Engine
# ==========================================
if __name__ == "__main__":
    print(f"Loading data from {SWEEP_DATA}...")
    try:
        sweep_df = pd.read_csv(SWEEP_DATA)
    except FileNotFoundError:
        print(f"Error: {SWEEP_DATA} not found.")
        exit()

    print(f"Loading mapping from {MAPPING_FILE}...")
    try:
        mapping_df = pd.read_csv(MAPPING_FILE)
    except FileNotFoundError:
        try:
            # Fallback to excel if mapping is .xlsx
            mapping_df = pd.read_excel(MAPPING_FILE.replace('.csv', '.xlsx'))
        except FileNotFoundError:
            print(f"Error: Could not find mapping file for {STUDY_MODE}.")
            exit()

    # 1. Extract 'Site' from mapping and merge to sweep dataframe
    # A single scan belongs to one site, so we drop duplicates
    site_map = mapping_df[['scan_id', 'Site']].drop_duplicates()
    df = pd.merge(sweep_df, site_map, on='scan_id', how='left')

    # Verify we mapped sites successfully
    if df['Site'].isna().any():
        missing = df[df['Site'].isna()]['scan_id'].unique()
        print(f"Warning: Missing 'Site' mapping for {len(missing)} scans.")
        df = df.dropna(subset=['Site'])

    target_classes = df['Target_Class'].unique()
    all_optimizations = []

    for target in target_classes:
        if target not in CONFIGS:
            print(f"Skipping {target}: No config defined.")
            continue

        print(f"\nOptimizing Site Thresholds for: {target}")
        class_df = df[df['Target_Class'] == target].copy()

        # 2. Average the reviewer ground truth AND AI ROI percentages per scan
        # Because Scan_Percent_at_Xum is identical for multiple reviewers on the same scan,
        # mean() safely collapses it without changing the value.
        scan_grouped = class_df.groupby(['scan_id', 'Site']).mean(numeric_only=True).reset_index()

        sites = scan_grouped['Site'].unique().tolist()
        print(f"  -> Found {len(sites)} sites: {sites}")
        print(f"  -> {len(scan_grouped)} total unique scans collapsed into mean points.")

        # 3. Build the threshold grid
        cfg = CONFIGS[target]
        # Use np.arange but round to 2 decimal to avoid floating point precision issues
        th_range = [round(t, 2) for t in np.arange(cfg['min'], cfg['max'] + cfg['step'], cfg['step'])]

        # Create all possible combinations of thresholds for the sites
        combinations = list(itertools.product(th_range, repeat=len(sites)))
        print(f"  -> Evaluating {len(combinations)} threshold combinations...")

        for combo in combinations:
            # Map the current threshold combination to the sites
            site_th_map = dict(zip(sites, combo))

            # Construct the synthetic arrays for both Full Scan and ROI
            y_scan_ai = []
            y_roi_ai = []
            x_truth = []

            for _, row in scan_grouped.iterrows():
                site = row['Site']
                th = site_th_map[site]

                # Retrieve the exact column name for this threshold
                scan_col = f"Scan_Percent_at_{th:.1f}um"
                roi_col = f"ROI_Percent_at_{th:.1f}um"

                if scan_col in row and pd.notna(row[scan_col]) and pd.notna(row['percent_in_roi_rev']):
                    x_truth.append(row['percent_in_roi_rev'])  # X Axis: Mean of reviewers' ROI %
                    y_scan_ai.append(row[scan_col])  # Y Axis 1: Scan AI % at target threshold
                    y_roi_ai.append(row[roi_col])  # Y Axis 2: Mean of ROI AI % at target threshold

            x_truth = np.array(x_truth)
            y_scan_ai = np.array(y_scan_ai)
            y_roi_ai = np.array(y_roi_ai)

            if len(x_truth) < 2:
                continue

            # 4. Run Deming & Correlation for Full Scan vs Reviewer ROI Mean
            scan_slope, scan_int = deming_regression_pe(x_truth, y_scan_ai, lambda_=1)
            scan_r = calculate_r(x_truth, y_scan_ai)

            # Distance from ideal (Slope=1, R=1)
            scan_dist = np.sqrt((1 - scan_slope) ** 2 + (1 - scan_r) ** 2)

            # 5. Run Deming & Correlation for AI ROI Mean vs Reviewer ROI Mean
            roi_slope, roi_int = deming_regression_pe(x_truth, y_roi_ai, lambda_=1)
            roi_r = calculate_r(x_truth, y_roi_ai)

            roi_dist = np.sqrt((1 - roi_slope) ** 2 + (1 - roi_r) ** 2)

            # 6. Log the results
            res_dict = {
                'Target_Class': target,
            }
            # Dynamically add the site thresholds to the output
            for site, th in site_th_map.items():
                res_dict[f'Threshold_{site}'] = th

            res_dict.update({
                # Full Scan Metrics (Main interest)
                'Scan_Slope': round(scan_slope, 4),
                'Scan_Pearson_R': round(scan_r, 4),
                'Scan_Ideal_Distance': round(scan_dist, 4),
                # ROI Metrics (Secondary interest)
                'ROI_Slope': round(roi_slope, 4),
                'ROI_Pearson_R': round(roi_r, 4),
                'ROI_Ideal_Distance': round(roi_dist, 4)
            })

            all_optimizations.append(res_dict)

    # ---------------------------------------------------------
    # Final Output Formatting
    # ---------------------------------------------------------
    if all_optimizations:
        opt_df = pd.DataFrame(all_optimizations)

        # Sort the output by the best Scan Ideal Distance so the best options are at the top
        opt_df = opt_df.sort_values(by=['Target_Class', 'Scan_Ideal_Distance'], ascending=[True, True])

        os.makedirs(OUTPUT_DIR, exist_ok=True)
        out_path = os.path.join(OUTPUT_DIR, f"site_threshold_optimization_{STUDY_MODE}.csv")
        opt_df.to_csv(out_path, index=False)
        print(f"\nOptimization complete! All regressions saved to {out_path}")
