import os
import pandas as pd
import numpy as np
import ast

# ==========================================
# Configuration
# ==========================================
STUDY_MODE = "PLT_SIZE"  # "PLT_SIZE" or "RBC_SIZE"

# File paths
CELLS_DATA = rf"side_results/full_scan_cells_{STUDY_MODE}.csv"
ROI_CELLS_DATA = rf"side_results/behavioral_cell_data_{STUDY_MODE}.csv"
RESULTS_DATA = rf"results/{STUDY_MODE}_results.csv"

# Exclusions mapping (Matches blobs_analysis.py logic for accurate full-scan tracking)
MACROCYTE_EXCLUSIONS = {
    "Sickle", "Bite", "Blister", "Spherocyte", "Acanthocyte", "Helmet",
    "Teardrop", "Schistocytes", "Poikilocyte", "Elliptocyte"
}
MICROCYTE_EXCLUSIONS = MACROCYTE_EXCLUSIONS.union({"Stomatocyte"})


def safe_eval_list(val):
    try:
        return ast.literal_eval(val) if pd.notna(val) else []
    except:
        return []


# ==========================================
# Sweep Configurations
# ==========================================
if STUDY_MODE == "PLT_SIZE":
    # Use only Max Size for Platelets
    size_col = 'Cell_Size_Max_um'
    sweep_configs = [
        {
            'name': 'PLT Large (includes Giant)',
            'results_morph': 'PLT large&giant',  # Matches the Morphology column in results_df
            'targets': ["Marked as PLT large", "Marked as PLT giant"],
            'direction': 'greater',
            'min': 4.0, 'max': 8.0,
            'exclusions': set()
        }
    ]
elif STUDY_MODE == "RBC_SIZE":
    # Use only Mean Size for RBCs
    size_col = 'Cell_Size_Mean_um'
    sweep_configs = [
        {
            'name': 'RBC Macrocyte',
            'results_morph': 'RBC macrocyte',
            'targets': ["Marked as RBC macrocyte"],
            'direction': 'greater',
            'min': 7.0, 'max': 12.0,
            'exclusions': MACROCYTE_EXCLUSIONS
        },
        {
            'name': 'RBC Microcyte',
            'results_morph': 'RBC microcyte',
            'targets': ["Marked as RBC microcyte"],
            'direction': 'less',
            'min': 4.0, 'max': 8.0,
            'exclusions': MICROCYTE_EXCLUSIONS
        }
    ]

# ==========================================
# Execution
# ==========================================
if __name__ == "__main__":
    print(f"Loading datasets for {STUDY_MODE}...")
    scan_cells_df = pd.read_csv(CELLS_DATA)
    roi_cells_df = pd.read_csv(ROI_CELLS_DATA)
    results_df = pd.read_csv(RESULTS_DATA)

    # Convert string representation of lists back to Python lists if the column exists
    if 'base_morph_list' in scan_cells_df.columns:
        scan_cells_df['base_morph_list'] = scan_cells_df['base_morph_list'].apply(safe_eval_list)

    step = 0.5
    all_configs_results = []

    print("Sweeping thresholds and calculating percentages...")
    for config in sweep_configs:
        print(f"  -> Processing {config['name']} ({config['direction']} than threshold)...")
        thresholds = np.arange(config['min'], config['max'] + step, step)

        # ---------------------------------------------------------
        # 1. Compute Full-Scan Percentages
        # ---------------------------------------------------------
        # Apply shape exclusions to the scan cells if applicable
        if 'base_morph_list' in scan_cells_df.columns and config['exclusions']:
            valid_scan_cells = scan_cells_df[
                ~scan_cells_df['base_morph_list'].apply(lambda morphs: any(m in config['exclusions'] for m in morphs))
            ]
        else:
            valid_scan_cells = scan_cells_df

        scan_summaries = []
        for scan_id, group in valid_scan_cells.groupby('scan_id'):
            # The denominator is always the total UNFILTERED base class cells in the scan
            total_scan_cells = len(scan_cells_df[scan_cells_df['scan_id'] == scan_id])
            if total_scan_cells == 0: continue

            sizes = group[size_col].values
            row = {'scan_id': scan_id}

            for th in thresholds:
                th_round = round(th, 2)
                if config['direction'] == 'greater':
                    count = (sizes > th).sum()
                else:
                    count = (sizes < th).sum()
                row[f'Scan_Percent_at_{th_round}um'] = round((count / total_scan_cells) * 100, 3)
            scan_summaries.append(row)

        scan_df = pd.DataFrame(scan_summaries)

        # ---------------------------------------------------------
        # 2. Compute ROI-Level AI Percentages
        # ---------------------------------------------------------
        roi_summaries = []
        # Group by both scan and reviewer to handle multiple reviewers correctly
        for (scan_id, reviewer), group in roi_cells_df.groupby(['scan_id', 'Reviewer_Name']):
            total_roi_cells = len(group)
            if total_roi_cells == 0: continue

            sizes = group[size_col].values
            row = {'scan_id': scan_id, 'Reviewer_Name': reviewer}

            for th in thresholds:
                th_round = round(th, 2)
                if config['direction'] == 'greater':
                    count = (sizes > th).sum()
                else:
                    count = (sizes < th).sum()
                row[f'ROI_Percent_at_{th_round}um'] = round((count / total_roi_cells) * 100, 3)
            roi_summaries.append(row)

        roi_df = pd.DataFrame(roi_summaries)

        # ---------------------------------------------------------
        # 3. Grab Reviewer Ground Truth & Merge
        # ---------------------------------------------------------
        # Filter the main results CSV for the specific morphology ground truth
        truth_subset = results_df[results_df['Morphology'] == config['results_morph']]
        if not truth_subset.empty:
            truth_subset = truth_subset[['scan_id', 'Reviewer Name', 'percent_in_roi_rev']].copy()
            # Standardize column name to match our ROI sweeps dataframe
            truth_subset.rename(columns={'Reviewer Name': 'Reviewer_Name'}, inplace=True)
            truth_subset.drop_duplicates(subset=['scan_id', 'Reviewer_Name'], inplace=True)

            # Merge Reviewer Truth with ROI sweeps
            merged_roi = pd.merge(truth_subset, roi_df, on=['scan_id', 'Reviewer_Name'], how='inner')
        else:
            print(f"    [Warning] Could not find '{config['results_morph']}' in {RESULTS_DATA}. Skipping merge.")
            merged_roi = roi_df

        # Merge everything together with the Scan sweeps
        if not scan_df.empty and not merged_roi.empty:
            final_config_df = pd.merge(merged_roi, scan_df, on='scan_id', how='left')
            final_config_df.insert(0, 'Target_Class', config['name'])
            all_configs_results.append(final_config_df)

    # ---------------------------------------------------------
    # Final Output
    # ---------------------------------------------------------
    if all_configs_results:
        final_regression_df = pd.concat(all_configs_results, ignore_index=True)

        # Ensure the output directory exists
        os.makedirs("side_results", exist_ok=True)
        output_name = rf"side_results/regression_sweep_{STUDY_MODE}.csv"

        final_regression_df.to_csv(output_name, index=False)
        print(f"\nDone! Regression data saved to {output_name}.")
    else:
        print("No data processed. Check if files exist and paths are correct.")

