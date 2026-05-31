import pandas as pd
import os

# Math Engine
from processing import process_repeatability, process_reproducibility, export_results
# Printer Engine
from report_creation import generate_docx

# --- CONFIGURATION ---
FILE_REPEATABILITY = r"raw/bma_repeatability.csv"
FILE_REPRODUCIBILITY = r"raw/bma_reproducibility.csv"
OUTPUT_DIR = "results"

def main():
    print("Starting Precision Analysis...")
    rep_results = []
    repro_results = []

    # Ensure output directory exists for our CSVs
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---------------------------------------------------------
    # 1. PROCESS REPEATABILITY (Single-Site)
    # ---------------------------------------------------------
    if os.path.exists(FILE_REPEATABILITY):
        print(f"Loading {FILE_REPEATABILITY}...")
        df_rep = pd.read_csv(FILE_REPEATABILITY)

        # Identify parameters by excluding known metadata columns
        exclude_cols = ['Sample', 'Day', 'Run', 'Scan', 'totalWBC',
                        'Unclassified', 'megakaryocyte', 'stripped',
                        'Particles', 'Scan_UUID']
        rep_params = [c for c in df_rep.columns if c not in exclude_cols]

        print("Calculating Repeatability variance components...")
        rep_results = process_repeatability(df_rep, rep_params)

    if rep_results:
        export_results(rep_results, os.path.join(OUTPUT_DIR, "Repeatability_Results.csv"))

    else:
        print(f"Skipping Repeatability: {FILE_REPEATABILITY} not found.")

    # ---------------------------------------------------------
    # 2. PROCESS REPRODUCIBILITY (Multi-Site)
    # ---------------------------------------------------------
    if os.path.exists(FILE_REPRODUCIBILITY):
        print(f"\nLoading {FILE_REPRODUCIBILITY}...")
        df_repro = pd.read_csv(FILE_REPRODUCIBILITY)

        # Identify parameters (Note 'Machine' instead of 'Run')
        exclude_cols = ['Sample', 'Machine', 'Day', 'Scan', 'totalWBC',
                        'Unclassified', 'megakaryocyte', 'stripped',
                        'Particles', 'Scan_UUID']
        repro_params = [c for c in df_repro.columns if c not in exclude_cols]

        print("Calculating Reproducibility variance components...")
        repro_results = process_reproducibility(df_repro, repro_params)

        if repro_results:
            export_results(repro_results, os.path.join(OUTPUT_DIR, "Reproducibility_Results.csv"))

    else:
        print(f"Skipping Reproducibility: {FILE_REPRODUCIBILITY} not found.")

    # ---------------------------------------------------------
    # 3. GENERATE REPORT
    # ---------------------------------------------------------
    if rep_results or repro_results:
        print("\nGenerating Consolidated Word Document...")
        generate_docx(rep_results, repro_results)
        print("Done! Analysis complete.")
    else:
        print("\nNo data was processed. Exiting.")


if __name__ == "__main__":
    main()

