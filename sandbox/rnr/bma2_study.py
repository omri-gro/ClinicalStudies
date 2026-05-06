import pandas as pd
import numpy as np
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from processing import process_repeatability, process_reproducibility

warnings.simplefilter('ignore', ConvergenceWarning)

# --- CONFIGURATION ---
FILE_REPEATABILITY = r"raw/bma_repeatability.csv"
FILE_REPRODUCIBILITY = r"raw/bma_reproducibility.csv"

def main():
    # 1. Process Repeatability CSV (if exists)
    try:
        df_rep = pd.read_csv(FILE_REPEATABILITY)
        # Extract parameter columns (assuming they start after the 4 metadata columns)
        param_cols = [c for c in df_rep.columns if c not in ['Sample', 'Day', 'Run', 'Scan']]
        process_repeatability(df_rep, param_cols)
    except FileNotFoundError:
        print(f"File {FILE_REPEATABILITY} not found. Skipping single-site analysis.")

    # 2. Process Reproducibility CSV (if exists)
    try:
        df_repro = pd.read_csv(FILE_REPRODUCIBILITY)
        # Extract parameter columns
        param_cols = [c for c in df_repro.columns if c not in ['Sample', 'Machine', 'Day', 'Scan']]
        process_reproducibility(df_repro, param_cols)
    except FileNotFoundError:
        print(f"File {FILE_REPRODUCIBILITY} not found. Skipping multisite analysis.")


if __name__ == "__main__":
    main()

