# clinstudtools/preprocessing/ingestion.py

import re
import os
import pandas as pd

from clinstudtools.utils import read_to_df

DEFAULT_PRINT_ORDER = ["Variable", "Site", "SampleID", "Investigator"]


""" Data Entry & Formatting Functions """
""" Pre-MethodComparator data reading """
def standardize_sample_ids(df, id_col="SampleID", no_dup=True):
    """
    Normalize sample IDs by removing site prefixes and leading zeros.

    Args:
        df (DataFrame): Input DataFrame with raw SampleIDs.
        id_col (str): Name of the sample ID column.

    Returns:
        DataFrame: DataFrame with standardized SampleIDs.
    """

    def clean_id(raw_id):
        if pd.isna(raw_id):
            return raw_id  # Preserve missing IDs

        # Extract trailing number (with optional leading zeros)
        match = re.search(r"(\d+(?:\.\d+)?)$", str(raw_id))
        if match:
            # numeric_part = match.group(1).lstrip("0") or "0"  # Preserve 0 if all zeros - delete this if the next line works for stripping 0s
            numeric_part = match.group(1)
            return str(int(float(numeric_part)))
        else:
            print(f"\033[91mCould not parse SampleID: {raw_id}\033[0m")
            return raw_id  # fallback: return raw

    df[id_col] = df[id_col].apply(clean_id)

    # check for sampleIDs duplicates
    if no_dup:
        duplicates = df["SampleID"][df["SampleID"].duplicated()]
        if not duplicates.empty:
            print(f"Duplicate SampleIDs after cleaning: {duplicates.unique()}")

    df["SampleID"] = df["SampleID"].str.zfill(5)  # 45 → 00045
    return df

def stnd_names(df, alias_map):
    rename_dict = {
        col: alias_map[col] for col in df.columns if col in alias_map
    }
    df = df.rename(columns=rename_dict)
    return df


def raw_to_df(file_name, site=None, method=None, sheet_name='Sheet1', dir=None):
    df = read_to_df(file_name, sheet_name, dir)

    # Standardize the SampleID
    possible_id_cols = ["SampleID", "Sample", "Sample ID", "ID", "Barcode", "barcode", "Case", "Anonymised no.", "Case ID"]
    id_col = next((col for col in possible_id_cols if col in df.columns), None)
    if not id_col:
        raise ValueError(f"No sample ID column found in {file_name}")
    df.rename(columns={id_col: "SampleID"}, inplace=True)

    if "SampleID" not in df.columns:
        raise ValueError(f"Missing 'SampleID' column in {file_name}")

    # use appropriate ID standardization method according to study type
    # if 'Mast cell' in df.columns and df['SampleID'][0][-2].isalpha():

    # If investigator/reviewer column exists, count number of investigators
    possible_inv_cols = ["Investigator", "Reviewer", "Investigator's Name", "Reviewer's Name", "Reviewer's full name"]
    inv_col = next((col for col in possible_inv_cols if col in df.columns), None)
    if inv_col:
        inv_str = ', '.join(df[inv_col].unique())
        num_inv = df[inv_col].nunique()
        print(f'{num_inv} investigators in dataframe: {inv_str}')
        df.rename(columns={inv_col: "Investigator"}, inplace=True)

    df = standardize_sample_ids(df, id_col="SampleID", no_dup=False)

    # Add metadata columns
    if isinstance(site, str):
        df["Site"] = site
    if isinstance(method, str):
        df["Method"] = method
    df["FileName"] = os.path.basename(file_name)

    # Strip leading/trailing spaces from column names
    df.columns = df.columns.str.strip()

    # to do: lowercase column names?
    return df


