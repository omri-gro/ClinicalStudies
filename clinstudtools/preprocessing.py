# preprocessing.py

import re
import os
import pandas as pd

from .transforms import robust_dup as _robust_dup
from .core import MetadataBundle
from .utils import read_to_df



DEFAULT_PRINT_ORDER = ["Variable", "Site", "SampleID", "Investigator"]


""" transforms wrappers """
def robust_dup(df, *args, **kwargs):
    """Project-aware version of robust_dup."""
    sort_cols = [c for c in ["Variable", "SampleID"] if c in df.columns]
    return _robust_dup(df, *args, sort_cols=sort_cols, **kwargs)


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

# to do - split and clean into pipelines
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

    # removal of duplicates only if no investigator column (need to improve handling of investigator existence in main pipeline)
    if inv_col:  # consider adding option for removing multiple rows with same SampleID and investigator
        df = standardize_sample_ids(df, id_col="SampleID", no_dup=False)
    else:
        df = standardize_sample_ids(df, id_col="SampleID")


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



# to do!! - clean and organize all of the following (split module)
# to do - see if metadata can be removed to convert to function into general use + wrapper
def fill_nans(df, metadata, target, fill_value=None, condition_mask=None):
    """
    Fill NaNs in specified columns or group of columns, optionally conditioned on a mask.

    Args:
        df (DataFrame): DataFrame to operate on.
        metadata (MetadataBundle): Context object with metadata.
        target (str): Column name or group name.
        fill_value (any): If None, use metadata to determine value.
        condition_mask (pd.Series, optional): Boolean mask of rows where filling applies.
    Returns:
        DataFrame: Modified DataFrame.
    """
    # Resolve target to list of variables
    if target in df.columns:
        columns = [target]
    else:
        columns = metadata.variable_groups.get(target, [])
        columns = [col for col in columns if col in df.columns]

    if not columns:
        print(f"No matching columns found for target '{target}'")
        return df

    for col in columns:
        dtype = metadata.variables[col].get("data_type", "numeric")
        val = fill_value
        if val is None:
            if dtype == "numeric":
                val = 0
            elif dtype == "binary":
                val = False  # or "Not present"
            else:
                val = None  # fallback, don't fill

        if val is not None:
            if condition_mask is not None:
                rows_to_fill = df[col].isna() & condition_mask
            else:
                rows_to_fill = df[col].isna()

            num_filled = rows_to_fill.sum()
            df.loc[rows_to_fill, col] = val

            if num_filled > 0:
                print(f"Filled {num_filled} NaNs in '{col}' with {val} (conditioned: {target is not None})")
    return df


def handle_fill_nans_by_rule(df, metadata, rules):
    """
    Apply multiple NaN filling rules to a DataFrame.

    Iterates through a list of rules, each specifying a target column or group,
    an optional fill value, and an optional condition. Calls `fill_nans` for each
    rule to fill NaNs accordingly.

    Args:
        df (pd.DataFrame): The DataFrame to modify.
        metadata (MetadataBundle): Object with metadata and variable groups.
        rules (list): List of dicts defining fill rules. Each rule may include:
            - target (str): Column name or group name.
            - fill_value (any, optional): Value to fill NaNs with.
            - condition_column (str, optional): Column to evaluate condition on.
            - condition_values (list, optional): Values to match for conditional fill.

    Returns:
        pd.DataFrame: Modified DataFrame with NaNs filled as specified.

    Example:
        rules = [
            {"target": "WBC diff"},
            {"target": "RBC morphology", "condition_column": "Review Status",
             "condition_values": ["Reviewed", "Unremarkable"]}
        ]
        handle_fill_nans_by_rule(df, context, rules)
    """
    for rule in rules:
        target = rule["target"]
        fill_value = rule.get("fill_value")
        condition_mask = None

        # if a condition was defined
        if "condition_column" in rule and "condition_values" in rule:
            condition_mask = df[rule["condition_column"]].isin(rule["condition_values"])
            print(f"Applying conditional fill for '{target}' where {rule['condition_column']} in {rule['condition_values']}")

        df = fill_nans(df, metadata, target, fill_value, condition_mask)
    return df


def apply_src_fixes(df, metadata, src):
    """
    Apply site or method specific corrections using specific rules.
    Args:
        df (DataFrame): DataFrame to operate on.
        metadata (MetadataBundle): Object containing source fixes.
        src (str): Name of site and/or method (which appears in the metadata's source fixes).

    - Note "raw_given_as_grade" is handled during add_grade_column, and is ignored here.
    Returns:
        DataFrame: Modified DataFrame.
    """
    # get the source-specific fix directions
    fixes = metadata.src_fixes.get(src, {})
    if not fixes:
        print(f"No source-specific fixes for {src}")
        return df

    for rules_group, rules in fixes.items():
        # to do: add function for converting values (words to numeric grades or to True/False)
        # to do: add function that raises messages for unexpected values (like other words in "RBC & PLT Morphology")
        # to do: raise messages for unexpected data types (words in a numeric variable)
        if rules_group == "fill_nans_rules":
            df = handle_fill_nans_by_rule(df, metadata, rules)
        elif rules_group == "rename_rules":
            # df = handle_rename_rules(df, rules)
            pass
        elif rules_group == "raw_given_as_grade":
            pass
        else:
            print(f"\033[91Unknown fix rules group: {rules_group}\033[0m")

    return df



def one_to_one_hundered(df, metadata, diff_cells='percent', diff_like_cells='percent-like'):
    """
    If differential values provided as decimal (all <=1), convert to a percent-like (0-100).

    Args:
        df(pd.DataFrame): The DataFrame to check.
        metadata (MetadataBundle): Object with variable groups metadata.
        diff_cells (str, optional): Name of variables group (appearing in metadata) to be converted.
        diff_like_cells (str, optional): Name of variables group allowed to be >1 but will still be converted (e.g., can be 110%).
    """
    diff_vars = metadata.variable_groups.get(diff_cells, [])
    diff_vars = [var for var in diff_vars if var in df.columns]  # consider converting this row to validation method within metadatabundle and applying in all functions

    diff_like_vars = metadata.variable_groups.get(diff_like_cells, [])
    diff_like_vars = [var for var in diff_like_vars if var in df.columns]
    vars_to_convert = diff_vars + diff_like_vars

    # Identify rows diff_vars are not <=1
    all_above_1 = (df[diff_vars] < 1).all(axis=0).all()

    if all_above_1:
        df[vars_to_convert] = df[vars_to_convert] * 100
        print(f"[INFO] {diff_cells} values appeared to be fractions. Converted to percentages by multiplying by 100.")

    return df
