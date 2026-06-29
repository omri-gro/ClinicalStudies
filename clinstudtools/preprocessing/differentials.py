# clinstudtools/preprocessing/differentials.py
""" Quantitative Math Preprocessing Functions """

import pandas as pd
import numpy as np
from clinstudtools.core.metadata import MetadataBundle

def calc_diff(df, metadata, diff_cells="WBC diff", total_var="Total WBC", additional_cells=None, to_100=True):
    """
    Convert cell counts to cell differential, with
    Args:
        df (pd.DataFrame): The DataFrame to calculate for.
        metadata (MetadataBundle): Object with variable groups metadata.
        diff_cells (str or list): Group of cells included in the differential.
        total_var (str): Name of new column that will include the total of the diff_cells variables.
        additional_cells (str or list, optional): Group of cells to be converted to percentage according to the totalcells in diff_cells.
        to_100 (bool, optional): If true will calculate percetange values by 100 (so 100% is 100 instead of 1).
        # consider adding option when diff_cells are already percentages and only additional cells need conversion to percentages (given the total WBC variable).

    Returns:
        DataFrame: Modified DataFrame.
    """
    if additional_cells is None:
        additional_cells = []

    # Resolve diff_cells and additional_cells to lists of variable names
    if isinstance(diff_cells, str):
        diff_vars = metadata.variable_groups.get(diff_cells, [])
    else:
        diff_vars = diff_cells

    if isinstance(additional_cells, str):
        additional_vars = metadata.variable_groups.get(additional_cells, [])
    else:
        additional_vars = additional_cells

    # Check that all variables exist in the DataFrame
    diff_vars = [v for v in diff_vars if v in df.columns]
    additional_vars = [v for v in additional_vars if v in df.columns]

    if not diff_vars:
        print(f"Warning: No {diff_cells} variables found in DataFrame.")
        return df

    # Combine all variables used for validation
    all_vars = diff_vars + additional_vars

    # try converting to numeric
    df_orig = df
    df[all_vars] = df[all_vars].apply(pd.to_numeric, errors="coerce")

    # Check that all variables have exact counts (non-negative integers)
    invalid_vars = []
    for var in all_vars:
        if not pd.api.types.is_numeric_dtype(df[var]):
            invalid_vars.append(var)
            continue
        if not ((df[var].dropna() >= 0).all() and
                (df[var].dropna() == df[var].dropna().astype(int)).all()):
            invalid_vars.append(var)

    if invalid_vars:
        print(f"Warning: The following variables do not contain exact counts "
              f"(non-negative integers): {invalid_vars}. Aborting calculation.")
        return df_orig

    # Calculate total WBC count (sum of diff cells)
    df[total_var] = df[diff_vars].sum(axis=1, skipna=True)

    # Replace diff cells with percentages
    multiplier = 100 if to_100 else 1
    for var in diff_vars:
        df[var] = multiplier * df[var] / df[total_var]
        # print(f"Updated '{var}' to percentage of total WBC.")

    # Replace additional cells with percentages relative to total WBC
    for var in additional_vars:
        df[var] = multiplier * df[var] / df[total_var]
        # print(f"Updated '{var}' to percentage of total WBC (additional cell).")

    return df


def diff_from_total(df, metadata, diff_cells="WBC diff", total_count="Total WBC", to_100=True):
    """
    Convert cell counts to cell differential when given a TotalCount column.
    Number of cells counted for categories does not need to add up to the TotalCount
    (multi-classing and/or normals not included).
    Args:
        df (pd.DataFrame): The DataFrame to calculate for.
        metadata (MetadataBundle): Object with variable groups metadata.  # can change script to turn optional when diff_cells is list
        diff_cells (str or list): Group of cells included in the differential.
        total_count (str): Name of column to be used as denominator.
        to_100 (bool, optional): If true will calculate percetange values by 100 (so 100% is 100 instead of 1).

    Returns:
        DataFrame: Modified DataFrame.
    """
    # Resolve diff_cells and additional_cells to lists of variable names
    if isinstance(diff_cells, str):
        diff_vars = metadata.variable_groups.get(diff_cells, [])
    else:
        diff_vars = diff_cells

    # Check that all variables exist in the DataFrame
    diff_vars = [v for v in diff_vars if v in df.columns]

    if not diff_vars:
        print(f"\033[91mWarning: No {diff_cells} variables found in DataFrame.\033[0m")
        return df
    if total_count not in df.columns:
        print(f"\033[91mWarning: {total_count} variable not found in DataFrame.\033[0m")
        return df

    # for handling cases where total_count is <1
    no_total = ~(df[total_count] > 1)

    # Replace exact values of diff cells with percentages
    multiplier = 100 if to_100 else 1
    for var in diff_vars:
        df[var] = multiplier * df[var] / df[total_count]
        df.loc[no_total, var] = pd.NA

    return df


def check_diff_sum(df, metadata, diff_cells="WBC diff", tolerance=5, auto_convert=True, drop_flagged=True, force_num=True):
    """
    Check that differential variables (WBC, NDC, etc.) sum up to approximately 100%.
    If the sum is ~1 for all or nearly all rows and auto_convert=True, convert all diff_vars to percentage scale.

    Args:
        df (pd.DataFrame): The DataFrame to check.
        metadata (MetadataBundle): Object with variable groups metadata.
        diff_cells (str, optional): Name of variables group (appearing in metadata) to be summed up to 100%.
        tolerance (float, optional): Acceptable deviation from 100%. Default is ±5.
        auto_convert (bool, optional): Whether to automatically multiply by 100 when sums are ~1.
        drop_flagged (bool, optional): Whether to remove from dataframe the flagged (not adding up to 100%) rows.
        force_num (bool, optional): Whether to remove from dataframe rows where diff values are non-numeric (notice this will mess with rows that are fully grade-based).

    Prints:
        A warning for each sample where the sum of WBC diff variables
        falls outside (100 ± tolerance), unless they are NaN.

    Returns:
        df: possibly modified DataFrame (with scaled WBC diff variables)
    """
    # Get WBC diff variables from context metadata
    wbc_diff_vars = metadata.variable_groups.get(diff_cells, [])
    wbc_diff_vars = [var for var in wbc_diff_vars if var in df.columns]

    if not wbc_diff_vars:
        print(f"\033[91mNo {diff_cells} variables found in DataFrame.\033[0m")
        return df

    # Compute row-wise sums over diff

    # Treat empty strings / pure whitespace as NaN before coercion
    raw = df[wbc_diff_vars].replace(r"^\s*$", np.nan, regex=True)

    # check if values are numeric, and if not print which ones are strings
    # Coerce to numeric; anything non-numeric becomes NaN
    numeric_view = raw.apply(pd.to_numeric, errors="coerce")
    # Mask rows with any truly non-numeric values
    non_numeric_mask = numeric_view.isna() & raw.notna()
    offending_rows = non_numeric_mask.any(axis=1)
    offending_idx = df.index[offending_rows]
    # If all values are non-numeric
    if numeric_view.notna().sum().sum() == 0:
        print("WARNING: All values in selected columns are non-numeric. Nothing to sum.")
        return df, df.index, numeric_view
    # Print problematic rows
    if not offending_idx.empty:
        print("Dropping rows with non-numeric values (showing first 5):")
        print(raw.loc[offending_idx].head(5))
    # Keep safe rows only
    df = df.loc[~offending_rows]
    df[wbc_diff_vars] = df[wbc_diff_vars].apply(pd.to_numeric, errors="coerce")

    # Compute row-wise sum of WBC diff variables
    wbc_sum = df[wbc_diff_vars].sum(axis=1, skipna=True)


    # Identify rows where all WBC diff variables are NaN
    all_nan = df[wbc_diff_vars].isna().all(axis=1)

    # Auto-convert if nearly all rows are close to 1.0
    close_to_1 = wbc_sum.between(1 - tolerance / 100, 1 + tolerance / 100)
    n_valid = (~all_nan).sum()
    n_close_to_1 = close_to_1[~all_nan].sum()

    if auto_convert and n_valid > 0 and n_close_to_1 / n_valid > 0.9:
        # currently conversion applied automatically to all percent and percent-like variables
        # need to add section where absolute count of percent-like turned divided by wbc_sum before next step
        # need to add section where other diffs/grades (RBCs) are not included in next function if given in different format than WBC
        df = one_to_one_hundred(df, metadata, diff_cells='percent', diff_like_cells='percent-like')
        wbc_sum = wbc_sum * 100


    # Determine rows where sum is outside acceptable range
    lower_bound = 100 - tolerance
    upper_bound = 100 + tolerance
    out_of_range = ~wbc_sum.between(lower_bound, upper_bound)

    # Final flagged rows (excluding all-NaN rows)
    flagged = out_of_range & ~all_nan

    # if most are flagged, try calculating differential (numbers might be absolute counts)
    if flagged.sum() > flagged.size * 0.9:
        df = calc_diff(df, metadata, diff_cells="WBC diff", additional_cells="percent-like", to_100=True)


    # Report
    elif flagged.any():
        sample_ids = df.loc[flagged, "SampleID"]
        for sid, total in zip(sample_ids, wbc_sum[flagged]):
            if drop_flagged:
                print(f"\033[91mSample {sid} dropped: {diff_cells} sum = {total:.1f}%\033[0m")
            else:
                print(f"\033[91mSample {sid}: {diff_cells} sum = {total:.1f}% (expected ~100%)\033[0m")
        if drop_flagged:
            df = df[~flagged]
    else:
        print(f"All {diff_cells} sums within acceptable range.")

    return df



def one_to_one_hundred(df, metadata, diff_cells='percent', diff_like_cells='percent-like'):
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

