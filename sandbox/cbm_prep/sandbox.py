import pandas as pd
import numpy as np
import os
import re
import yaml
import matplotlib.pyplot as plt
from typing import List, Iterable, Mapping, Optional, Sequence, Union, Tuple
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

"""
creating derived_df that is a pivoted dataframe with row for each measurement (parameter, sample, site, method combination)
"""


class MetadataBundle:  # to do: add attribute for thresholds & grades?
    #  Holds all metadata and configuration for the pipeline.
    def __init__(self, metadata_path):
        # consider a simpler init without metadata_path, and move current content into from_yaml method
        context = load_yaml(metadata_path)
        self.variables = context.get("variables", {})
        self.alias_map = self.build_alias_map(self.variables)  # not the standard way of creating variables, could have just defined self.alias_map within the function
        self.variable_groups = self.build_variable_groups(self.variables)
        self.crit_points = self.build_lists_map(self.variables, keyword="crit_points")  # to do: add functionality where critical points are defined by grading thresholds if one is provided and the other isn't.
        self.normal_ranges = self.build_lists_map(self.variables, keyword="normal_range")  # notice this currently doesn't directly require exactly 2 values in normal_range
        self.grading_specs = self.building_grading_specs(self.variables)

        self.src_fixes = self.build_src_fixes(context)

        self.pregraded_index = self.build_pregraded_index()


    def build_alias_map(self, variables):
        alias_map = {}
        for canonical, props in variables.items():
            alias_map[canonical] = canonical  # map canonical to itself
            for alias in props.get("aliases", []):
                alias_map[alias] = canonical
        return alias_map

    def build_variable_groups(self, variables):
        """
        Precompute all variable groups as a dict: group_name → list of variable names.
        """
        group_map = {}
        for var, props in variables.items():
            for group in props.get("groups", []):
                group_map.setdefault(group, []).append(var)
        return group_map

    def build_lists_map(self, variables, keyword):
        """
        Create dictionary of all variables where the keyword holds a list (of critical points, ranges, etc.)
        """
        points_map = {}
        for var, props in variables.items():
            crit_points = props.get(keyword)
            if isinstance(crit_points, list):
                points_map[var] = crit_points
        return points_map

    def building_grading_specs(self, variables):
        """
        Create a dictionary of variable_name -> dict(thresholds, grades, right_closed, clamp_out_of_range),
        where only variables that have thresholds+grades appear.
        """
        grading_specs = {
            name: spec for name, spec in variables.items()
            if isinstance(spec, dict)
            and spec.get("thresholds") is not None
            and spec.get("grades") is not None
        }
        return grading_specs

    def build_src_fixes(self, context):
        """
        Create a dictionary of (site, method) -> {"rule name": ...}
        """
        fix_dict = context.get("Source fixes", {})
        src_fixes = {}
        for k, v in fix_dict.items():
            site, method = self.parse_site_method_key(k)
            src_fixes[(site, method)] = v or {}
        return src_fixes

    def build_pregraded_index(self):  # to do: consider applying dimensionless approach for indices as well
        """
        Create a Series where True if index (Site, Method, Variable) is provided as grade already
        Within the src_fixes attribute, this can be specified for variables or groups of variables
        """
        rows = []
        for (site, method), rule in self.src_fixes.items():
            for var in (rule.get("raw_given_as_grade") or []):
                if var in self.variables.keys():
                    rows.append((site, method, var))
                else:
                    for var_name in self.variable_groups.get(var, []):
                        rows.append((site, method, var_name))
        if rows:
            idx = pd.MultiIndex.from_tuples(rows, names=["Site", "Method", "Variable"])
            return pd.Series(True, index=idx, dtype=bool)
        else:
            return None

    @staticmethod
    def parse_site_method_key(key: str):
        if isinstance(key, tuple):
            site, method = key
        else:
            k = key.strip()
            if k.startswith("(") and k.endswith(")"):
                k = k[1:-1]
            site, method = [p.strip() for p in k.split(",", 1)]
        return site, method


"""General tools"""
def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def write_df_to_file(df: pd.DataFrame,
                     out_path: Union[str, Path]):
    """ Write dataframe to either Excel or csv"""
    format = out_path.split(".")[-1]
    if format.lower() == "csv":
        df.to_csv(out_path, index=False)
    elif format.lower() in ("xlsx", "excel"):
        df.to_excel(out_path, index=False)
    else:
        raise ValueError("Format must be 'csv' or 'excel'.")

def _as_df(obj_or_df) -> pd.DataFrame:
    """
        For functions that can treat either a MethodComparator or a DataFrame, except either.
    """
    df = getattr(obj_or_df, "df", None)   # change or set fallback for different attribute name ('data' instead of 'df')
    if isinstance(df, pd.DataFrame):
        return df.copy()
    if isinstance(obj_or_df, pd.DataFrame):
        return obj_or_df.copy()
    raise TypeError("Expected a MethodComparator or a pandas DataFrame.")


def _ensure_list(x):
    """If argument is a tuple/series, convert to list.
    If argument is anything else but a list, return list containing that one object."""
    if isinstance(x, list):
        return x
    elif isinstance(x, tuple):
        return list(x)
    elif isinstance(x, pd.Series):
        return x.to_list()
    else:
        return [x]

def _round_df(df: pd.DataFrame, decimals: int = 2) -> pd.DataFrame:
    out = df.copy()
    if isinstance(decimals, int):
        # global rounding of numerics
        out = out.apply(pd.to_numeric, errors="ignore")
        for c in out.columns:
            if pd.api.types.is_numeric_dtype(out[c]):
                out[c] = out[c].round(decimals)
    else:
        print(f'{decimals} not integer. Rounding not performed.')
    return out


def robust_dup(
        df: pd.DataFrame,
        key_cols: Sequence[str],
        on_duplicates: str = "error"  # "error" | "first" | "last" | "mean" | "median" | "sum" | "count" | "any" | "all"
) -> pd.DataFrame:
    """
    If duplicates exist according to key_cols, resolves them per on_duplicates.
    Notice functionality previously covered by aggfunc in pandas' pivot_table not added yet.
    """
    dup_mask = df.duplicated(subset=key_cols, keep=False)

    if not dup_mask.any():
        return df

    sort_col = ["Variable", "SampleID"] if "Variable" in df.columns and "SampleID" in df.columns else "n"
    dup_keys = (
        df.loc[dup_mask, key_cols]
        .value_counts()
        .rename("n")
        .reset_index()
        .sort_values(sort_col, ascending=False)
    )

    dup_report_str = "Duplicate entries found for pivot keys.\n" \
                     f"Examples:\n{dup_keys.head(5).to_string(index=False)}"

    if on_duplicates == "error":
        raise ValueError(dup_report_str)
    else:
        print(dup_report_str)

    if on_duplicates in {"first", "last"}:  # to do - add option for 'none' using keep=True in drop_duplicates
        return df.drop_duplicates(subset=key_cols, keep=on_duplicates)


def safe_pivot(
        df: pd.DataFrame,
        index: Sequence[str],
        columns: Sequence[str],
        values: str,
        sort_by: Optional[Sequence[str]] = None,   # e.g., ["ReviewDate","Investigator"]
        on_duplicates: str = "error",   # "error" | "first" | "last" | "mean" | "median" | "sum" | "count" | "any" | "all"
        ascending: Union[bool, Sequence[bool]] = True
) -> pd.DataFrame:
    """ Pivot with robust duplicate handling. """
    if sort_by:
        df = df.sort_values(sort_by, ascending=ascending, kind="mergesort")

    # duplicate probe
    key_cols = list(index) + list(columns)

    post_dup_df = robust_dup(df, key_cols, on_duplicates)

    """
    dup_mask = df.duplicated(subset=key_cols, keep=False)

    # replace this section with the new robust_dup function
    if not dup_mask.any():
        return df.pivot(index=index, columns=columns, values=values)

    dup_keys = (
        df.loc[dup_mask, key_cols]
        .value_counts()
        .rename("n")
        .reset_index()
        .sort_values("n", ascending=False)
    )

    dup_report_str = "Duplicate entries found for pivot keys.\n" \
                     f"Examples:\n{dup_keys.head(5).to_string(index=False)}"

    if on_duplicates == "error":
        raise ValueError(dup_report_str)
    else:
        print(dup_report_str)

    if on_duplicates in {"first", "last"}:
        # keep first/last row within each key group (deterministic due to pre-sort)
        dedup = (df.drop_duplicates(subset=key_cols, keep=on_duplicates))
        return dedup.pivot(index=index, columns=columns, values=values)

    return df.pivot_table(index=index, columns=columns, values=values, aggfunc=on_duplicates, observed=True)
    """
    if on_duplicates in {"error", "first", "last"}:
        return post_dup_df.pivot(index=index, columns=columns, values=values)
    else:
        return post_dup_df.pivot_table(index=index, columns=columns, values=values, aggfunc=on_duplicates, observed=True)


"""Data import and preparation tools"""
def read_to_df(file_name, sheet_name='Sheet1', file_dir=None, **kwargs):
    if file_dir == None:
        # dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        file_dir = os.path.abspath(os.path.dirname(__file__))
        file_dir = os.path.join(file_dir, r'raw')
    filepath = os.path.join(file_dir, file_name)
    if not os.path.exists(filepath):
        filepath = file_name
    elif not os.path.exists(filepath):
        raise FileNotFoundError(f"\033[91mFailed to find {file_name}\033[0m")
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filepath, sheet_name=sheet_name)
            print(f"Loaded Excel file: {filepath}")
        elif ext == '.csv':
            df = pd.read_csv(filepath, **kwargs)
            print(f"Loaded CSV file: {filepath}")
        else:
            raise ValueError(f"\033[93mUnsupported file format: {ext}\033[0m")
    except Exception as e:
        raise RuntimeError(f"\033[91mFailed to load {filepath}: {e}\033[0m")

    # Basic sanity check
    if df.empty:
        raise ValueError(f"\033[91mFile {filepath} is empty.\033[0m")

    print(f"Loaded {filepath} with shape {df.shape}")
    return df


def raw_to_df(file_name, site=None, method=None, sheet_name='Sheet1', dir=None):
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
        inv_str = ', '.join(df[inv_col].dropna().unique())
        num_inv = df[inv_col].nunique()
        print(f'{num_inv} investigators in dataframe: {inv_str}')
        df.rename(columns={inv_col: "Investigator"}, inplace=True)

    # removal of duplicates only if no investigator column (need to improve handling of investigator existence in main pipeline)
    if inv_col or site is None:  # consider adding option for removing multiple rows with same SampleID and investigator
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




def stnd_names(df, alias_map):
    rename_dict = {
        col: alias_map[col] for col in df.columns if col in alias_map
    }
    df = df.rename(columns=rename_dict)
    return df


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


def safe_row_sum(df: pd.DataFrame, cols: List[str], verbose: bool = True) -> Tuple[Optional[pd.Series], pd.Index, pd.DataFrame]:
    """
    Compute row-wise sums over `cols`, handling non-numeric values.

    - Converts values to numeric (treats empty strings/whitespace as NaN).
    - If all values are non-numeric, warns and returns None.
    - Otherwise, drops rows with truly non-numeric values (prints first 5 if verbose).
    - Returns:
        sum_series: pd.Series or None
        dropped_idx: index of rows dropped
        numeric_view: coerced numeric DataFrame of cols
    """
    if not cols:
        raise ValueError("`cols` must be a non-empty list of column names.")

    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Columns not in df: {missing_cols}")

    # Replace empty strings / whitespace with NaN before coercion
    raw = df[cols].replace(r"^\s*$", np.nan, regex=True)

    # Coerce to numeric
    numeric_view = raw.apply(pd.to_numeric, errors="coerce")

    # Mask rows with any truly non-numeric values
    non_numeric_mask = numeric_view.isna() & raw.notna()
    offending_rows = non_numeric_mask.any(axis=1)
    offending_idx = df.index[offending_rows]

    # If all values are non-numeric
    if numeric_view.notna().sum().sum() == 0:
        if verbose:
            print("WARNING: All values in selected columns are non-numeric. Nothing to sum.")
        return None, df.index, numeric_view

    # Print problematic rows (first 5 only)
    if not offending_idx.empty and verbose:
        print("Dropping rows with non-numeric values (showing first 5):")
        print(raw.loc[offending_idx].head(5))

    # Keep safe rows only
    safe_numeric = numeric_view.loc[~offending_rows]

    # Compute row-wise sums
    sum_series = safe_numeric.sum(axis=1, skipna=True)

    return sum_series, offending_idx, numeric_view


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
    diff_vars = [var for var in diff_vars if var in df.columns]

    diff_like_vars = metadata.variable_groups.get(diff_like_cells, [])
    diff_like_vars = [var for var in diff_like_vars if var in df.columns]
    vars_to_convert = diff_vars + diff_like_vars

    # Identify rows diff_vars are not <=1
    all_above_1 = (df[diff_vars] < 1).all(axis=0).all()

    if all_above_1:
        df[vars_to_convert] = df[vars_to_convert] * 100
        print(f"[INFO] {diff_cells} values appeared to be fractions. Converted to percentages by multiplying by 100.")

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
        df = one_to_one_hundered(df, metadata, diff_cells='percent', diff_like_cells='percent-like')
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


def calc_diff(df, metadata, diff_cells="WBC diff", additional_cells=None, to_100=True):
    """
    Convert cell counts to cell differential, with
    Args:
        df (pd.DataFrame): The DataFrame to calculate for.
        metadata (MetadataBundle): Object with variable groups metadata.
        diff_cells (str or list): Group of cells included in the differential.
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
    df["_total_WBC"] = df[diff_vars].sum(axis=1, skipna=True)

    # Replace diff cells with percentages
    multiplier = 100 if to_100 else 1
    for var in diff_vars:
        df[var] = multiplier * df[var] / df["_total_WBC"]
        print(f"Updated '{var}' to percentage of total WBC.")

    # Replace additional cells with percentages relative to total WBC
    for var in additional_vars:
        df[var] = multiplier * df[var] / df["_total_WBC"]
        print(f"Updated '{var}' to percentage of total WBC (additional cell).")

    # Drop helper column
    df.drop(columns="_total_WBC", inplace=True)

    return df

def diff_from_total(df, metadata, diff_cells="WBC diff", total_count="TotalWBC", to_100=True):
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


# to do: add step/function before/after this function where other shapes of same grade are converted to uniform look (e.g., "x" changed to "Negative")
# to do: integrate this somehow with checking positive cases (will behave like grade with single threshold, or checking where not nan/0 when pregraded)
def add_grade_column(df_long: pd.DataFrame, meta: "MetadataBundle"):
    """
    Expects df_long columns: SampleID, Variable, Value, Method, Site (case-sensitive).
    Adds column for grades based on number in value column.
    Uses MetaDataBundle, which is expected to include grading_specs dict with thresholds and grades for variable names.
    """
    # currently data is converted to grades but not to pd.CategoricalDtype - if want to change add the _coerce_grade_categorical function

    df = df_long.copy()
    req = {"Variable", "Value", "Method", "Site"}
    if not req.issubset(df.columns):
        raise ValueError(f"Expected columns: {req}")

    df["Grade"] = pd.NA
    df["Grade_from"] = pd.NA

    # Copy-as-is for non-numeric rows
    value_num = pd.to_numeric(df["Value"], errors="coerce")
    non_numeric = value_num.isna() & df["Value"].notna()  # strings like "N/A", "Not evaluable", etc.
    if non_numeric.any():
        df.loc[non_numeric, "Grade"] = df.loc[non_numeric, "Value"]
        df.loc[non_numeric, "Grade_from"] = "provided"
        # in these cases, the "values" were really just grades, so no need to include them in continuous values analysis
        df.loc[non_numeric, "Value"] = pd.NA


    # mask indices where MetaDataBundle claims raw was given as grade
    if meta.pregraded_index is not None:
        keys = pd.MultiIndex.from_frame(df[["Site", "Method", "Variable"]])
        pregraded = meta.pregraded_index.reindex(keys, fill_value=False).to_numpy()
    else:
        pregraded = np.zeros(len(df), dtype=bool)

    # copy when source says raw values are already grades  (notice that currently nan of some types are also copied if belonging to these variables)
    mask_pregraded_numeric = pregraded & ~non_numeric
    if mask_pregraded_numeric.any():
        df.loc[mask_pregraded_numeric, "Grade"] = df.loc[mask_pregraded_numeric, "Value"]
        df.loc[mask_pregraded_numeric, "Grade_from"] = "provided"
        # in these cases, the "values" were really just grades, so no need to include them in continuous values analysis
        df.loc[mask_pregraded_numeric, "Value"] = pd.NA


    # derive grade from gradable numeric values (which were not already provided as grades)
    gradable_vars = set(meta.grading_specs.keys())
    need_convert = (~pregraded) & ~non_numeric & df["Variable"].isin(gradable_vars)
    if need_convert.any():
        for var in sorted(gradable_vars):
            mv = need_convert & df["Variable"].eq(var)
            if not mv.any():
                continue
            spec = meta.grading_specs[var]
            try:
                df.loc[mv, "Grade"] = cut_series_to_categorical(
                    df.loc[mv, "Value"].astype(float),
                    thresholds=spec["thresholds"],
                    grades=spec["grades"],
                    right_closed=spec.get("right_closed", True),
                    clamp_out_of_range=spec.get("clamp_out_of_range", True))
                df.loc[mv, "Grade_from"] = "derived"
            except ValueError as e:
                print(f"\033[93mGrading specs error for {var}: {e}. Values not converted to grades.\033[0m")
                df.loc[mv, "Grade_from"] = "conversion error"
    return df


def add_pos_column(df_long: pd.DataFrame, meta: "MetadataBundle",
                   normal_grades=[0, "0", "Normal", "Negative", "normal", "negative"]):
    """
    Adds boolean column for positivity based on normal ranges in MetaDataBundle.
    If grade already exists, treat values in normal_vals as False (negative) and rest as True,
    Else use normal ranges and the "Value" column if possible.
    """
    # create empty boolean column
    df = df_long.copy()
    df["Positive"] = np.nan
    df["Positive"] = df["Positive"].astype('boolean')

    # where grade exists, use it for positivity  -  this section might need changing if grade ever not related to positivity
    grade_negative = df["Grade"].isin(normal_grades)
    grade_positive = df["Grade"].notna() & ~df["Grade"].isin(normal_grades)
    df.loc[grade_negative, "Positive"] = False
    df.loc[grade_positive, "Positive"] = True

    # find values where positivity will be based on normal ranges
    value_num = pd.to_numeric(df["Value"], errors="coerce")
    numeric = value_num.notna()
    normal_ranges = getattr(meta, "normal_ranges", {})
    norm_range_vars = set(normal_ranges.keys())
    need_convert = df["Positive"].isna() & numeric & df["Variable"].isin(norm_range_vars)

    # convert based on normal ranges
    if need_convert.any():
        for var in sorted(norm_range_vars):
            mv = need_convert & df["Variable"].eq(var)
            if not mv.any():
                continue
            norm_range = normal_ranges[var]
            if len(norm_range) == 2:
                df.loc[mv, "Positive"] = ~df.loc[mv, "Value"].between(norm_range[0], norm_range[1], inclusive="both")
            else:
                print(f'{norm_range} is not an appropriate normal range for {var}')
    return df


def cut_series_to_categorical(x: pd.Series,
                              thresholds,
                              grades,
                              *,
                              right_closed: bool = True,
                              clamp_out_of_range: bool = False) -> pd.Categorical:
    """
    Bin continuous values into grades.

    Args:
        x: numeric Series.
        thresholds: list-like of bin edges (e.g., [0,5,10,20,101]).
        grades: list-like of labels (e.g., [0,1,2,3] or ["none","mild","mod","sev"]).
        right_closed: if True, bins are right-closed ( (a,b] ).
        clamp_out_of_range: clip x to [min(thresholds), max(thresholds)] before cutting. Only used when len(thresholds) == len(grades) + 1.

    - If len(thresholds) == len(grades) - 1: treated as *interior cut points*; we pad with -inf/+inf
    - If len(thresholds) == len(grades) + 1: treated as full bin edges; must already bracket all bins.

    Returns:
        pandas Series with categories exactly as in `grades`.
    """
    y = x.astype(float)
    if len(thresholds) == len(grades) - 1:
        edges = np.concatenate(([-np.inf], thresholds, [np.inf]))
    elif len(thresholds) == len(grades) + 1:
        edges = thresholds
        if clamp_out_of_range:
            y = y.clip(lower=edges[0], upper=edges[-1])
    else:
        raise ValueError(f"\033[93mGrades must have one more or one less member than thresholds\033[0m")

    # Validate monotonicity (required by pandas.cut)
    if not np.all(np.diff(edges) > 0):
        raise ValueError("thresholds/edges must be strictly increasing.")

    out = pd.cut(y, bins=edges, labels=grades, right=right_closed, include_lowest=True)
    return out


def curate_df(df, metadata, src=None, wbcs_as_counts=False):
    # to do: check the type of results' source (scopio, OMR from specific site, CRF, etc.)
    if not src:
        src = (df["Site"][0], df["Method"][0])

    # standardize column names
    df = stnd_names(df, metadata.alias_map)

    # apply source-specific fixes
    df = apply_src_fixes(df, metadata, src)

    if wbcs_as_counts:
        # convert wbc and wbc-like variables into percetnages
        df = calc_diff(df, metadata, diff_cells="WBC diff", additional_cells="WBC-like")
    else:
        # print warning if WBCs in differential don't add up to ~100
        df = check_diff_sum(df, metadata, tolerance=5)
        df = check_diff_sum(df, metadata, tolerance=5)

    return df

def add_mean_investigator(df, mthd='ClV', min_inv=0, mean_inv_name="Mean Investigator"):
    # for the observations in mthd, calculate the average of all investigator's reviews
    # observations with less than min_inv investigators' reviews will not have mean calculated
    # returns modified df

    val = "Value"
    subset = ['SampleID', 'Site', 'Variable']
    inv_subset = subset + ['Investigator']

    # take only rows to be used for mean calculation
    invs_df = df.query(f"Method=='{mthd}' and {val}.notna()")

    # only rows with numeric values
    invs_df[val] = pd.to_numeric(invs_df[val], errors='coerce')
    invs_df = invs_df[invs_df[val].notnull()]

    # error if same investigator reviewed same sample twice
    invs_df = robust_dup(invs_df, key_cols=inv_subset, on_duplicates='error')

    if min_inv:
        # check number of investigators that reviewed each sample
        inv_counts = invs_df.groupby(subset)["SampleID"].transform('count')
        mask_keep = inv_counts >= min_inv

        # report on dropped samples
        dropped = invs_df.loc[~mask_keep, ['SampleID', 'Site']].drop_duplicates()
        dropped_report_str = f"Removed {dropped.shape[0]} samples with <{min_inv} investigators' reviews.\n" \
                             f"Examples:\n{dropped.head(5).to_string(index=False)}"
        print(dropped_report_str)

        invs_df = invs_df.loc[mask_keep].copy()

    # Compute mean for each group
    means_df = invs_df.groupby(subset, as_index=False)[val].mean()

    # add to original
    means_df['Investigator'] = mean_inv_name
    means_df['Method'] = mthd
    means_df['ValueOrigin'] = 'Mean'

    return pd.concat([df, means_df])


def pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName"], value_vars=[]):
    """
    Convert wide DataFrame to long format.

    Args:
        df (DataFrame): Wide format DataFrame.
        id_vars (list): Columns to keep as identifiers.
        piv_vars (list): Variables to get rows in the long format DataFrame (if empty, take all variables not in id_vars)

    Returns:
        DataFrame: Long format DataFrame.
    """
    if value_vars == []:
        value_vars = [col for col in df.columns if col not in id_vars]
    else:
        value_vars = [col for col in df.columns if col not in id_vars and col in value_vars]

    long_df = pd.melt(
        df,
        id_vars=id_vars,
        value_vars=value_vars,
        var_name="Variable",
        value_name="Value"
    )
    return long_df


def create_derived_variables_long(df, metadata):
    """
    Create derived variables in a long-format DataFrame using metadata.

    Args:
        df (pd.DataFrame): Long-format DataFrame (concatenated from all sites).
        metadata (MetadataBundle): Metadata containing derived variable recipes.

    Returns:
        pd.DataFrame: DataFrame with derived variable rows appended.
    """
    derived_rows = []

    # Possible grouping columns to preserve context
    possible_group_cols = ["SampleID", "Site", "Method", "Investigator", "Unit"]

    for derived_var, props in metadata.variables.items():
        derived_from = props.get("derived_from")
        if not derived_from:
            continue  # Skip if no derivation defined

        components = derived_from.get("components", [])
        operation = derived_from.get("operation", "sum")

        # Filter to component rows
        comp_df = df[df["Variable"].isin(components)].copy()
        if len(comp_df) == 0:
            print(f'Dataframe does not include any of {components}. The derived variable {derived_var} will not be calculated.')
            continue
        comp_df["Value"] = pd.to_numeric(comp_df["Value"], errors='raise')

        # Group by sample and site identifiers
        group_cols = [col for col in possible_group_cols if col in df.columns]

        aggregated = comp_df.groupby(group_cols, dropna=False).agg({"Value": operation})

        # Build derived rows
        aggregated = aggregated.reset_index()
        aggregated["Variable"] = derived_var
        aggregated["ValueOrigin"] = "Derived"

        derived_rows.append(aggregated)

    if derived_rows:
        derived_df = pd.concat(derived_rows, ignore_index=True)
        df = pd.concat([df, derived_df], ignore_index=True)
    else:
        print("No derived variables created (no recipes in metadata).")

    return df


""" functions mostly used as MethodComparator helpers but can also be used independently"""
def long_pivoted_to_side_by_side(df: pd.DataFrame, comp_cols='Method', remain_index='SampleID', values='Measurment'):
    """
    Convert a pivoted row-per-measurement single-column-for-all-numbers dataframe to a readable single-sample-per-row.
    Args:
        df (pd.DataFrame): Long-format DataFrame (concatenated from all sites)
        comp_cols (str or list, optional): Variable or variables that will now be column names (values in last variable name in list will be side-by-side).
        remain_index (str or list, optional): Row identifiers for new dataframe (old column names that stay as before)
        values (str or list, optional): Current column with numbers that will appear in the table (if list, they will be shown side-by-side)
    Returns:
        pd.Dataframe
    """
    # to do: add checks of inputs? add print of new table?
    return df.pivot(index=remain_index, columns=comp_cols, values=values)


def _default_col_order(df: pd.DataFrame, comparison_dims) -> list[list]:
    orders = []
    for col in comparison_dims:
        if col in df.columns:
            vals = df[col].dropna().unique()
            # sort values, but preserve natural ordering if mixed types
            try:
                vals = sorted(vals)
            except Exception:
                vals = list(vals)  # fall back to insertion order if un-sortable
            orders.append(list(vals))
        else:
            orders.append([])
    return orders


def to_comparison_matrix(
        obj_or_df: Union["MethodComparator", pd.DataFrame],
        metadata=None,
        # identifiers that define a datapoint row (only the ones present will be used) - remember to move "Investigator" to here if each investigator in own row
        row_identifiers: Sequence[str] = ("SampleID", "Site"),
        # which columns define the side-by-side comparison blocks
        comparison_dims: Sequence[str] = ("Variable", "Method", "Investigator"),
        value_col: str = "Value",  # or "Grade"
        needed_vars: Optional[Iterable[str]] = None,   # if None and metadata provided, pull from metadata.variable_groups['percent']
        require_complete_cases: bool = True,           # only keep rows where all comparison cells are present (no NaNs)
        drop_na_mode: str = "all",  # or "any" if even a single NaN in row is enough to drop it
        decimals: int = 3,   # will attempt to round any output column possible. Change to non-integer to avoid rounding.
        column_order: Optional[Sequence[Sequence]] = None,  # e.g., [list_of_methods, list_of_investigators]; order applied where dims exist
        flatten_columns: bool = True,  # flatten MultiIndex columns like ('MethodA','Inv1') -> "MethodA|Inv1"
        sep: str = "|",  # for column name flattening
        # arguments for handling duplicates (see safe_pivot)
        on_duplicates: str = "error",
        sort_by: Optional[Sequence[str]] = None,   # e.g., ["ReviewDate","Investigator"]
        ascending: Union[bool, Sequence[bool]] = True,
        **kwargs
) -> pd.DataFrame:
    """
        Build and save a wide, readable matrix of only the datapoints actually used in comparisons.

        Returns the wide DataFrame (also saved to Excel/CSV per out_path).
    """
    df = _as_df(obj_or_df)

    # --- choose variables ---
    if needed_vars is None and metadata is not None:
        # future adjustment - change tuple to order of preference on variables group names
        for key in ("Evaluation parameter", "percent"):  # future adjustment - provide choice of non-percent if value_col=Grade or specifically requested
            if getattr(metadata, "variable_groups", None) and key in metadata.variable_groups:
                needed_vars = metadata.variable_groups[key]
                break
    if needed_vars is not None:
        df = df[df["Variable"].isin(needed_vars)]

    # --- keep only columns that actually exist ---
    present_identifiers = [c for c in row_identifiers if c in df.columns]
    present_comp_dims = [c for c in comparison_dims if c in df.columns]
    if not present_identifiers:
        raise ValueError(f"No the columns {row_identifiers} were found in the data.")
    if not present_comp_dims:
        raise ValueError(f"None of the columns {comparison_dims} exist in the data.")

    # --- pivot to side-by-side ---
    wide = safe_pivot(df, index=present_identifiers, columns=present_comp_dims, values=value_col,
                      on_duplicates=on_duplicates, sort_by=sort_by, ascending=ascending)
    wide = df.pivot(index=present_identifiers, columns=present_comp_dims, values=value_col)

    # Name the row index levels using the row identifiers used
    wide.index.names = list(present_identifiers)

    # --- enforce complete cases if requested (ensure true ‘used’ datapoints only) ---
    # to do: correct require_complete_cases so it requires data in all methods (getting nans in all columns of a single method is enough for dropping) - for simplicity, use a single column that must exist in all methods and sites for this filtering
    if require_complete_cases:
        # drop rows with any missing across all comparison cells
        if isinstance(wide.columns, pd.MultiIndex):
            subset_cols = wide.columns
        else:
            subset_cols = list(wide.columns)
        wide = wide.dropna(axis=0, how=drop_na_mode, subset=subset_cols)


    # --- column reordering by provided order ---
    if not column_order:
        column_order = _default_col_order(df, present_comp_dims)

    # Build a MultiIndex product for the dims we actually have.
    # If user supplied lengths mismatch (e.g., only methods list but no Investigator in data), use what exists.
    # Example: column_order = [methods, investigators]
    # If investigators aren't present, only use 'methods' order.
    usable_orders = []
    for dim_name, order in zip(present_comp_dims, column_order):
        # Only include if this dim is present
        if dim_name in present_comp_dims and order is not None:
            usable_orders.append(order)

    # Reindex by product when we still have at least 1 reorder list
    if usable_orders:
        try:
            if len(present_comp_dims) == 1:
                # Single level columns
                new_cols = pd.Index(usable_orders[0], name=present_comp_dims[0])
            else:
                # Multi level columns
                new_cols = pd.MultiIndex.from_product(usable_orders, names=list(present_comp_dims))
            wide = wide.reindex(columns=new_cols)
        except Exception:
            # keep current order if mismatch
            pass

    wide = _round_df(wide, decimals)

    # --- optionally flatten columns for friendlier Excel/CSV ---
    if flatten_columns and isinstance(wide.columns, pd.MultiIndex):
        wide.columns = [sep.join(map(str, col)).strip(sep) for col in wide.columns]
        wide.columns = pd.Index(wide.columns, name=sep.join(present_comp_dims))

    # drop empty columns - turn into boolean argument based if needed
    wide = wide.dropna(axis=1, how='all')

    return wide.reset_index()


"""Statistical tools"""
def calc_bias_at_points(reg_res, crit_points: List):
    """
    Calculate bias for specific points on a single regression line.
    Args:
        reg_res: a RegressionResult object
        crit_points: a list of integers and/or floats

    Returns:
        dict: dict with the ints/floats being keys, and values being dicts which include the keys abs_bias, abs_bias_ci, rel_bias, rel_bias_ci
    """
    if not isinstance(crit_points, list):
        return TypeError('crit_points must be list of numbers')
    bias_dict = {}
    def bias_at_point(cp, slope, intercept):
        obs_val = cp * slope + intercept
        abs_err = obs_val - cp
        rel_err = np.nan if cp == 0 else abs_err / (cp / 100)
        return abs_err, rel_err
    for point in crit_points:
        abs_err, rel_err = bias_at_point(point, reg_res.slope, reg_res.intercept)
        abs_err_btm, rel_err_btm = bias_at_point(point, reg_res.slope_ci[0], reg_res.intercept_ci[0])
        abs_err_top, rel_err_top = bias_at_point(point, reg_res.slope_ci[1], reg_res.intercept_ci[1])
        sngl_point_dict = {
            "abs_bias": abs_err,
            "abs_bias_ci": (abs_err_btm, abs_err_top),
            "rel_bias": rel_err,
            "rel_bias_ci": (rel_err_btm, rel_err_top)
        }
        bias_dict[point] = sngl_point_dict
    return bias_dict



"""Plotting tools - need to move them to their own script file"""
def plot_ver_reg(x, y, cls_var=None, reg_ser=None, eq_line=True, fig=None, ax=None):
    if ax is None or fig is None:
        fig, ax = plt.subplots()
    ax.plot(x, y, marker='o', linestyle='', markersize=4)
    if isinstance(reg_ser, pd.Series) or isinstance(reg_ser, dict):
        # Plot the regression line

        x_range = np.linspace(min(x) * 0.9, max(x) * 1.1, 100)
        y_line = reg_ser['slope'] * x_range + reg_ser['intercept']
        ax.plot(x_range, y_line, label="Deming Regression Line", color="red")
        #
        # Plot the confidence interval
        y_low = reg_ser['slope_ci_bottom'] * x_range + reg_ser['intercept_ci_bottom']
        y_top = reg_ser['slope_ci_top'] * x_range + reg_ser['intercept_ci_top']
        ax.fill_between(x_range, y_low, y_top, color="pink", alpha=0.3, label="95% Confidence Interval")
    if eq_line:
        # plot y=x line (perfect agreement) for comparison
        line = np.arange(0, 100, 1)
        ax.plot(line, line, color="grey")
    x_max = max(x)
    y_max = max(y)
    plt.xlim(-0.05 * x_max, 1.05 * x_max)
    plt.ylim(-0.05 * y_max, 1.05 * y_max)
    # tick_spacing = 5
    # ax.xaxis.set_major_locator(MultipleLocator(tick_spacing))
    # ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
    # ax.legend()
    ax.grid(True)
    return fig, ax


def set_equal_limits_and_scale(fig=None, ax=None):
    # to do - remove, as already exists in plotting.py
    """
    Sets the x and y axes to have:
    - The same scale (1 unit in x = 1 unit in y)
    - The same limits (xmin == ymin, xmax == ymax)

    Parameters:
    - fig: matplotlib.figure.Figure (optional)
    - ax: matplotlib.axes.Axes (optional)

    If both are None, uses the current active figure.
    """
    import matplotlib.pyplot as plt

    if ax is not None:
        axes = [ax]
    elif fig is not None:
        axes = fig.get_axes()
    else:
        axes = plt.gcf().get_axes()

    for a in axes:
        x_min, x_max = a.get_xlim()
        y_min, y_max = a.get_ylim()

        data_min = min(x_min, y_min)
        data_max = max(x_max, y_max)

        a.set_xlim(data_min, data_max)
        a.set_ylim(data_min, data_max)
        a.set_aspect('equal')


""" BMA study tools """
def raw_bma_to_df(file_name, site, method, sheet_name='Sheet1', dir=None):
    # temporary form, works when semi-automatic data prep used
    # to do: create pointer function that chooses which raw-to-df function to use
    df = read_to_df(file_name, sheet_name, dir)
    # this entire section plus the header=None can be its own transpose with multiple identifiers function
    df = pd.read_csv(os.path.join(dir, file_name), header=None)
    sample_ids = df.iloc[0]
    reviewers = df.iloc[1]
    df.columns = pd.MultiIndex.from_arrays([sample_ids, reviewers])
    df = df.iloc[2:].copy()
    df = df.T.reset_index()
    df.columns = df.iloc[0]
    df = df[1:]

    # Standardize the SampleID
    possible_id_cols = ["SampleID", "Sample", "Sample ID", "ID", "Barcode", "barcode", "Case", "Anonymised no.", "Case ID"]
    id_col = next((col for col in possible_id_cols if col in df.columns), None)
    if not id_col:
        raise ValueError(f"No sample ID column found in {file_name}")
    df.rename(columns={id_col: "SampleID", "Signed Off By": "Investigator"}, inplace=True)

    if "SampleID" not in df.columns:
        raise ValueError(f"Missing 'SampleID' column in {file_name}")

    # Add metadata columns
    df["Site"] = site
    df["Method"] = method
    df["FileName"] = os.path.basename(file_name)

    # Strip leading/trailing spaces from column names
    df.columns = df.columns.str.strip()

    # to do: lowercase column names?
    return df


