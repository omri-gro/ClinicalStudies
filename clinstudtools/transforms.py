# transforms.py - functions reshape or filter dataframe generically

import pandas as pd
from typing import List, Iterable, Mapping, Optional, Sequence, Union, Tuple
import numpy as np
from .utils.coercion import ensure_list
from .utils.io import read_to_df



def _round_df(df: pd.DataFrame, decimals: int = 2) -> pd.DataFrame:
    """
        Round all numeric or numeric-like columns by N decimal places.
    """
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


def filter_by_reference(
        df: pd.DataFrame,
        filtering_source: Union[str, pd.DataFrame],
        filtering_cols: Optional[Sequence[str]] = None,
        include_rows: bool = False,
        target_vars: Optional[Union[str, Sequence[str]]] = None,
        metadata=None
) -> pd.DataFrame:
    """
    Filter a DataFrame based on a reference list of samples to exclude/include.
    filtering_cols - Sequence of column names to use for filtering (e.g., ['Site', 'SampleID', 'Investigator', 'Method']),
                         else infer from the filtering df.
    include_rows - Change to True so resulting df will have only samples included in the filtering_source, instead of those that don't  (overwrites both).
    target_vars - If provided, filtering only drops/keeps rows for these specific variables.
                      Data for other variables in the same sample is left untouched.
    metadata - Used to unpack variable groups (e.g., 'WBC diff') if target_vars is provided.
    """
    # Use your existing safe read function
    df2 = filtering_source if isinstance(filtering_source, pd.DataFrame) else read_to_df(filtering_source)

    if df2.empty:
        print('\033[34mNothing to filter out!\033[0m')
        return df.copy()

    if not filtering_cols:
        filtering_cols = set(df2.columns)
        filtering_cols.discard("FileName")

    df1 = df.copy()
    common = list(set(df1.columns) & set(filtering_cols))

    df1_com = df1[common]
    df2_com = df2[common]

    # Handle SampleID type mismatches
    if pd.api.types.is_numeric_dtype(df2_com['SampleID'].dtype) and not pd.api.types.is_numeric_dtype(
            df1_com['SampleID'].dtype):
        df2_com = df2_com[df2_com["SampleID"].notna()]
        num0s = len(df1_com["SampleID"].astype(str).iloc[0])
        ids_in_format = df2_com["SampleID"].astype(int).astype(str).str.zfill(num0s)
        df2_com = df2_com.assign(SampleID=ids_in_format)

    # Build multi-index mask
    keys_df2 = pd.MultiIndex.from_frame(df2_com.drop_duplicates())
    mask = pd.MultiIndex.from_frame(df1_com).isin(keys_df2)

    # Target variable filtering
    if target_vars is not None:
        target_vars = ensure_list(target_vars)
        if metadata:
            expanded_vars = []
            for tv in target_vars:
                expanded_vars.extend(metadata.variable_groups.get(tv, [tv]))
            target_vars = expanded_vars

        is_target_var = df1['Variable'].isin(target_vars)
        final_mask = ~(is_target_var & ~mask) if include_rows else ~(is_target_var & mask)
        return df1.loc[final_mask].copy()

    # Global filtering
    mask = mask if include_rows else ~mask
    return df1.loc[mask].copy()


def filter_by_condition(df: pd.DataFrame, condition: str) -> pd.DataFrame:
    """Filter cases for which a condition (on specific variables) is true."""
    try:
        return df.query(condition).copy()
    except (SyntaxError, ValueError) as e:
        raise ValueError(f"\033[91mQuery expression {condition} raised error: {e}") from e


def filter_samples_by_condition(df: pd.DataFrame, condition: str, filtering_cols=None) -> pd.DataFrame:
    """
    Keep all rows for samples (e.g., SampleID + Site) where a specific condition is met.

    Example:
        filter_samples_by_condition(df, "Variable == 'Total WBC' and Value >= 100")
    """
    needed_cases = filter_by_condition(df, condition)
    return filter_by_reference(df, needed_cases, filtering_cols, include_rows=True)

