# transforms.py - functions reshape or filter dataframe generically

import pandas as pd
from typing import List, Iterable, Mapping, Optional, Sequence, Union, Tuple
import numpy as np




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


def robust_dup(
        df: pd.DataFrame,
        key_cols: Sequence[str],
        on_duplicates: str = "error",  # "error" | "first" | "last" | "mean" | "median" | "sum" | "count" | "any" | "all",
        sort_cols: Union[Sequence[str], None] = None,  # domain-specific order for printing duplicates
) -> pd.DataFrame:
    """
    If duplicates exist according to key_cols, resolves them per on_duplicates.
    Notice functionality previously covered by aggfunc in pandas' pivot_table not added yet.
    """
    dup_mask = df.duplicated(subset=key_cols, keep=False)

    if not dup_mask.any():
        return df

    # If duplicates exist according to key_cols, resolves them per on_duplicates.
    print_order = sort_cols or ["n"] + key_cols
    dup_keys = (
        df.loc[dup_mask, key_cols]
        .value_counts()
        .rename("n")
        .reset_index()
        .sort_values(print_order, ascending=False)
    )

    dup_report_str = "Duplicate entries found for pivot keys.\n" \
                     f"Examples:\n{dup_keys.head(5).to_string(index=False)}"

    if on_duplicates == "error":
        raise ValueError(dup_report_str)
    else:
        print(dup_report_str)

    if on_duplicates in {"first", "last"}:  # to do - add option for 'none' using keep=True in drop_duplicates
        return df.drop_duplicates(subset=key_cols, keep=on_duplicates)

# to do: find a way to merge this and previous function (& refactor legacy calls)
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

    if on_duplicates in {"error", "first", "last"}:
        return post_dup_df.pivot(index=index, columns=columns, values=values)
    else:
        return post_dup_df.pivot_table(index=index, columns=columns, values=values, aggfunc=on_duplicates, observed=True)


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
