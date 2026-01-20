from typing import Tuple, Optional, Union, Sequence
import pandas as pd


def safe_pivot(df: pd.DataFrame,
               index: Sequence[str],
               columns: str,
               values: str,
               on_duplicates: str = "raise",
               sort_by: Optional[Sequence[str]] = None,
               ascending: Union[bool, Sequence[bool]] = True,
               ) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """ Pivot with explicit, deterministic duplicate handling. """
    if sort_by:
        df = df.sort_values(sort_by, ascending=ascending, kind="mergesort")

    # Verify uniqueness per ID, and dimension
    if isinstance(columns, str):
        columns = [columns]
    key_cols = list(index) + columns
    df_clean = robust_dup(df, key_cols=key_cols, on_duplicates=on_duplicates)

    wide = df_clean.pivot(
        index=index,
        columns=columns,
        values=values,
    )
    return wide


def robust_dup(df: pd.DataFrame,
               key_cols: Sequence[str],
               on_duplicates: str = "raise",  # "raise" | "first" | "last" | "none"
               ) -> pd.DataFrame:
    """ Resolve duplicate rows according to key_cols. """
    if on_duplicates == "error":   # backwards compatibility
        on_duplicates = "raise"
    if on_duplicates not in {"raise", "first", "last", "none"}:
        raise ValueError(f"Invalid on_duplicates option: {on_duplicates}")

    report = dup_report(df, key_cols)
    if report is None or report.empty:
        return df

    preview = report.head(10).to_string(index=False)
    n_dup_keys = len(report)
    dup_str = f"[dedup] on_duplicates='{on_duplicates}': {n_dup_keys} duplicate keys detected.\nExamples:\n{preview}"

    if on_duplicates == "raise":
        raise ValueError(dup_str)

    if on_duplicates in {"first", "last"}:
        df_out = df.drop_duplicates(subset=key_cols, keep=on_duplicates)


    elif on_duplicates == "none":
        # drop all rows whose key appears in dup_keys
        dup_keys = report[key_cols]
        df_out = df.merge(dup_keys.assign(_dup=1), on=key_cols, how="left")
        df_out = df_out[df_out["_dup"].isna()].drop(columns="_dup")

    # inform on dropped rows
    rows_before = len(df)
    rows_after = len(df_out)
    rows_dropped = rows_before - rows_after
    print(f"{dup_str}\nDropped {rows_dropped} rows.")

    return df_out


def dup_report(df: pd.DataFrame,
               key_cols: Sequence[str],
               ) -> Optional[pd.DataFrame]:
    """
    Return keys with counts for rows that are duplicates across all key_cols, or None if no duplicates exist.
    key_cols should be names of columns in df.
    """
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    if not dup_mask.any():
        return None

    return (df.loc[dup_mask, key_cols]
            .value_counts(dropna=False)
            .rename("count")
            .reset_index()
            .sort_values("count", ascending=False))


def filter_to_ids(df, ids, id_cols):
    idx = df.set_index(list(id_cols)).index
    return df[idx.isin(ids)]


# can add safe_pivot_agg later (likely to different module) if was want to integrate aggfunc from pd's pivot_table
