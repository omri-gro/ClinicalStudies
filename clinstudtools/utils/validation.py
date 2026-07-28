# utils/validation.py

from typing import Union
import pandas as pd

def expect_single(
    items,
    *,
    what: str = "item",
    context: Union[str, None] = None
):
    """
    Assert that `items` contains exactly one element and return it.
    """
    n = len(items)
    if n != 1:
        msg = f"Expected exactly one {what}, found {n}"
        if context:
            msg += f" ({context})"
        raise ValueError(msg)
    return items[0]


def validate_numeric_columns(df, columns):
    """
    Checks specified columns for non-numeric values, prints informative warnings
    directing the user to the exact rows and values, and coerces the columns to numeric.

    Args:
        df (pd.DataFrame): The DataFrame to check.
        columns (list or str): A column name or list of column names to validate.

    Returns:
        pd.DataFrame: The DataFrame with the specified columns coerced to numeric.
    """
    if isinstance(columns, str):
        columns = [columns]

    # Filter to only columns that actually exist to avoid KeyErrors
    columns = [col for col in columns if col in df.columns]

    found_errors = False

    for col in columns:
        # Only check columns that Pandas doesn't already recognize as strictly numeric
        if not pd.api.types.is_numeric_dtype(df[col]):
            numeric_view = pd.to_numeric(df[col], errors="coerce")

            # Find rows that became NaN but weren't NaN originally
            bad_mask = numeric_view.isna() & df[col].notna()

            if bad_mask.any():
                found_errors = True
                bad_indices = df.index[bad_mask].tolist()
                print(f"\033[91m[TYPE ERROR] Column '{col}' contains non-numeric data.\033[0m")
                print(f"  -> Problematic row indices (showing up to 5): {bad_indices[:5]}")
                print(f"  -> Problematic values: {df.loc[bad_mask, col].head().tolist()}")

            # Coerce the column to numeric so future math operations don't crash
            df[col] = numeric_view

    if found_errors:
        print("Coerced non-numeric values to NaN to proceed safely with calculations.")

    return df
