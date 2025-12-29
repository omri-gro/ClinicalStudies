# utils.py

import yaml
from pathlib import Path
from typing import List, Iterable, Mapping, Optional, Sequence, Union, Tuple
import pandas as pd
import os


""" Validation / Data handling """
def _ensure_list(x):
    """If argument is a tuple, convert to list.
    If argument is anything else but a list, return list containing that one object."""
    if isinstance(x, list):
        return x
    elif isinstance(x, tuple):
        return list(x)
    elif isinstance(x, pd.Series):
        return x.to_list()
    else:
        return [x]



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




""" I/O and config """
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

def read_to_df(file_name, sheet_name='Sheet1', file_dir=None):
    if file_dir == None:
        # dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        file_dir = os.path.abspath(os.path.dirname(__file__))
        file_dir = os.path.join(file_dir, r'raw')
    filepath = os.path.join(file_dir, file_name)
    if not os.path.exists(filepath):
        filepath = file_name
    elif not os.path.exists(filepath):
        raise FileNotFoundError(f"\033[91mFailed to find {file_name}: {e}\033[0m")
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filepath, sheet_name=sheet_name)
            print(f"Loaded Excel file: {filepath} (sheet={sheet_name})")
        elif ext == '.csv':
            df = pd.read_csv(filepath)
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
