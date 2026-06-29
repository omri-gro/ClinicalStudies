# utils/io.py

import yaml
import os
import pandas as pd
from pathlib import Path
from typing import Union


def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def write_df_to_file(df: pd.DataFrame, out_path: Union[str, Path]):
    """ Write dataframe to either Excel or csv"""
    format = str(out_path).split(".")[-1]
    if format.lower() == "csv":
        df.to_csv(out_path, index=False)
    elif format.lower() in ("xlsx", "excel"):
        df.to_excel(out_path, index=False)
    else:
        raise ValueError("Format must be 'csv' or 'excel'.")


def read_to_df(file_name, sheet_name='Sheet1', file_dir=None, encodings=['utf-8', 'windows-1252', 'latin1'], **kwargs):
    if file_dir is None:
        file_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "raw"))

    filepath = os.path.join(file_dir, file_name)
    if not os.path.exists(filepath):
        filepath = file_name
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"\033[91mFailed to find {file_name}\033[0m")

    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filepath, sheet_name=sheet_name, **kwargs)
        elif ext == '.csv':
            for encoding in encodings:
                try:
                    # Attempt to read the CSV with the current encoding
                    df = pd.read_csv(filepath, encoding=encoding, **kwargs)
                    return df
                except (UnicodeDecodeError, LookupError):
                    # Catching UnicodeDecodeError for bad bytes
                    # and LookupError in case an invalid encoding name is provided
                    continue
            raise ValueError(
                f"Could not decode the file '{filepath}' with any of the attempted encodings: {encodings}"
            )
        else:
            raise ValueError(f"\033[93mUnsupported file format: {ext}\033[0m")
    except Exception as e:
        raise RuntimeError(f"\033[91mFailed to load {filepath}: {e}\033[0m")

    # Basic sanity check
    if df.empty:
        raise ValueError(f"\033[91mFile {filepath} is empty.\033[0m")

    print(f"Loaded {filepath} with shape {df.shape}")
    return df
