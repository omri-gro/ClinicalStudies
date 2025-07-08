import pandas as pd
import os


# creating derived_df that is a pivoted dataframe with row for each measurement (parameter, sample, site, method combination)
def read_to_df(file_name, sheet_name=None, dir=None):
    if dir == None:
        dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        dir = os.path.join(dir, r'raw')
    filepath = os.path.join(dir, file_name)
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filepath, sheet_name=sheet_name)
            print(f"Loaded Excel file: {filepath} (sheet={sheet_name})")
        elif ext == 'csv':
            df = pd.read_csv(filepath)
            print(f"Loaded CSV file: {filepath}")
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")

    # Basic sanity check
    if df.empty:
        raise ValueError(f"File {filepath} is empty.")

    return df

def read_to_df(file_name, site, method, sheet_name=None, dir=None):
    if dir == None:
        dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        dir = os.path.join(dir, r'raw')
    filepath = os.path.join(dir, file_name)
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filepath)
        elif ext == 'csv':
            df = pd.read_csv(filepath)
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")

    # Basic sanity check
    if df.empty:
        raise ValueError(f"File {filepath} is empty.")

    # Standardize the SampleID
    possible_id_cols = ["SampleID", "Sample", "Sample ID", "ID", "Barcode", "barcode", "Case"]
    id_col = next((col for col in possible_id_cols if col in df.columns), None)
    if not id_col:
        raise ValueError(f"No sample ID column found in {filepath}")
    df.rename(columns={id_col: "SampleID"}, inplace=True)

    if "SampleID" not in df.columns:
        raise ValueError(f"Missing 'SampleID' column in {filepath}")

    # to do sanity check for empty files, looking for the SampleID column
    # to do strip leading/trailing spaces from column names, lowercase column names, print of table shape
    return df



def clean_omr():
    pass






if __name__ == "__main__":