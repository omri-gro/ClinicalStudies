import pandas as pd
import numpy as np
import os
import re
import yaml
import matplotlib.pyplot as plt
from typing import List
from matplotlib.backends.backend_pdf import PdfPages


# creating derived_df that is a pivoted dataframe with row for each measurement (parameter, sample, site, method combination)
class MetadataBundle:
    #  Holds all metadata and configuration for the pipeline.
    def __init__(self, metadata_path):
        context = self.load_yaml(metadata_path)
        self.variables = context["variables"]
        self.alias_map = self.build_alias_map(self.variables)
        self.variable_groups = self.build_variable_groups(self.variables)
        self.crit_points = self.build_crit_points(self.variables)

        self.src_fixes = context.get("Source fixes", {})

    def load_yaml(self, path):
        with open(path, "r") as f:
            return yaml.safe_load(f)

    def build_alias_map(self, variables):
        alias_map = {}
        for canonical, props in variables.items():
            alias_map[canonical] = canonical  # map canonical to itself
            for alias in props.get("aliases", []):
                alias_map[alias] = canonical
        return alias_map

    def get_variable_groups(self, group_name):
        return [var for var, props in self.variables.items()
                if group_name in props.get("groups", [])]

    def build_variable_groups(self, variables):
        """
        Precompute all variable groups as a dict: group_name → list of variable names.
        """
        group_map = {}
        for var, props in variables.items():
            for group in props.get("groups", []):
                group_map.setdefault(group, []).append(var)
        return group_map

    def build_crit_points(self, variables):
        """
        Create dictionary of quantitative variables and lists of critical points.
        """
        points_map = {}
        for var, props in variables.items():
            crit_points = props.get("crit_points")
            if isinstance(crit_points, list):
                points_map[var] = crit_points
        return points_map


def read_to_df(file_name, sheet_name='Sheet1', file_dir=None):
    if file_dir == None:
        # dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        file_dir = os.path.abspath(os.path.dirname(__file__))
        file_dir = os.path.join(file_dir, r'raw')
    filepath = os.path.join(file_dir, file_name)
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


def raw_to_df(file_name, site, method, sheet_name='Sheet1', dir=None):
    def standardize_sample_ids(df, id_col="SampleID"):
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
            match = re.search(r"(\d+)$", str(raw_id))
            if match:
                numeric_part = match.group(1).lstrip("0") or "0"  # Preserve 0 if all zeros
                return numeric_part
            else:
                print(f"\033[91Could not parse SampleID: {raw_id}\033[0m")
                return raw_id  # fallback: return raw

        df[id_col] = df[id_col].apply(clean_id)

        # check for sampleIDs duplicates
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
    df = standardize_sample_ids(df, id_col="SampleID")

    # Add metadata columns
    df["Site"] = site
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
        src (str): Name of site or method (which appears in the metadata's source fixes).

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
        else:
            print(f"\033[91Unknown fix rules group: {rules_group}\033[0m")

    return df


def check_diff_sum(df, metadata, diff_cells="WBC diff", tolerance=5, auto_convert=True, drop_flagged=True):
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

    Prints:
        A warning for each sample where the sum of WBC diff variables
        falls outside (100 ± tolerance), unless they are NaN.

    Returns:
        df: possibly modified DataFrame (with scaled WBC diff variables)
    """
    # TO do: if WBC diff needs dividing by 100, do same for percentage-like
    # Get WBC diff variables from context metadata
    wbc_diff_vars = metadata.variable_groups.get(diff_cells, [])
    wbc_diff_vars = [var for var in wbc_diff_vars if var in df.columns]

    if not wbc_diff_vars:
        print(f"\033[91mNo {diff_cells} variables found in DataFrame.\033[0m")
        return

    # Compute row-wise sum of WBC diff variables
    wbc_sum = df[wbc_diff_vars].sum(axis=1, skipna=True)

    # Identify rows where all WBC diff variables are NaN
    all_nan = df[wbc_diff_vars].isna().all(axis=1)

    # Auto-convert if nearly all rows are close to 1.0
    close_to_1 = wbc_sum.between(1 - tolerance / 100, 1 + tolerance / 100)
    n_valid = (~all_nan).sum()
    n_close_to_1 = close_to_1[~all_nan].sum()

    if auto_convert and n_valid > 0 and n_close_to_1 / n_valid > 0.9:
        df[wbc_diff_vars] = df[wbc_diff_vars] * 100
        wbc_sum = df[wbc_diff_vars].sum(axis=1, skipna=True)
        print(f"[INFO] {diff_cells} values appeared to be fractions. Converted to percentages by multiplying by 100.")

    # Determine rows where sum is outside acceptable range
    lower_bound = 100 - tolerance
    upper_bound = 100 + tolerance
    out_of_range = ~wbc_sum.between(lower_bound, upper_bound)

    # Final flagged rows (excluding all-NaN rows)
    flagged = out_of_range & ~all_nan

    # Report
    if flagged.any():
        sample_ids = df.loc[flagged, "SampleID"]
        for sid, total in zip(sample_ids, wbc_sum[flagged]):
            if drop_flagged:
                print(f"\033[91mSample {sid} dropped: {diff_cells} sum = {total:.1f}%\033[0m")
            else:
                print(f"\033[91mSample {sid}: {diff_cells} sum = {total:.1f}% (expected ~100%)\033[0m")
    else:
        print(f"All {diff_cells} sums within acceptable range.")

    if drop_flagged:
        df = df[~flagged]

    return df


def calc_diff(df, metadata, diff_cells="WBC diff", additional_cells=None):
    """
    Convert cell counts to cell differential, with
    Args:
        df (pd.DataFrame): The DataFrame to calculate for.
        metadata (MetadataBundle): Object with variable groups metadata.
        diff_cells (str or list): Group of cells included in the differential.
        additional_cells (str or list, optional): Group of cells to be converted to percentage according to the totalcells in diff_cells.
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
        return df

    # Calculate total WBC count (sum of diff cells)
    df["_total_WBC"] = df[diff_vars].sum(axis=1, skipna=True)

    # Replace diff cells with percentages
    for var in diff_vars:
        df[var] = df[var] / df["_total_WBC"]
        print(f"Updated '{var}' to percentage of total WBC.")

    # Replace additional cells with percentages relative to total WBC
    for var in additional_vars:
        df[var] = df[var] / df["_total_WBC"]
        print(f"Updated '{var}' to percentage of total WBC (additional cell).")

    # Drop helper column
    df.drop(columns="_total_WBC", inplace=True)

    return df


def curate_df(df, metadata, src=None, wbcs_as_counts=False):
    # to do: check the type of results' source (scopio, OMR from specific site, CRF, etc.)
    if not src:
        if df['Method'][0] in ['CBM', 'Scopio', 'Test', 'BMA', 'ClV']:
            src = df['Method'][0]
        else:
            src = df['Site'][0]

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

    return df


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


""" pipeline functions - will change later """
def short_pipe(df, metadata):
    # standardize column names
    df = stnd_names(df, metadata.alias_map)

    # print warning if WBCs in differential don't add up to ~100
    check_diff_sum(df, metadata, tolerance=5)

    df = pivot_long(df)
    df = df.dropna(subset=["Value"])

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    return df
