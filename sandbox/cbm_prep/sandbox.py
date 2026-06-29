import pandas as pd
import numpy as np
import os
import re
import yaml
import matplotlib.pyplot as plt
from typing import List, Iterable, Mapping, Optional, Sequence, Union, Tuple, Literal
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools.core.metadata import MetadataBundle
from clinstudtools import safe_pivot, robust_dup
# backwards compatibility
from clinstudtools.preprocessing import calc_diff, fill_nans, apply_src_fixes, raw_to_df, stnd_names, check_diff_sum, add_grade_column, add_pos_column
from clinstudtools.utils import load_yaml, write_df_to_file, _as_df, _ensure_list, read_to_df



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


"""Data import and preparation tools"""



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


def assign_dynamic_roles(df, group_cols=['Site', 'SampleID'], inv_col='Investigator', prefix='Rev'):
    """
    Dynamically assigns relative reviewer roles per sample to avoid hardcoded mappings.
    Optionally preserves the original name in a new column for traceability.
    """
    df = df.copy()
    df = df.reset_index(drop=True)

    # Save the original name just in case you need it for debugging/arbitration
    if 'Original_Investigator' not in df.columns:
        df['Original_Investigator'] = df[inv_col]

    # Assign the dynamic role using a list comprehension to avoid NumPy string errors
    df[inv_col] = (
        df.groupby(group_cols)[inv_col]
        .transform(lambda x: [f"{prefix}{i + 1}" for i in pd.factorize(x)[0]])
    )

    return df


def min_inv_filt(df, mthd, min_inv=2, exact=False):
    """
    Keep only (SampleID, Site, Variable) groups with >= min_inv unique investigators.

    Reporting: print only (SampleID, Site) pairs where the fraction of variables kept
    is < report_if_kept_frac_lt (default: majority of variables fail).
    """
    report_if_kept_frac_lt = 0.5
    subset = ['SampleID', 'Site', 'Variable']
    inv_subset = subset + ['Investigator']

    invs_df = df.query(f"Method=='{mthd}'").copy()

    # error if same investigator reviewed same sample twice
    invs_df = robust_dup(invs_df, key_cols=inv_subset, on_duplicates='raise')

    # check number of investigators that reviewed each sample
    inv_counts = invs_df.groupby(subset)["Investigator"].transform('nunique')
    if exact:
        mask_keep = (inv_counts == min_inv)
        sign = '=='
    else:
        mask_keep = inv_counts >= min_inv
        sign = '>='
    out_df = invs_df.loc[mask_keep].copy()

    # ---- Reporting at (SampleID, Site) level ----
    # total variables present for each (SampleID, Site) in this method/value-filtered frame
    total_vars = (
        invs_df.drop_duplicates(subset=subset)
        .groupby(["SampleID", "Site"])["Variable"]
        .count()
        .rename("n_total_vars")
    )

    # variables that survive (i.e., had >= min_inv investigators)
    kept_vars = (
        out_df.drop_duplicates(subset=subset)
        .groupby(["SampleID", "Site"])["Variable"]
        .count()
        .rename("n_kept_vars")
    )

    summary = total_vars.to_frame().join(kept_vars, how="left").fillna({"n_kept_vars": 0})
    summary["n_kept_vars"] = summary["n_kept_vars"].astype(int)
    summary["kept_frac"] = summary["n_kept_vars"] / summary["n_total_vars"]

    flagged = summary[summary["kept_frac"] < report_if_kept_frac_lt].sort_values(
        ["kept_frac", "n_total_vars", "n_kept_vars"]
    )

    if len(flagged):
        print(
            f"Flagging {len(flagged)} (SampleID, Site) pairs where kept variables fraction "
            f"< {report_if_kept_frac_lt:.2f} after requiring {sign} {min_inv} investigators per variable.\n"
            f"Examples:\n{flagged.head(30).reset_index().to_string(index=False)}"
        )


    return out_df


def add_mean_investigator(df, mthd='ClV', min_inv=0, mean_inv_name="Mean Investigator"):
    # for the observations in mthd, calculate the average of all investigator's reviews
    # observations with less than min_inv investigators' reviews will not have mean calculated
    # returns modified df

    val = "Value"
    subset = ['SampleID', 'Site', 'Variable']
    inv_subset = subset + ['Investigator']

    # take only rows to be used for mean calculation
    invs_df = df.query(f"Method=='{mthd}' and {val}.notna()").copy()

    # only rows with numeric values
    invs_df.loc[:, val] = pd.to_numeric(invs_df.loc[:, val], errors='coerce')
    invs_df = invs_df[invs_df[val].notnull()]

    # error if same investigator reviewed same sample twice
    invs_df = robust_dup(invs_df, key_cols=inv_subset, on_duplicates='raise')

    if min_inv:
        # # check number of investigators that reviewed each sample
        # inv_counts = invs_df.groupby(subset)["SampleID"].transform('count')
        # mask_keep = inv_counts >= min_inv
        #
        # # report on dropped samples
        # dropped = invs_df.loc[~mask_keep, ['SampleID', 'Site']].drop_duplicates()
        # dropped_report_str = f"Removed {dropped.shape[0]} samples with <{min_inv} investigators' reviews.\n" \
        #                      f"Examples:\n{dropped.head(5).to_string(index=False)}"
        # print(dropped_report_str)
        #
        # invs_df = invs_df.loc[mask_keep].copy()
        invs_df = min_inv_filt(invs_df, mthd, min_inv=min_inv)

    # Compute mean for each group
    means_df = invs_df.groupby(subset, as_index=False)[val].mean()

    # add to original
    means_df['Investigator'] = mean_inv_name
    means_df['Method'] = mthd
    means_df['ValueOrigin'] = 'Mean'
    out_df = pd.concat([df, means_df])
    return out_df


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


def create_derived_variables_long(df, metadata, *, on_existing: str = "skip"):
    """
    Create derived variables in a long-format DataFrame using metadata.

    Args:
        df (pd.DataFrame): Long-format DataFrame (concatenated from all sites).
        metadata (MetadataBundle): Metadata containing derived variable recipes.
        on_existing (str): What to do if a derived variable already exists for the same
            group_cols + Variable in df.
            Options:
              - "skip": do not add the derived row for those cases
              - "overwrite": delete existing provided rows and replace with derived rows
              - "keep_both": keep both (previous behavior)

    Returns:
        pd.DataFrame: DataFrame with derived variable rows appended.
    """
    if on_existing not in {"skip", "overwrite", "keep_both"}:
        raise ValueError("on_existing must be one of: 'skip', 'overwrite', 'keep_both'")

    derived_rows = []

    # Possible grouping columns to preserve context
    possible_group_cols = ["SampleID", "Site", "Method", "Investigator", "Unit"]

    # Get grouping columns actually present
    group_cols = [col for col in possible_group_cols if col in df.columns]
    if not group_cols:
        raise ValueError(f"None of the {possible_group_cols} exist in df; cannot safely derive.")

    key_cols = group_cols + ["Variable"]
    existing_keys = df[key_cols].drop_duplicates()


    for derived_var, props in metadata.variables.items():
        derived_from = props.get("derived_from")
        if not derived_from:
            continue  # Skip if no derivation defined

        components = derived_from.get("components", [])
        operation = derived_from.get("operation", "sum")

        # Filter to component rows
        comp_df = df[df["Variable"].isin(components)].copy()
        if len(comp_df) == 0:
            # print(f'Dataframe does not include any of {components}. The derived variable {derived_var} will not be calculated.')
            continue
        comp_df["Value"] = pd.to_numeric(comp_df["Value"], errors='raise')

        aggregated = comp_df.groupby(group_cols, dropna=False).agg({"Value": operation})

        # Build derived rows
        aggregated = aggregated.reset_index()
        aggregated["Variable"] = derived_var
        aggregated["ValueOrigin"] = "Derived"

        # Handle collisions with existing provided rows
        if on_existing in {"skip", "overwrite"}:
            # Build keys for the would-be derived rows
            der_keys = aggregated[key_cols].drop_duplicates()

            # Mark derived rows that already exist in df
            # (merge with indicator is an easy "is in" for multiple columns)
            m = der_keys.merge(existing_keys, on=key_cols, how="left", indicator=True)["_merge"].eq("both")

            if on_existing == "skip":
                # Keep only derived rows whose key is NOT already present
                keep_keys = der_keys.loc[~m, key_cols]
                aggregated = aggregated.merge(keep_keys, on=key_cols, how="inner")

            elif on_existing == "overwrite":
                # Remove any existing rows in df that match derived keys
                # Do it once per derived variable (cheap enough; still vectorized)
                to_remove_keys = der_keys.loc[m, key_cols]
                if not to_remove_keys.empty:
                    df = df.merge(to_remove_keys.assign(_drop=1), on=key_cols, how="left")
                    df = df[df["_drop"].isna()].drop(columns="_drop")

        derived_rows.append(aggregated)

    if not derived_rows:
        print("No derived variables created (no recipes in metadata).")
        return df

    derived_df = pd.concat(derived_rows, ignore_index=True)

    # At this point:
    # - keep_both: df and derived_df are untouched
    # - skip: derived_df has been filtered
    # - overwrite: df has been filtered
    df = pd.concat([df, derived_df], ignore_index=True)
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


def _filter_rows_by_completeness(wide: pd.DataFrame, comparison_dims: list[str], mode: str = "all_cells",):
    # comparison_dims here only for future proofing in case we add more modes
    if mode == "any":
        # keep rows with at least one non-NA anywhere
        return wide.dropna(axis=0, how="all")

    if mode == "all_cells":
        # requires all variables to exists (unlikely to be used)
        return wide.dropna(axis=0, how="any")

    if mode == "all_methods":
        if not isinstance(wide.columns, pd.MultiIndex):
            # single comparison dimension → equivalent to all_cells
            return wide.dropna(axis=0, how="any")

        method_level = wide.columns.names.index("Method")

        # for each method, check row-wise if any column is non-NA
        per_method_present = (
            wide
            .groupby(level=method_level, axis=1)
            .apply(lambda df_: df_.notna().any(axis=1))
        )

        # keep rows where ALL methods are present
        keep_mask = per_method_present.all(axis=1)
        return wide.loc[keep_mask]

    raise ValueError(f"Unknown row completeness mode: {mode}")


def to_comparison_matrix(
        obj_or_df: Union["MethodComparator", pd.DataFrame],
        metadata=None,
        # identifiers that define a datapoint row (only the ones present will be used) - remember to move "Investigator" to here if each investigator in own row
        row_identifiers: Sequence[str] = ("Site", "SampleID"),
        # which columns define the side-by-side comparison blocks
        comparison_dims: Sequence[str] = ("Variable", "Method", "Investigator"),
        value_col: str = "Value",  # or "Grade"
        needed_vars: Optional[Iterable[str]] = None,   # if None and metadata provided, pull from metadata.variable_groups['percent']
        row_completeness: Literal[
            "any",  # keep if any comparison value exists
            "all_cells",  # even a single NaN in row is enough to drop it
            "all_methods"  # require at least one value per Method
            "none"   # no completeness filtering, keep all rows
        ] = "all_methods",
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
    comparison_dims = _ensure_list(comparison_dims)
    row_identifiers = _ensure_list(row_identifiers)

    if value_col == 'Positive':
        df['Positive'] = df['Positive'].astype(str)

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
    if row_completeness != 'none':
        wide = _filter_rows_by_completeness(wide, comparison_dims=present_comp_dims, mode=row_completeness,)

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
    # need to figure out why this section makes values from 'Positive' column disappear
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

    # wide = _round_df(wide, decimals)

    # --- optionally flatten columns for friendlier Excel/CSV ---
    if flatten_columns and isinstance(wide.columns, pd.MultiIndex):
        wide.columns = [sep.join(map(str, col)).strip(sep) for col in wide.columns]
        wide.columns = pd.Index(wide.columns, name=sep.join(present_comp_dims))

    # drop empty columns - turn into boolean argument based if needed
    wide = wide.dropna(axis=1, how='all')

    # sort depends row identifiers argument
    wide = wide.sort_values(by=present_identifiers)

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
    # site=None if Site column already exists in raw data
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
    if site is not None:
        df["Site"] = site

    df["Method"] = method
    df["FileName"] = os.path.basename(file_name)

    # Strip leading/trailing spaces from column names
    df.columns = df.columns.str.strip()

    # to do: lowercase column names?
    return df


