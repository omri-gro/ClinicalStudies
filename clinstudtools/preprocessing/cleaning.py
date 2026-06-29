# clinstudtools/preprocessing/cleaning.py

""" Imputation & Source Corrections """

import pandas as pd

from clinstudtools.table_integrity import robust_dup as _robust_dup


# wrapper to table_integrity's robust_dup
def robust_dup(df, *args, **kwargs):
    """Project-aware version of robust_dup."""
    sort_cols = [c for c in ["Variable", "SampleID"] if c in df.columns]
    return _robust_dup(df, *args, sort_cols=sort_cols, **kwargs)


# finds values that do not have a map and print them as message
def careful_map(sris, map_dict):
    unique_values = set(sris.dropna())
    missing_values = unique_values - set(map_dict.keys())
    if missing_values:
        print(f"\033[91mWARNING: The following values are missing from the dictionary and will become NaN: {missing_values}\033[0m")
    return sris.map(map_dict)


# to do - see if metadata can be removed to convert to function into general use + wrapper
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



