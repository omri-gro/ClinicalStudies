from objects import MethodComparator
from sandbox import MetadataBundle
import sandbox as sb
import os
from itertools import product
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Union, Iterable, Optional, Sequence


"""
For developing the following MethodComparator methods:
1. Creating new method comparator with only specific variables and only samples that can be compared (value/grade exist for both methods, numeric if necessary) - might want to integrate into to_comparison_matrix
2. Positive cases count (using the grades/normal ranges) - saved as another attribute (need to decide if print right there or as part of regression results)
3. Flagging? (see if better here or in wide dataframe)
4. Exclusion of samples according to list
5. Calculation of mean investigator (possibly keeping only the mean reviewer, and possibly only samples with N reviewers)
6. Way to deal with arbitration (keep only arbitrator samples if exist, otherwise keep both regular reviewers)
7. Create contingency tables for grades (calculate metrics and export as part of combined method?)
8. Sensitivity/specificity (with CI filtering added as different method later)
9. Use contingency tables or sensitivity-specificity to export table with variables&sites as rows and TP/FP/TN/FN as columns

Plotting:
1. Inter-user bars to datapoints
"""

def add_pos_column(df_long: pd.DataFrame, meta: "MetadataBundle",
                   normal_grades=[0, "0", "Normal", "Negative", "normal", "negative"]):
    """
    Adds boolean column for positivity based on normal ranges in MetaDataBundle.
    If grade already exists, treat values in normal_vals as False (negative) and rest as True,
    Else use normal ranges and the "Value" column if possible.
    """
    # create empty boolean column
    df = df_long.copy()
    df["Positive"] = np.nan
    df["Positive"] = df["Positive"].astype('boolean')

    # where grade exists, use it for positivity  -  this section might need changing if grade ever not related to positivity
    grade_negative = df["Grade"].isin(normal_grades)
    grade_positive = df["Grade"].notna() & ~df["Grade"].isin(normal_grades)
    df.loc[grade_negative, "Positive"] = False
    df.loc[grade_positive, "Positive"] = True

    # find values where positivity will be based on normal ranges
    value_num = pd.to_numeric(df["Value"], errors="coerce")
    numeric = value_num.notna()
    normal_ranges = getattr(meta, "normal_ranges", {})
    norm_range_vars = set(normal_ranges.keys())
    need_convert = df["Positive"].isna() & numeric & df["Variable"].isin(norm_range_vars)

    # convert based on normal ranges
    if need_convert.any():
        for var in sorted(norm_range_vars):
            mv = need_convert & df["Variable"].eq(var)
            if not mv.any():
                continue
            norm_range = normal_ranges[var]
            if len(norm_range) == 2:
                df.loc[mv, "Positive"] = ~df.loc[mv, "Value"].between(norm_range[0], norm_range[1], inclusive="both")
            else:
                print(f'{norm_range} is not an appropriate normal range for {var}')
    return df










if __name__ == "__main__":
    use_omr = False
    meta_path = r'config.yaml'
    # save_name = 'all_sites_combined_omr_lym'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    if use_omr:
        # build list of files to read from
        sites = ['BWH', 'CPG', 'LMU', 'SYN', 'TASMC']
        mthds = ['OMR', 'CBM']
        srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    else:
        srcs = {('TASMC', 'CBM'): f'TASMC_CBM.csv',
                ('TASMC', 'ClV'): f'TASMC_ClV%.csv',
                ('BWH', 'CBM'): f'BWH_CBM.csv',
                ('BWH', 'ClV'): f'BWH_ClV%_means_decimal.csv'}
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')

    df = methd_comp.df
    # df = add_pos_column(df, metadata)

    need_arb_path = r'C:\Users\omrig\Downloads\exclude_from_ClV.xlsx'
    # no_arb_meth_comp = methd_comp.filter_by_df(need_arb_path)
    # arb_meth_comp = methd_comp.filter_by_df(need_arb_path, include_rows=True)





    print(methd_comp)
