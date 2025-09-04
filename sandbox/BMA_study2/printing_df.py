""" Converting pivoted dataframe to something presentable. an later be implemented as general functions."""
import pandas as pd
import os
import sys
from itertools import product
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle

if __name__ == "__main__":
    meta_path = r'config_BMA.yaml'
    site = 'OHSU'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)
    srcs = {(site, 'REF'): f'{site}_CRF_REF.csv',
            (site, 'TEST'): f'{site}_CRF_TEST.csv'}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=cur_dir, bma=True)



    # start cleaning for representation
    # will be arguments in future function
    needed_vars = metadata.variable_groups['percent']  # + metadata.variable_groups['Evaluation parameter']
    row_identifier = ['SampleID', 'Variable']  # would include 'Site' in many future use cases
    comp_cols = ['Investigator', 'Method']
    values = 'Value'
    decimals = 2  # default to be None, only set this if all values are numeric
    col_ord = [list(methd_comp.df[col].unique()) for col in comp_cols]
    # col_ord = [list(methd_comp.df['Method']),  list(methd_comp.df['Investigator'])]  # order of columns appearance in output - can be list, list of tuples (if multiindex) or list of lists with size of first list same as comp_cols

    needed_rows = methd_comp.df[methd_comp.df['Variable'].isin(needed_vars)]

    # take only SampleID-Variable combinations with right number of rows/filled CRFs - might not always be needed in future function
    m = needed_rows.groupby(row_identifier)['Value'].transform('size').eq(4)
    needed_rows = needed_rows.loc[m]

    out_wide = needed_rows.pivot(index=row_identifier,
                                 columns=comp_cols,
                                 values=values)

    # to do: order of columns - depending on format of input to col_ord
    if col_ord:
        new_col_index = pd.MultiIndex.from_product(col_ord)
        out_wide = out_wide.reindex(columns=new_col_index)

    if isinstance(decimals, int):
        try:  # might want to use different check that number is integer (in case of full number float)
            out_wide = out_wide.astype(float).round(decimals=decimals)
        except ValueError:
            print('Rounding not performed.')



    out_wide.to_csv(r'results/OHSU_wide_57_samples_NDC.csv')


