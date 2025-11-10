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



if __name__ == "__main__":
    use_both_inv = True
    use_bma = False
    use_omr = True
    meta_path = r'config.yaml'
    save_name = 'omr_cbm_sandbox'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    if use_both_inv:
        clv_raw_df = sb.raw_to_df('BWH_ClV_all_revs.csv', 'BWH', 'ClV', dir=r'raw/cbm_method_comparison')
        df = sb.stnd_names(clv_raw_df, metadata.alias_map)
        print(df)

    elif use_bma:
        pass
    elif use_omr:
        # build list of files to read from
        sites = ['BWH', 'CPG', 'LMU', 'SYN', 'TASMC']
        test_arm = 'CBM'
        ref_arm = 'OMR'
        mthds = ['OMR', 'CBM']
        srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    else:
        sites = ['BWH', 'TASMC']
        srcs = {('TASMC', 'CBM'): f'TASMC_CBM.csv',
                ('TASMC', 'ClV'): f'TASMC_ClV%.csv',
                ('BWH', 'CBM'): f'BWH_CBM.csv',
                ('BWH', 'ClV'): f'BWH_ClV%_means_decimal.csv'}
        test_arm = 'CBM'
        ref_arm = 'ClV'
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')

    df = methd_comp.df

    need_arb_path = r'C:\Users\omrig\Downloads\exclude_from_ClV.xlsx'
    # no_arb_meth_comp = methd_comp.filter_by_df(need_arb_path)
    # arb_meth_comp = methd_comp.filter_by_df(need_arb_path, include_rows=True)

    # regression per site and regression overall, plotted
    vars_to_test = 'Monocyte'  # can be string or list
    methd_comp.batch_fit(ref_arm, test_arm, vars_to_test)
    methd_comp.batch_fit(ref_arm, test_arm, vars_to_test, site_filters=sites)
    methd_comp.plot_all_regressions(f'results/sandbox_results/{save_name}_reg.pdf')







    print(methd_comp)
