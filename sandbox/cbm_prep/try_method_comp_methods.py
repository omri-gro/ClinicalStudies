from objects import MethodComparator
from sandbox import MetadataBundle
import sandbox as sb
import os
from itertools import product
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Union, Iterable, Optional, Sequence
from pipelines import clv_pipe, medium_pipe


"""
For developing the following MethodComparator methods:
1. Creating new method comparator with only specific variables and only samples that can be compared (value/grade exist for both methods, numeric if necessary) - might want to integrate into to_comparison_matrix
    a. can be done by filtering NaNs in value/grade (or checking if they are numeric) and counting number of unique methods(/investigators) per [SampleID, Site, Variable] combination, keeping only if >1
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
    save_name = 'clv_cbm_both_inv'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    if use_both_inv:
        test_arm = 'CBM'
        ref_arm = 'ClV'
        sites = ['BWH', 'TASMC']
        id_vars_clv = ["SampleID", "Site", "Method", "FileName", 'Investigator']
        id_vars_cbm = ["SampleID", "Site", "Method", "FileName"]
        df_srcs_list = []
        for site in sites:   # this could be used as framework for new pipeline
            df = clv_pipe(f'{site}_ClV_all_revs.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
            df_srcs_list.append(df)

            df = medium_pipe(f'{site}_CBM.csv', site, 'CBM', metadata, dir=r'raw/cbm_method_comparison', id_vars=id_vars_cbm)
            df["Investigator"] = "CBM"
            df_srcs_list.append(df)

        all_dfs = pd.concat(df_srcs_list)
        methd_comp = MethodComparator(all_dfs)
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


    need_arb_path = r'C:\Users\omrig\Downloads\exclude_from_ClV.xlsx'
    methd_comp = methd_comp.filter_by_df(need_arb_path)
    vars_to_test = metadata.variable_groups['RBC morphology'] + metadata.variable_groups['PLT morphology'] + metadata.variable_groups['RBC combinations']



    methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}.csv',
                                        row_identifiers=["SampleID", "Site"],
                                        needed_vals=vars_to_test,
                                        needed_grades=["scan_id"])



    methd_comp.batch_fit(ref_arm, test_arm, vars_to_test)

    methd_comp.save_results(rf'results/{save_name}_all_scans_reg.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_all_scans_reg.pdf')



    print(methd_comp)
