from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
import sandbox as sb
import os
from itertools import product
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Union, Iterable, Optional, Sequence
from pipelines import clv_pipe, medium_pipe


"""
Taken from try_method_comp_methods.
"""

if __name__ == "__main__":

    remove_cases_by_list = True
    save_mtrx = False
    plot_reg = True
    # invs_names = ['Christine Lavoie', 'Christopher Wright', 'Thu Tran', 'Ebikebuna Rufus',
    #               'YAEL SAYEGH', 'Sarah Pereira Rodrigues']
    # invs_str = ','.join(invs_names)
    save_name = f'clv_cbm_agg_without_scratched'


    meta_path = r'special_config_files/config_agg.yaml'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)
    vars_to_test = metadata.variable_groups['RBC shape'] + metadata.variable_groups['RBC combinations']

    test_arm = 'CBM'
    ref_arm = 'ClV'
    sites = ['BWH', 'TASMC']
    id_vars_clv = ["SampleID", "Site", "Method", "FileName", 'Investigator']
    id_vars_cbm = ["SampleID", "Site", "Method", "FileName"]
    df_srcs_list = []
    for site in sites:   # this could be used as framework for new pipeline
        df = clv_pipe(f'{site}_ClV_all_revs.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=True)
        df_srcs_list.append(df)

    df = medium_pipe(f'5sites_CBM.csv', None, 'CBM', metadata, dir=r'raw/cbm_method_comparison',
                     id_vars=id_vars_cbm)
    df["Investigator"] = "CBM"
    df = df.query("Site in @sites")
    df_srcs_list.append(df)
    all_dfs = pd.concat(df_srcs_list)
    methd_comp = MethodComparator(all_dfs)


    # vars_to_test = metadata.variable_groups['RBC morphology'] + metadata.variable_groups['PLT morphology'] + metadata.variable_groups['RBC combinations']
    # just for games

    if save_mtrx:
        methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}.csv',
                                            row_identifiers=["SampleID", "Site"],
                                            needed_vals=vars_to_test,
                                            needed_grades=["scan_id"])

    if remove_cases_by_list:
        rmv_file = 'slides_to_remove_long.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        rmv_df["SampleID"] = rmv_df["SampleID"].astype(str).str.zfill(5)

        methd_comp = methd_comp.filter_by_df(rmv_df)

    # only wanted investigators
    # save_name = f'clv_agg_just_Ebi'
    # df = methd_comp.df
    # tasmc_rmvd_inv = 'Ebikebuna Rufus'
    # df = df.query(f"Investigator == '{tasmc_rmvd_inv}' or Investigator == 'CBM' and Site == 'BWH'")
    # methd_comp.df = df




    methd_comp.batch_fit(ref_arm, test_arm, vars_to_test)
    methd_comp.batch_fit(ref_arm, test_arm, vars_to_test, site_filters=sites)

    methd_comp.save_results(rf'results/side_games/{save_name}_reg.csv')

    if plot_reg:
        methd_comp.plot_all_regressions(f'results/side_games/{save_name}_reg.pdf')



