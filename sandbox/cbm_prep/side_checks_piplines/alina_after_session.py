import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from itertools import *
from pipelines import mean_manual_pipe, medium_pipe



if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    save_name = 'manual_after_Alina_session'
    meta_path = r'config.yaml'
    site = 'LMU'
    investigator = 'Alina'

    metadata = MetadataBundle(meta_path)

    srcs = {('LMU', 'CBM'): f'LMU_CBM.csv',
            ('LMU', 'Manual'): f'LMU_Manual.csv'}

    df_srcs_list = []
    df = mean_manual_pipe(f'{site}_{save_name}.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
    df_srcs_list.append(df)
    df = medium_pipe(f'{site}_CBM.csv', site, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    df["Investigator"] = "CBM"
    df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', f"Investigator=='{investigator}' or Investigator=='CBM'", inplace=False)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']

    methd_comp.export_comparison_matrix(
        out_path=fr'comp_tables/{save_name}.csv',
        row_identifiers=["SampleID"],
        comparison_dims=("Variable", "Method"),
        needed_vars=vars_to_test)

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.save_results(rf'results/{save_name}_reg.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_reg.pdf')

