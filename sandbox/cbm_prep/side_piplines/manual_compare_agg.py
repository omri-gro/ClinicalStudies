import os
import sys
import pandas as pd
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import mean_manual_pipe, medium_pipe
from itertools import product


if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = 'config_agg.yaml'
    save_name = 'manual_LMU_CPG'
    sites = ['LMU', 'CPG']

    metadata = MetadataBundle(meta_path)

    sites_investigators = {'LMU': ['Alina', 'Sladana'], 'CPG': ['Aubrey B Charlton', 'Deborah Swearingen']}

    df_srcs_list = []
    for site in sites:
        investigators = sites_investigators[site]

        # import and pre-process each method separately
        df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
        df_srcs_list.append(df)
        df = medium_pipe(f'{site}_CBM.csv', site, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
        df["Investigator"] = "CBM"
        df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'", inplace=False)


    # aggregated only - mean investigator
    vars_to_test = metadata.variable_groups['derived']
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/{save_name}_agg_reg.csv')

    # separate classes - mean investigator
    methd_comp.clean_calculations()
    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/{save_name}_classes_reg.csv')
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.plot_all_regressions(f'results/{save_name}_classes_reg.pdf')






