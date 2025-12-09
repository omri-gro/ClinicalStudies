import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    sites = ['CPG', 'LMU', 'TASMC']
    save_name = 'mnl'
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'

    exprt_mtrx = False
    remove_cases_by_list = True

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    df_srcs_list = []
    for site in sites:
        df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
        df_srcs_list.append(df)
    cbm_file_name = '5sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'
    df_srcs_list.append(cbm_df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                                inplace=False)


    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = ['scan_id'] + metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_means.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method"),
            needed_vals=vars_to_test,
            needed_grades=['scan_id'])

        methd_comp_all_inv.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_all_inv.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=vars_to_test,
            needed_grades=['scan_id'])

        methd_comp_all_inv.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_grades_all_inv.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=morph_vals_to_test,
            needed_grades=grades_to_test)

    if remove_cases_by_list:
        rmv_file = 'slides_to_remove_long.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/{save_name}_reg.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_reg.pdf')
