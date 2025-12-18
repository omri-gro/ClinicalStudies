import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    sites = ['CPG', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'

    bin_params = False

    exprt_mtrx = True
    plot_reg = True

    remove_cases_by_list = True
    min_unclass = 10  # number (0-100) or False   currently doesn't matter, as no such slides sent to manual
    min_wbc = 40  # number or False   currently doesn't matter, as no such slides sent to manual
    crf_ssn = 'all'  # 'all', 'pre' or 'post'
    diff500 = True

    save_name = f'mnl_{crf_ssn}-ssn_minuncls_{min_unclass}_minwbc_{min_wbc}_500diff-{diff500}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    invstigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1', 'Sarah Pereira Rodrigues': 'Rev1',
                        'Sladana': 'Rev2', 'Deborah Swearingen': 'Rev2', 'Tony Omigie': 'Rev2',
                        'YAEL ASYEGH': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Yael S': 'Rev2', 'Yael Sayegh': 'Rev2',
                        'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}

    df_srcs_list = []
    for site in sites:
        df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
        df_srcs_list.append(df)
    cbm_file_name = '5sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'
    df_srcs_list.append(cbm_df)

    all_dfs = pd.concat(df_srcs_list)
    all_dfs['Investigator'] = all_dfs['Investigator'].map(invstigators_map)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                                inplace=False)

    if min_unclass:
        methd_comp.df['Value'] = pd.to_numeric(methd_comp.df['Value'], errors='coerce')
        methd_comp = methd_comp.filter_by_unclass(threshold=min_unclass)

    if min_wbc:
        methd_comp = methd_comp.only_when_cond(f"Variable == 'TotalWBC' and Value >= {min_wbc}")

    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = ['scan_id'] + metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    print_also = ['Unclassified WBC', "TotalWBC"]
    vals_to_print = vars_to_test + print_also
    morph_vals_to_print = morph_vals_to_test + print_also

    rmv_file = 'flt_lists/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    if diff500:
        diff500_file = 'flt_lists/500_WBC_mnl_cases.csv'
        diff500_df = read_to_df(diff500_file, file_dir=os.getcwd())
        diff500_comp_df = methd_comp.filter_by_df(diff500_df, include_rows=True).df.copy()

        cbm_df = methd_comp.only_when_cond("Method=='CBM'").df.copy()
        diff500_dfs = pd.concat([diff500_comp_df, cbm_df])
        methd_comp = MethodComparator(diff500_dfs)

    if crf_ssn == "post":
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_means.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method"),
            needed_vals=vals_to_print,
            needed_grades=['scan_id'])

        methd_comp_all_inv.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_all_inv.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=vals_to_print,
            needed_grades=['scan_id'])

        if min_wbc is False:
            methd_comp_all_inv.export_comparison_matrix(
                out_path=fr'comp_tables/{save_name}_grades_all_inv.csv',
                row_identifiers=["Site", "SampleID"],
                comparison_dims=("Variable", "Method", "Investigator"),
                needed_vals=morph_vals_to_test,
                needed_grades=grades_to_test)


    # binary parameters' sensitivity/specificity
    if bin_params:
        binary_vars = metadata.variable_groups["binary"]
        methd_comp.batch_compare('manual', 'CBM', binary_vars, site_filters=sites, function='binary')
        methd_comp.save_results(rf'results/mnl/{save_name}_bin.csv', result_type="binary")

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg.pdf')
