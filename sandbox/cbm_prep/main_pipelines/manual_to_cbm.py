import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    test_arm = 'CBM'
    ref_arm = 'manual'

    bin_params = True
    inter = False

    exprt_mtrx = True
    plot_reg = True

    min_inv = False  # False or number
    rmv_brd = False
    max_unclass = False  # number (0-100) or False   currently doesn't matter, makes not difference to any parameter
    min_wbc = False  # number or False   don't use use value>=100, currently TASMC raw data is percentages
    diff500 = False
    crf_ssn = 'all'  # 'all' or 'post'
    aftr_2nd_ssn = True

    diff500str = '_diff500' if diff500 else ''
    scnd_ssn_str = '_aftr2ndssn' if aftr_2nd_ssn else ''
    save_name = f'mnl_{crf_ssn}ssn_maxuncls{max_unclass}_minwbc{min_wbc}_mininv{min_inv}_brdrmv-{rmv_brd}{diff500str}{scnd_ssn_str}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    investigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1',
                        'Sarah Pereira Rodrigues': 'Rev1', 'Maria Buen Viana De Perio': 'Rev1',
                        'Christine Lavoie': 'Rev1', 'Ebikebuna Rufus': 'Rev1',
                        'Sladana': 'Rev2', 'Deborah Swearingen': 'Rev2', 'Tony Omigie': 'Rev2',
                        'Joy Arthur': 'Rev2', 'Tiffany I Highsmith': 'Rev2', 'Tiffany I. Highsmith': 'Rev2',
                        'YAEL ASYEGH': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Yael S': 'Rev2', 'Yael Sayegh': 'Rev2',
                        'Yael S ': 'Rev2',
                        'Christopher Wright': 'Rev2', 'Thu Tran': 'Rev2',
                        'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}

    df_srcs_list = []
    for site in sites:
        raw_df = raw_to_df(f'{site}_manual.csv', site, ref_arm, dir=r'raw/cbm_method_comparison')
        df = stnd_names(raw_df, metadata.alias_map)
        df = calc_diff(df, metadata, additional_cells="WBC-like")
        df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
        df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp = MethodComparator(all_dfs)


    if max_unclass:
        methd_comp.df['Value'] = pd.to_numeric(methd_comp.df['Value'], errors='coerce')
        methd_comp = methd_comp.filter_by_unclass(threshold=max_unclass)

    if min_wbc:
        methd_comp = methd_comp.only_when_cond(f"Variable == 'Total WBC' and Value >= {min_wbc}")

    if diff500:
        diff500_file = 'flt_lists/500_WBC_mnl_cases.csv'
        diff500_df = read_to_df(diff500_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(diff500_df, include_rows=True)

    if crf_ssn == "post":
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if aftr_2nd_ssn:
        rmv_file = 'flt_lists/pre_2nd_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    # calculate mean investigator
    df = methd_comp.df

    df = df.query("Value!='--------'")

    df['Investigator'] = df['Investigator'].map(investigators_map)
    df = add_mean_investigator(df, mthd=ref_arm, min_inv=min_inv)

    binary_vars = metadata.variable_groups["PLT morphology"] + metadata.variable_groups["RBC arrangement"] + metadata.variable_groups["WBC morphology"]
    raw_grade_cond=lambda d: (
        d["Method"].isin([ref_arm])
        & d["Variable"].isin(binary_vars)
    )

    df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value nor grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)


    cbm_file_name = '6sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, test_arm, metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'

    all_dfs = pd.concat([df, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    # cases to always remove - waiting for arbitration, horrible scans, etc.
    rmv_file = 'flt_lists/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    # cases of borderline quality - dirty, investigators' comments on quality, etc.
    if rmv_brd:
        rmv_file = 'flt_lists/slides_to_remove_borderline.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology']
    grades_to_print = grades_to_test + ['ScanID']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also
    morph_vals_to_print = morph_vals_to_test + print_also


    if inter:
        if bin_params:
            methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=binary_vars,
                                     dim_col='Investigator', comp_func='binary')
            methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=binary_vars,
                                     dim_col='Investigator', split_by='Site', comp_func='binary')
            methd_comp.save_results(rf'results/mnl/{save_name}_bin_inter.csv', result_type="binary")

        methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=vars_to_test,
                                 dim_col='Investigator')
        methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=vars_to_test,
                                 dim_col='Investigator', split_by='Site')
        methd_comp.save_results(rf'results/mnl/{save_name}_reg_inter.csv')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg_inter.pdf')
        methd_comp.clean_calculations()

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_all_inv.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=vals_to_print,
            needed_grades=['ScanID'])

        if min_wbc is False:
            methd_comp.export_comparison_matrix(
                out_path=fr'comp_tables/{save_name}_grades_all_inv.csv',
                row_identifiers=["Site", "SampleID"],
                comparison_dims=("Variable", "Method", "Investigator"),
                needed_vals=morph_vals_to_test,
                needed_grades=grades_to_print)

    if bin_params and exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_bin.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vars=binary_vars,
            value_col="Positive")
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_bin_vals.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=binary_vars,
            needed_grades=binary_vars + ['ScanID'],
            row_completeness="none")

    # keep only mean investigator
    methd_comp = methd_comp.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                        inplace=False)

    # binary parameters' sensitivity/specificity
    if bin_params:
        methd_comp.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=binary_vars, comp_func='binary')
        methd_comp.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=binary_vars, split_by='Site', comp_func='binary')
        methd_comp.save_results(rf'results/mnl/{save_name}_bin.csv', result_type="binary")

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_means.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method"),
            needed_vals=vals_to_print,
            needed_grades=['ScanID'])



    methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test)
    methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg.pdf')
