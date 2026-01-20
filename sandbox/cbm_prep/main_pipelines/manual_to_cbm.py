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

    bin_params = True

    exprt_mtrx = True
    plot_reg = False

    remove_cases_by_list = True
    min_unclass = 10  # number (0-100) or False   currently doesn't matter, as no such slides sent to manual
    min_wbc = 40  # number or False   currently doesn't matter, as no such slides sent to manual
    crf_ssn = 'all'  # 'all', 'pre' or 'post'
    diff500 = False

    save_name = f'mnl_{crf_ssn}-ssn_minuncls_{min_unclass}_minwbc_{min_wbc}_500diff-{diff500}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    invstigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1',
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
        raw_df = raw_to_df(f'{site}_manual.csv', site, "manual", dir=r'raw/cbm_method_comparison')
        df = stnd_names(raw_df, metadata.alias_map)
        df = calc_diff(df, metadata, additional_cells="WBC-like")
        df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
        df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp = MethodComparator(all_dfs)

    """
    if min_unclass:
        methd_comp.df['Value'] = pd.to_numeric(methd_comp.df['Value'], errors='coerce')
        methd_comp = methd_comp.filter_by_unclass(threshold=min_unclass)

    if min_wbc:
        methd_comp = methd_comp.only_when_cond(f"Variable == 'Total WBC' and Value >= {min_wbc}")
    """


    if diff500:
        diff500_file = 'flt_lists/500_WBC_mnl_cases.csv'
        diff500_df = read_to_df(diff500_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(diff500_df, include_rows=True)

    if crf_ssn == "post":
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    # calculate mean investigator
    df = methd_comp.df

    df = df.query("Value!='--------'")

    df['Investigator'] = df['Investigator'].map(invstigators_map)
    df = add_mean_investigator(df, 'manual')

    binary_vars = metadata.variable_groups["PLT morphology"] + metadata.variable_groups["RBC arrangement"] + metadata.variable_groups["WBC morphology"]
    raw_grade_cond=lambda d: (
        d["Method"].isin(['manual'])
        & d["Variable"].isin(binary_vars)
    )

    df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value or grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)


    cbm_file_name = '6sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'

    all_dfs = pd.concat([df, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    rmv_file = 'flt_lists/slides_to_remove.csv'
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

    # binary parameters' sensitivity/specificity
    if bin_params:
        if exprt_mtrx:
            methd_comp.df['Positive'] = methd_comp.df['Positive'].astype(str)
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

        methd_comp.batch_compare(levels_a='manual', levels_b='CBM', variables=binary_vars, comp_func='binary')
        methd_comp.batch_compare(levels_a='manual', levels_b='CBM', variables=binary_vars, split_by='Site', comp_func='binary')
        methd_comp.save_results(rf'results/mnl/{save_name}_bin.csv', result_type="binary")

    # keep only mean investigator
    methd_comp = methd_comp.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                                inplace=False)

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_means.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method"),
            needed_vals=vals_to_print,
            needed_grades=['ScanID'])



    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg.pdf')
