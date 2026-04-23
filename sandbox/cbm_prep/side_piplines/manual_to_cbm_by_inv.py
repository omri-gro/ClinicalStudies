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

    bin_params = False

    exprt_long = False
    exprt_mtrx = False
    plot_reg = True

    min_inv = 2  # False or number
    max_unclass = False  # number (0-100) or False   currently doesn't matter, makes not difference to any parameter
    min_wbc = False  # number or False   don't use use value>=100, currently TASMC raw data is percentages
    diff500 = False
    aftr_2nd_ssn = True

    suffix = '_byInv'
    diff500str = '_diff500' if diff500 else ''
    scnd_ssn_str = '_aftr2ndssn' if aftr_2nd_ssn else ''
    save_name = f'mnl_maxuncls{max_unclass}_minwbc{min_wbc}_mininv{min_inv}_{diff500str}{scnd_ssn_str}{suffix}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)


    investigators_map = {'Alina': 'Alina', 'Aubrey B Charlton': 'Aubrey', 'Thomas Muddiman': 'Thomas',
                         'Sarah Pereira Rodrigues': 'Sarah', 'Maria Buen Viana De Perio': 'Buen',
                         'Christine Lavoie': 'Christine', 'Ebikebuna Rufus': 'Ebi', 'Donald': 'Donald',
                         'Sladana': 'Sladana', 'Deborah Swearingen': 'Deborah', 'Tony Omigie': 'Tony',
                         'Joy Arthur': 'Joy', 'Tiffany I Highsmith': 'Tiffany', 'Tiffany I. Highsmith': 'Tiffany',
                         'Harsha Hirani': 'Harsha',
                         'YAEL ASYEGH': 'Yael', 'YAEL SAYEGH': 'Yael', 'Yael S': 'Yael', 'Yael Sayegh': 'Yael',
                         'Yael S ': 'Yael',
                         'Christopher Wright': 'Chris', 'Thu Tran': 'Thu',
                         'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
    investigators = ['Christine', 'Chris', 'Ebi', 'Thu',
                     'Joy', 'Buen', 'Tiffany',
                     'Aubrey', 'Deborah',
                     'Alina', 'Sladana',
                     'Thomas', 'Tony', 'Donald', 'Harsha',
                     'Sarah', 'Yael']

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

    if aftr_2nd_ssn:
        rmv_file = 'flt_lists/pre_2nd_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    df = methd_comp.df
    df = df.query("Value!='--------'")
    df['Investigator'] = df['Investigator'].map(investigators_map)


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
    rmv_file = 'flt_lists/low_quality.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    rmv_file = 'flt_lists/for_arbitration.csv'
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

    if exprt_long:
        include_in_export = vals_to_print + grades_to_print
        df_long = methd_comp.df.query(f"Variable in @include_in_export and Method=='{ref_arm}' and Investigator!='Mean Investigator'")[['SampleID', 'Site', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long.csv')


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


    for var in vars_to_test:
        for inv in investigators:
            methd_comp.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=var, row_filters={'Investigator': [inv, test_arm]}, split_by='Site')
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg.pdf')
