import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *



if __name__ == "__main__":
    max_unclass = 1000
    max_smudge = 100000

    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'

    manual = True  # else OMR
    flag_type = 'or'  # 'or' | 'and'

    exprt_mtrx = True
    plot_reg = True
    bin_params = True

    remove_cases_by_list = True

    ref_arm = 'manual' if manual else 'OMR'
    save_name = f'{ref_arm}_cbm_maxuncls{max_unclass}_{flag_type}_maxsmdg{max_smudge}_post'
    rslts_dir = r'results/side_games'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    cbm_file_name = '6sites_CBM.csv'
    df = raw_to_df(cbm_file_name, None, 'CBM', dir=r'raw/cbm_method_comparison')
    df = stnd_names(df, metadata.alias_map)
    df = df.query(f"`Unclassified WBC`<={max_unclass} {flag_type} `Smudge Cell`<={max_smudge}")
    df = check_diff_sum(df, metadata, tolerance=5)
    df = pivot_long(df)
    df = create_derived_variables_long(df, metadata)  # calculate derived variables
    # prepare graded and boolean values
    grd_params = metadata.variable_groups.get("binary") + metadata.variable_groups.get("grade")
    df = add_grade_column(df, metadata)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    cbm_df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID

    if manual:
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        filtering_source = read_to_df(rmv_file, file_dir=os.getcwd())
    else:
        filtering_source = None

    # gather omr or manual as usual
    srcs = {(site, ref_arm): f'{site}_{ref_arm}.csv' for site in sites}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison', inv=manual, filtering_source=filtering_source)
    ref_df = methd_comp.df

    # combine
    df = pd.concat([cbm_df, ref_df])
    methd_comp = MethodComparator(df)

    if remove_cases_by_list:
        rmv_file = 'flt_lists/slides_to_remove.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)


    vars_to_test = metadata.variable_groups['WBC diff']
    grades_to_test = metadata.variable_groups['grade']
    bin_to_test = metadata.variable_groups['binary'] + grades_to_test

    """
    grades_to_print = grades_to_test + ['ScanID']
    bin_to_print = bin_to_test + ['ScanID']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also
    """

    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology']
    grades_to_print = grades_to_test + ['ScanID']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also
    morph_vals_to_print = morph_vals_to_test + print_also

    # binary parameters' sensitivity/specificity
    if bin_params:
        if exprt_mtrx:
            methd_comp.export_comparison_matrix(
                out_path=fr'comp_tables/{save_name}_bin.csv',
                row_identifiers=["Site", "SampleID"],
                comparison_dims=("Variable", "Method"),
                needed_vals=bin_to_test + grades_to_test,
                needed_grades=bin_to_test + grades_to_print,
                value_col="Positive")

        methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=bin_to_test, comp_func='binary')
        methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=bin_to_test, split_by='Site', comp_func='binary')
        methd_comp.save_results(rf'{rslts_dir}/{save_name}_bin.csv', result_type="binary")

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method"),
            needed_vals=vals_to_print,
            needed_grades=['ScanID'])

    methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=vars_to_test)
    methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=vars_to_test, split_by='Site')
    methd_comp.save_results(rf'{rslts_dir}/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'{rslts_dir}/{save_name}_reg.pdf')


