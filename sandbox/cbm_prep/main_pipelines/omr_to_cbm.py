import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import medium_pipe
from itertools import product


if __name__ == "__main__":
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    cbm_version = 'v325'  # v317 / v319 / v325

    exprt_mtrx = True
    plot_reg = True
    bin_params = False

    remove_cases_by_list = False
    diff500 = False   # only when manual
    only_good_sites = False

    manual = False  # if False use OMR as reference arm

    # regression settings
    reg_mthd = 'deming'   # 'deming' | 'passing'

    if diff500 and manual:
        ref_arm_wbcs = 1000
    elif manual:
        ref_arm_wbcs = 400
    else:
        ref_arm_wbcs = 200

    test_arm_wbcs = 1000
    lambda_ = ref_arm_wbcs / test_arm_wbcs   # only relevant for deming

    ref_arm = 'manual' if manual else 'OMR'
    short_ref_arm = 'mnl' if manual else 'omr'

    if only_good_sites:
        sites = ['BWH', 'HUP', 'LMU', 'TASMC']
    else:
        sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']

    mthd_str = f'{reg_mthd}{lambda_}' if (reg_mthd == 'deming') else f'{reg_mthd}'
    diff500_string = '_500wbc' if diff500 else ''
    good_sites_string = '_only_good_sites' if only_good_sites else ''

    save_name = f'{short_ref_arm}_cbm_{mthd_str}{diff500_string}{good_sites_string}_{cbm_version}'
    rslts_dir = rf'results/{short_ref_arm}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    cbm_file_name = f'all6_RGB_CBM_{cbm_version}.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    # gather omr as usual
    srcs = {(site, ref_arm): f'{site}_{ref_arm}.csv' for site in sites}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison', inv=manual)

    if diff500:
        diff500_file = 'flt_lists/500_WBC_mnl_cases.csv'
        diff500_df = read_to_df(diff500_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(diff500_df, include_rows=True)

    omr_df = methd_comp.df


    # combine
    df = pd.concat([cbm_df, omr_df])
    methd_comp = MethodComparator(df)

    if remove_cases_by_list:
        rmv_file = 'flt_lists/slides_to_remove.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = metadata.variable_groups['grade']
    bin_to_test = metadata.variable_groups['binary'] + grades_to_test

    grades_to_print = grades_to_test + ['ScanID']
    bin_to_print = bin_to_test + ['ScanID']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also

    # binary parameters' sensitivity/specificity
    if bin_params:
        if exprt_mtrx:
            methd_comp.export_comparison_matrix(
                out_path=fr'comp_tables/{save_name}_bin.csv',
                row_identifiers=["Site", "SampleID"],
                comparison_dims=("Variable", "Method"),
                needed_vals=bin_to_test,
                needed_grades=bin_to_print,
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

    methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=vars_to_test,
                             lambda_=lambda_, reg_method=reg_mthd)
    methd_comp.batch_compare(levels_a=ref_arm, levels_b='CBM', variables=vars_to_test, split_by='Site')
    methd_comp.save_results(rf'{rslts_dir}/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'{rslts_dir}/{save_name}_reg.pdf')
