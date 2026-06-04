import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools import careful_map, safe_pivot, robust_dup, apply_arbitration_override

if __name__ == "__main__":
    suffix = ''
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    test_arm = 'CBM'
    ref_arm = 'manual'

    bin_params = False
    inter = False

    exprt_long = True
    exprt_mtrx = True
    plot_reg = False

    min_inv = 2  # False or number  currently does not seem to make much of a difference
    rmv_brd = False
    max_unclass = False  # number (0-100) or False   currently doesn't matter, makes not difference to any parameter
    min_wbc = False  # number or False   don't use use value>=100, currently TASMC raw data is percentages
    diff500 = False
    crf_ssn = 'all'  # 'all' or 'post'
    aftr_2nd_ssn = False
    cbm_thresholding = False  # if True than CBM<1% changes to 0 for hairy cells, LGL and atypical
    after_last_ssn = True  # for lym types, Pelger, Auer Rods & RBC distributions, use only after mid-Jan session samples

    max_unclassstr = f'_maxuncls{max_unclass}' if max_unclass else ''
    min_wbcstr = f'_minwbc{min_wbc}' if min_wbc else ''
    rmv_brdstr = '_bdrrmv' if rmv_brd else ''
    diff500str = '_diff500' if diff500 else ''
    scnd_ssn_str = '_aftr2ndssn' if aftr_2nd_ssn else ''
    cbm_thres_str = '_cbm_thres' if cbm_thresholding else ''
    after_last_ssn_str = '_aftrlstssn' if after_last_ssn else ''
    save_name = f'mnl_{crf_ssn}ssn_{max_unclassstr}{min_wbcstr}{rmv_brdstr}_mininv{min_inv}{diff500str}{scnd_ssn_str}{cbm_thres_str}{after_last_ssn_str}{suffix}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    investigators_map = {
        # Standardize typos/variations
        'Christopher Wright': 'Chris',
        'Christine Lavoie': 'Christine',
        'Ebikebuna Rufus': 'Ebi',
        'Ebikebuna Rufus F.': 'Ebi',
        'Ebikebuna Rufus F': 'Ebi',
        'Thu Tran': 'Thu',
        'THU TRAN': 'Thu',
        'Aubrey B Charlton': 'Aubrey',
        'Deborah Swearingen': 'Deborah',
        'Maria Buen Viana De Perio': 'Buen',
        'Joy Arthur': 'Joy',
        'Madison Brooks': 'Madison',
        'Michelle Huynh': 'Michelle',
        'Michelle huynh': 'Michelle',
        'Tiffany I Highsmith': 'Tiffany',
        'Tiffany I. Highsmith': 'Tiffany',
        'Alina KÃƒÂ¼pper': 'Alina',
        'Alina': 'Alina',
        'Sladana Nikolic': 'Sladana',
        'Nikolic Sladana': 'Sladana',
        'Sladana': 'Sladana',
        'Ana Catarina Silva': 'Ana',
        'Harsha Hirani': 'Harsha',
        'Harsha HIrani': 'Harsha',
        'Thomas Muddiman': 'Thomas',
        'Tony Omigie': 'Tony',
        'Sarah Pereira Rodrigues': 'Sarah',
        'YAEL SAYEGH': 'Yael',
        'Yael Sayegh': 'Yael',
        'Yael S': 'Yael',
        'YAEL ASYEGH': 'Yael',

        # Explicitly tag Arbitrators
        'Jared Block': 'Arbitrator',
        'Jennifer Egan': 'Arbitrator',
        'Dr. med. Weigand, Michael': 'Arbitrator',
        'Dr Guy Hannah': 'Arbitrator',
        'Dan BENISTY': 'Arbitrator',

        # Preserve system/automated roles
        'CBM': 'CBM',
        'Mean Investigator': 'Mean Investigator'}

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
        # ONLY include these samples for Aberrant Lymphocyte and Plasma Cell.
        # Other variables (like Segmented Neutrophil) keep all samples.
        methd_comp = methd_comp.filter_by_df(
            diff500_file,
            include_rows=True,
            target_vars=['Aberrant Lymphocyte', 'Plasma Cell']
        )

    if crf_ssn == "post":
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if aftr_2nd_ssn:
        rmv_file = 'flt_lists/pre_2nd_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if after_last_ssn:
        last_ssn_file = 'flt_lists/after_last_session.xlsx'
        methd_comp = methd_comp.filter_by_df(
            last_ssn_file,
            include_rows=True,
            target_vars=['Aberrant Lymphocyte', 'Atypical Lymphocyte', 'LGL', 'Lymphocyte', 'Smudge Cell',
                         'RBC Agglutination', 'Rouleaux', 'Pelger Cell', 'Auer Rods',
                         'Aber&Atyp', 'Variant Lymphocyte']
        )


    df = methd_comp.df
    df = df.query("Value!='--------'").copy()
    # Standardize names and tag the Arbitrator
    df['Investigator'] = careful_map(df['Investigator'], investigators_map)
    # 2. Isolate the Arbitrator
    arb_mask = df['Investigator'] == 'Arbitrator'
    arb_df = df[arb_mask].copy()
    regular_df = df[~arb_mask].copy()

    # Dynamically assign roles (Rev1, Rev2, etc.) ONLY to regular reviewers
    regular_df = assign_dynamic_roles(regular_df, group_cols=['Site', 'SampleID'])

    # Calculate Mean Investigator ONLY on regular reviewers, returning the arbitrator lines afterwards
    regular_df = add_mean_investigator(regular_df, mthd=ref_arm, min_inv=min_inv)
    df = pd.concat([regular_df, arb_df], ignore_index=True)

    # Load arbitration rules (which samples/variables to override) and apply Quantitative Override (Overwrites 'Mean Investigator')
    arb_rules = read_to_df('flt_lists/for_arbitration.csv', file_dir=os.getcwd())
    df = apply_arbitration_override(df, arb_df, arb_rules, metadata)

    # still need to deal with arbitration of binary and graded parameters
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


    cbm_file_name = 'all6_RGB_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, test_arm, metadata, dir=r'raw/cbm_method_comparison')
    # cbm_file_name = 'BWH_newRGB_CBM.csv'
    # cbm_df = medium_pipe(cbm_file_name, 'BWH', test_arm, metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'

    if cbm_thresholding:
        # notice that values below threshold changed to 0% after calculation of percentages, so other percentages will not add up to 100%
        vars_to_thres = ['Hairy Cell', 'Atypical Lymphocyte', 'LGL']
        thres = 1
        cbm_df.loc[(cbm_df['Variable'].isin(vars_to_thres)) & (cbm_df['Value'] < thres), 'Value'] = 0

    all_dfs = pd.concat([df, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    # cases to always remove - horrible slides, horrible scans, etc.
    rmv_file = 'flt_lists/low_quality.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)


    # cases of borderline quality - dirty, investigators' comments on quality, etc.
    if rmv_brd:
        rmv_file = 'flt_lists/low_borderline_quality.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology'] + metadata.variable_groups['RBC arrangement']
    grades_to_print = grades_to_test + ['ScanID']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also
    morph_vals_to_print = morph_vals_to_test + print_also

    if diff500:
        vars_to_test = ['Aberrant Lymphocyte', 'Plasma Cell']
    elif cbm_thresholding:
        vars_to_test = vars_to_thres
    # elif aftr_2nd_ssn:
    #     vars_to_test = ['Aberrant Lymphocyte', 'Atypical Lymphocyte', 'LGL', 'Lymphocyte', 'Smudge Cell']

    if exprt_long:
        include_in_export = vals_to_print + grades_to_print
        df_long = methd_comp.df.query(f"Variable in @include_in_export and Investigator!='Mean Investigator'")[['SampleID', 'Site', 'Method', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long_all_revs.csv')
        df_long = methd_comp.df.query(f"Variable in @include_in_export and Investigator=='Mean Investigator'")[['SampleID', 'Site', 'Method', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long_final_values.csv')


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
