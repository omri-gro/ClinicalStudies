import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools import careful_map, apply_arbitration_override
from clinstudtools.transforms import filter_by_reference, filter_by_condition, filter_samples_by_condition, filter_samples_by_multiple_conditions
from clinstudtools.utils import read_to_df


if __name__ == "__main__":
    suffix = ''
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    test_arm = 'CBM'
    ref_arm = 'manual'
    cbm_version = 'v325'  # v317 / v319 / v325

    bin_params = True
    inter = False
    by_site = True

    only_hard_params = False

    exprt_long = True
    exprt_mtrx = True
    plot_reg = True

    min_inv = 2  # False or number  currently does not seem to make much of a difference
    rmv_brd = False
    max_unclass_cbm = False  # number (0-100) or False   currently doesn't matter, makes not difference to any parameter
    min_wbc_mnl = False  # number or False   don't use use value>=100, currently TASMC raw data is percentages
    min_wbc_cbm = False  # number or False
    diff500 = False
    crf_ssn = 'all'  # 'all' or 'post'
    aftr_2nd_ssn = False  # suggest not to use this one
    cbm_thresholding = False  # if True than CBM<1% changes to 0 for hairy cells, LGL and atypical
    after_last_ssn = False  # for lym types, smudge, Pelger, Auer Rods & RBC distributions, use only after mid-Jan session samples
    no_cpg = False
    by_inv = False
    no_arb_cands = False    # currently not in use - for checking arbitration request for additional slides
    rmv_tech_flgs = True
    max_unclass_cbm = False  # number or False
    lym3sites = True  # only ['BWH', 'SYN', 'TASMC'] for the lymphocytes subtypes


    max_unclassstr = f'_maxuncls{max_unclass_cbm}' if max_unclass_cbm else ''
    min_wbcmnlstr = f'_minmnlwbc{min_wbc_mnl}' if min_wbc_mnl else ''
    min_wbccbmstr = f'_mincbmwbc{min_wbc_cbm}' if min_wbc_cbm else ''
    rmv_brdstr = '_bdrrmv' if rmv_brd else ''
    min_inv_str = f'_mininv{min_inv}' if min_inv else ''
    diff500str = '_diff500' if diff500 else ''
    scnd_ssn_str = '_aftr2ndssn' if aftr_2nd_ssn else ''
    cbm_thres_str = '_cbm_thres' if cbm_thresholding else ''
    after_last_ssn_str = '_aftrlstssn' if after_last_ssn else '_all_ssns'
    cbm_version_str = f'_{cbm_version}' if cbm_version else ''
    no_cpg_str = '_no_CPG' if no_cpg else ''
    no_arb_cands_str = '_NoArbCands' if no_arb_cands else ''
    rmv_tech_flgs_str = '_NoTechFlgs' if rmv_tech_flgs else ''
    unclass_str = f'_{max_unclass_cbm}unclass' if max_unclass_cbm else ''
    lym3sites_str = '_lym3sites' if lym3sites else ''
    save_name = f'mnl_{max_unclassstr}{min_wbcmnlstr}{min_wbccbmstr}{rmv_brdstr}{min_inv_str}{diff500str}{scnd_ssn_str}{cbm_thres_str}{after_last_ssn_str}{no_cpg_str}{cbm_version_str}{no_arb_cands_str}{rmv_tech_flgs_str}{unclass_str}{lym3sites_str}{suffix}'


    # string for filtering raw CBM file (more can be added later)
    raw_cbm_cond = f"`Total WBC`>={min_wbc_cbm}" if min_wbc_cbm else None


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
        'thu tran': 'Thu',
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
        'Dr. med. Michael Weigand': 'Arbitrator',
        'Dr Guy Hannah': 'Arbitrator',
        'Andrew Doyle': 'Arbitrator',
        'Ben-Zion Katz': 'Arbitrator',
        'Dan BENISTY': 'Arbitrator',
        'Olga Pozdnyakova': 'Arbitrator',
        'Christopher Hergott': 'Arbitrator',
        'Robert P Hasserjian': 'Arbitrator',

        # Preserve system/automated roles
        'CBM': 'CBM',
        'Mean Investigator': 'Mean Investigator'}

    if no_cpg:
        by_site = False
        sites = ['BWH', 'HUP', 'LMU', 'SYN', 'TASMC']

    df_srcs_list = []
    for site in sites:
        raw_df = raw_to_df(f'{site}_manual.csv', site, ref_arm, dir=r'raw/cbm_method_comparison')
        df = stnd_names(raw_df, metadata.alias_map)
        df = calc_diff(df, metadata, additional_cells="WBC-like")
        df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
        df_srcs_list.append(df)

    df = pd.concat(df_srcs_list)


    if max_unclass_cbm:
        df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
        df = filter_samples_by_condition(df, f"Variable == 'Unclassified WBC' and Value < {max_unclass_cbm}")

    if min_wbc_mnl:
        df = filter_samples_by_condition(df, f"Variable == 'Total WBC' and Value >= {min_wbc_mnl}")

    if diff500:
        # ONLY filter for these samples for Aberrant Lymphocyte and Plasma Cell.
        # Other variables (like Segmented Neutrophil) keep all samples.
        df = filter_by_reference(
            df, 'flt_lists/500_WBC_mnl_cases.csv', include_rows=True,
            target_vars=['Aberrant Lymphocyte', 'Plasma Cell']
        )

    if crf_ssn == "post":
        df = filter_by_reference(df, 'flt_lists/pre_session_reviews.csv',)

    if aftr_2nd_ssn:
        df = filter_by_reference(df, 'flt_lists/pre_2nd_session_reviews.csv')

    if after_last_ssn:
        df = filter_by_reference(
            df, 'flt_lists/after_last_session.xlsx', include_rows=True,
            target_vars=['Aberrant Lymphocyte', 'Atypical Lymphocyte', 'LGL', 'Lymphocyte', 'Smudge Cell',
                         'RBC Agglutination', 'Rouleaux', 'Pelger Cell', 'Auer Rods', 'Hairy Cell',
                         'Aber&Atyp', 'Variant Lymphocyte', 'Aber&Hairy']
        )

    if lym3sites:
        lym_vars = ['Aberrant Lymphocyte', 'Atypical Lymphocyte', 'LGL', 'Lymphocyte', 'Hairy Cell',
                    'Aber&Atyp', 'Variant Lymphocyte', 'Aber&Hairy']
        lym_sites = ['BWH', 'SYN', 'TASMC']
        sites_ref_df = pd.DataFrame({'Site': lym_sites})
        df = filter_by_reference(df=df, filtering_source=sites_ref_df, include_rows=True, target_vars=lym_vars)


    df = df.query("Value!='--------'").copy()
    # Standardize names and tag the Arbitrator
    df['Investigator'] = careful_map(df['Investigator'], investigators_map)
    # Isolate the Arbitrator
    arb_mask = df['Investigator'] == 'Arbitrator'
    arb_df = df[arb_mask].copy()
    regular_df = df[~arb_mask].copy()

    # Dynamically assign roles (Rev1, Rev2, etc.) ONLY to regular reviewers
    regular_df = assign_dynamic_roles(regular_df, group_cols=['Site', 'SampleID'])

    # Calculate Mean Investigator ONLY on regular reviewers, returning the arbitrator lines afterwards
    regular_df = add_mean_investigator(regular_df, mthd=ref_arm, min_inv=min_inv)
    df = pd.concat([regular_df, arb_df], ignore_index=True)

    # Load arbitration rules (which samples/variables to override) and apply Quantitative Override (Overwrites 'Mean Investigator')
    arb_file = 'for_arbitration-including-candidates.csv' if no_arb_cands else 'for_arbitration.csv'
    arb_rules = read_to_df(f'flt_lists/{arb_file}', file_dir=os.getcwd())
    df = apply_arbitration_override(df, arb_df, arb_rules, metadata)

    # still need to deal with arbitration of binary and graded parameters
    # binary_vars = metadata.variable_groups["PLT morphology"] + metadata.variable_groups["RBC arrangement"] + metadata.variable_groups["WBC morphology"]
    binary_vars = ['PLT Clumps', 'Satellitisms']
    raw_grade_cond=lambda d: (
        d["Method"].isin([ref_arm])
        & d["Variable"].isin(binary_vars)
    )

    df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value nor grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)


    # cbm_file_name = f'all6_both_CBM_{cbm_version}{suffix}.csv'
    cbm_file_name = f'all6_both_CBM_{cbm_version}.csv'
    cbm_df = medium_pipe(cbm_file_name, None, test_arm, metadata, dir=r'raw/cbm_method_comparison', pre_cond=raw_cbm_cond)

    cbm_df['Investigator'] = test_arm
    cbm_df["Original_Investigator"] = test_arm

    if cbm_thresholding:
        # notice that values below threshold changed to 0% after calculation of percentages, so other percentages will not add up to 100%
        vars_to_thres = ['Hairy Cell', 'Aberrant Lymphocyte', 'Atypical Lymphocyte', 'LGL']
        thres = 5
        cbm_df.loc[(cbm_df['Variable'].isin(vars_to_thres)) & (cbm_df['Value'] < thres), 'Value'] = 0

    if rmv_tech_flgs:
        # maxmimal_rbc_fovs = 80
        maxmimal_rbc_fovs = 500
        cbm_df = filter_samples_by_condition(cbm_df, f"Variable == 'RBC Analysis Area' and Value <= {maxmimal_rbc_fovs}")
        minimal_mono_fovs = 500
        cbm_df = filter_samples_by_condition(cbm_df, f"Variable == 'Monolayer area' and Value >= {minimal_mono_fovs}")

        maximal_agran_plt = 40
        maximal_sphero = 3
        agran_sphero_conds = [
            f"Variable == 'Agranular Platelet' and Value > {maximal_agran_plt}",
            f"Variable == 'Spherocytes' and Value > {maximal_sphero}"
        ]
        cbm_df = filter_samples_by_multiple_conditions(cbm_df, agran_sphero_conds)


    if max_unclass_cbm:
        cbm_df = filter_samples_by_condition(cbm_df, f"Variable == 'Unclassified WBC' and Value <= {max_unclass_cbm}")


    # cbm_df = filter_by_condition(cbm_df, f"SampleID != '05394'")
    # cbm_df = filter_by_condition(cbm_df, f"SampleID != '05385'")


    """
    # RnA analysis areas that are incorrect
    df = filter_by_reference(
        df, 'flt_lists/500_WBC_mnl_cases.csv', include_rows=False,
        target_vars=['RBC Agglutination', 'Rouleaux']
    )
    """

    all_dfs = pd.concat([df, cbm_df])

    # cases to always remove - horrible slides, horrible scans, etc.
    rmv_df = read_to_df('flt_lists/fully_removed.csv', file_dir=os.getcwd())
    df = filter_by_reference(all_dfs, rmv_df)

    # cases of borderline quality - dirty, investigators' comments on quality, etc.
    if rmv_brd:
        rmv_df = read_to_df('flt_lists/low_borderline_quality.csv', file_dir=os.getcwd())
        df = filter_by_reference(all_dfs, rmv_df)

    methd_comp = MethodComparator(df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare'] + ['Aber&Hairy']
    # grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
    #     'PLT morphology'] + metadata.variable_groups['RBC arrangement']
    grades_to_test = ['PLT Clumps', 'Satellitisms']
    grades_to_print = grades_to_test + ['ScanID']
    # morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    morph_vals_to_test = ['PLT Clumps', 'Satellitisms']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also
    morph_vals_to_print = morph_vals_to_test + print_also

    if only_hard_params:
        vars_to_test = ['LGL', 'Atypical Lymphocyte', 'Aber&Hairy', 'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Blast', 'Plasma Cell', 'Lymphocyte']

    if diff500:
        vars_to_test = ['Plasma Cell']
    elif cbm_thresholding:
        vars_to_test = vars_to_thres


    if exprt_long:
        include_in_export = vals_to_print + grades_to_test + print_also
        df_long = methd_comp.df.query(f"Variable in @include_in_export and Investigator!='Mean Investigator'")[['SampleID', 'Site', 'Method', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long_all_revs.csv')
        df_long = methd_comp.df.query(f"Variable in @include_in_export and Investigator=='Mean Investigator'")[['SampleID', 'Site', 'Method', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long_final_values.csv')


    if inter:
        if bin_params:
            methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=binary_vars,
                                     dim_col='Investigator', comp_func='binary')
            if by_site:
                methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=binary_vars,
                                         dim_col='Investigator', split_by='Site', comp_func='binary')
            methd_comp.save_results(rf'results/mnl/{save_name}_bin_inter.csv', result_type="binary")

        methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=vars_to_test,
                                 dim_col='Investigator')
        if by_site:
            methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', variables=vars_to_test,
                                     dim_col='Investigator', split_by='Site')
        methd_comp.save_results(rf'results/mnl/{save_name}_reg_inter.csv')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg_inter.pdf')
        methd_comp.clean_calculations()

    if by_inv:
        methd_comp_by_rev = methd_comp.apply_to_df('query', f"Investigator!='Mean Investigator'", inplace=False)
        for rev in methd_comp_by_rev.df['Original_Investigator'].unique():
            if rev != test_arm and rev == rev:
                rev_row_filt = {"Original_Investigator": [rev, test_arm]}
                methd_comp_by_rev.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=vars_to_test,
                                                row_filters=rev_row_filt)
        methd_comp_by_rev.save_results(rf'results/mnl/{save_name}_by_rev_reg.csv')
        if plot_reg:
            methd_comp_by_rev.plot_all_regressions(f'results/mnl/{save_name}_by_rev_reg.pdf')

    if exprt_mtrx:
        methd_comp.export_comparison_matrix(
            out_path=fr'comp_tables/{save_name}_vals_all_inv.csv',
            row_identifiers=["Site", "SampleID"],
            comparison_dims=("Variable", "Method", "Investigator"),
            needed_vals=vals_to_print,
            needed_grades=['ScanID'])

        if min_wbc_mnl is False:
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
        if by_site:
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
    if by_site:
        methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')
    if plot_reg:
        methd_comp.plot_all_regressions(f'results/mnl/{save_name}_reg.pdf')
