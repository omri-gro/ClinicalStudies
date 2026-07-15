import os
import sys
import pandas as pd


sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import add_mean_investigator, create_derived_variables_long, assign_dynamic_roles
from pipelines import medium_pipe, clv_pipe


sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools import careful_map, apply_arbitration_override
from clinstudtools.core.metadata import MetadataBundle
from clinstudtools.utils import read_to_df, write_df_to_file
from clinstudtools.transforms import filter_by_reference, filter_samples_by_condition
from clinstudtools.preprocessing import add_grade_column, add_pos_column


if __name__ == "__main__":
    cbm_version = 'v319'  # currently v317 or v319
    inter = False
    comp_with_cbm = True
    sen_spec = True

    min_inv = 2  # False or number
    no_scrtch = False  # True to filter scratched slides out
    crf_ssn = 'all'  # 'all', 'pre' or 'post'
    rmv_brd = False

    """
    current filtering investigation
    min_inv=2 major improvement
    no_scrtch no major improvement
    crf_ssn not major difference between 'all' and 'post'
    rmv_brd no influence for clv (not enough borderline cases reviewed in ClV to make a difference)
    """

    save_name = f'clv_cbm_{crf_ssn}-ssn_mininv-{min_inv}_no_scrtch-{no_scrtch}_brdrmv-{rmv_brd}_{cbm_version}'

    intr_by_pair = False

    exprt_long = True
    exprt_mtrx = True
    plot_reg = True
    inv_names_in_export = False  # if False investigators will appear as Rev1 and Rev2 only
    by_rev_comp = False  # perform comparison for each reviewer separately
    rbc_agg_params = True  # parameters like Oval+Ellip, Acan+Echin
    with_morph_spec = True   # still need to implement

    sites = ['BWH', 'LMU', 'TASMC']

    inv_map = {
        # Standardize typos/variations
        'Christopher Wright': 'Chris',
        'Christine Lavoie': 'Christine',
        'Chris': 'Chris',
        'Christine': 'Christine',
        'Ebikebuna Rufus': 'Ebi',
        'Ebikebuna Rufus F.': 'Ebi',
        'Ebikebuna Rufus F': 'Ebi',
        'Ebi': 'Ebi',
        'Thu Tran': 'Thu',
        'THU TRAN': 'Thu',
        'Thu': 'Thu',
        'Alina KÃƒÂ¼pper': 'Alina',
        'Alina': 'Alina',
        'Sladana Nikolic': 'Sladana',
        'Nikolic Sladana': 'Sladana',
        'Sladana': 'Sladana',
        'Sarah Pereira Rodrigues': 'Sarah',
        'Sarah': 'Sarah',
        'YAEL SAYEGH': 'Yael',
        'Yael Sayegh': 'Yael',
        'Yael S': 'Yael',
        'YAEL ASYEGH': 'Yael',
        'Yael': 'Yael',

        # Explicitly tag Arbitrators
        'Dr. med. Weigand, Michael': 'Arbitrator',
        'Dan BENISTY': 'Arbitrator',

        # Preserve system/automated roles
        'ClV': 'ClV',
        'Mean Investigator': 'Mean Investigator'}

    names_map = inv_map
    pair_map = {'Alina': 'Alina&Sladana', 'Alina KÃ¼pper': 'Alina&Sladana', 'Christine Lavoie': 'Christine&Chris', 'Ebikebuna Rufus': 'Ebi&Thu', 'Sarah Pereira Rodrigues': 'Sarah&Yael',
               'Sladana': 'Alina&Sladana', 'Christopher Wright': 'Christine&Chris', 'Thu Tran': 'Ebi&Thu', 'YAEL SAYEGH': 'Sarah&Yael', 'Nikolic Sladana': 'Alina&Sladana',
               'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
    pairs = ['Christine&Chris', 'Ebi&Thu', 'Sarah&Yael', 'Alina&Sladana']
    by_list = pairs if intr_by_pair and inter else sites
    test_arm = 'CBM'
    ref_arm = 'ClV'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))

    if rbc_agg_params:
        meta_path = r'special_config_files/config_agg.yaml'
    else:
        meta_path = r'config.yaml'
    metadata = MetadataBundle(meta_path)

    # from previous attempts to quantify PLT morphologies with ClV
    vars_to_test = metadata.variable_groups.get('RBC morphology', []) + metadata.variable_groups.get('RBC combinations', [])


    vars_to_print = vars_to_test + ['TotalRBC'] + ['TotalPLT']
    id_vars_clv = ["SampleID", "Site", "Method", "FileName", 'Investigator']
    id_vars_cbm = ["SampleID", "Site", "Method", "FileName"]

    df_srcs_list = []
    for site in sites:
        # Currently cases with < min_inv investigators still appearing. Will be removed after mtrx_export.

        df = clv_pipe(f'{site}_ClV.csv', site, metadata, dir=r'raw/cbm_method_comparison',
                      min_inv=min_inv, mean_inv=False, drv_vars=False, only_mean=False)
        df_srcs_list.append(df)

    df = pd.concat(df_srcs_list)

    if crf_ssn in ['pre', 'post']:
        keep_pre = (crf_ssn == 'pre')
        df = filter_by_reference(df, r'flt_lists/pre_session_reviews.csv', include_rows=keep_pre)

    # Standardize names and tag the Arbitrator
    df['Investigator'] = careful_map(df['Investigator'], inv_map)
    # Isolate the Arbitrator
    arb_mask = df['Investigator'] == 'Arbitrator'
    arb_df = df[arb_mask].copy()
    regular_df = df[~arb_mask].copy()

    # Dynamically assign roles (Rev1, Rev2, etc.) ONLY to regular reviewers
    regular_df = assign_dynamic_roles(regular_df, group_cols=['Site', 'SampleID'])
    regular_df = add_mean_investigator(regular_df, mthd=ref_arm, min_inv=min_inv)
    df = pd.concat([regular_df, arb_df], ignore_index=True)

    # Load arbitration rules (which samples/variables to override) and apply Quantitative Override (Overwrites 'Mean Investigator')
    arb_rules = read_to_df('flt_lists/for_arbitration.csv', file_dir=os.getcwd())
    df_clv = apply_arbitration_override(df, arb_df, arb_rules, metadata)

    df_clv = add_grade_column(df_clv, metadata)
    df_clv = add_pos_column(df_clv, metadata)
    df_clv = df_clv.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    df_clv = create_derived_variables_long(df_clv, metadata)

    df_cbm = medium_pipe(f'all6_RGB_CBM_{cbm_version}.csv', None, 'CBM', metadata, dir=r'raw/cbm_method_comparison',
                     id_vars=id_vars_cbm, check_wbc_diff=False)
    df_cbm = add_grade_column(df_cbm, metadata)
    df_cbm = add_pos_column(df_cbm, metadata)
    df_cbm = df_cbm.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    df_cbm = create_derived_variables_long(df_cbm, metadata)
    df_cbm["Investigator"] = test_arm
    df_cbm["Original_Investigator"] = test_arm
    all_dfs = pd.concat([df_clv, df_cbm])


    df = filter_by_reference(all_dfs, r'flt_lists/low_quality.csv')

    # df = filter_by_reference(df, r'flt_lists/wrong_analysis_area.csv')

    if rmv_brd:
        df = filter_by_reference(df, r'flt_lists/slides_to_remove_borderline.csv')

    if no_scrtch:
        df = filter_by_reference(df, r'flt_lists/scratched.csv')


    methd_comp = MethodComparator(df)

    if exprt_long:
        arb_vars = metadata.variable_groups['RBC inclusions'] + metadata.variable_groups['RBC shape'] + metadata.variable_groups['RBC color']
        df_long = methd_comp.df.query(
            f"Variable in @arb_vars and Method=='{ref_arm}' and Investigator=='Mean Investigator'")[
            ['SampleID', 'Site', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        df_long['Investigator'] = careful_map(df_long['Investigator'], inv_map)
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long.csv')

    if inter:
        int_save_name = f'clv_inter_{crf_ssn}-ssn_no_scrtch-{no_scrtch}_arbrmv-{rmv_brd}_bypair-{intr_by_pair}'

        intr_df = methd_comp.only_when_cond(f"Investigator!='Mean Investigator' and Method=='{ref_arm}'").df.copy()
        if intr_by_pair:
            intr_df['Site'] = careful_map(intr_df['Investigator'], pair_map)
        # intr_df['Investigator'] = careful_map(intr_df['Investigator'], inv_map)
        intr_df['Method'] = intr_df['Investigator']
        inter_comp = MethodComparator(intr_df)
        if exprt_mtrx:
            # csv with column for each investigator
            inter_comp.export_comparison_matrix(out_path=fr'comp_tables/{int_save_name}.csv',
                                                row_identifiers=["Site", "SampleID"],
                                                comparison_dims=("Variable", "Method"),
                                                needed_vals=vars_to_print,
                                                row_completeness='any')

        inter_comp.batch_fit('Rev1', 'Rev2', vars_to_test)
        by_list = pairs if intr_by_pair else sites
        inter_comp.batch_fit('Rev1', 'Rev2', vars_to_test, site_filters=by_list)

        inter_comp.save_results(rf'results/clv/{int_save_name}_reg.csv')
        if plot_reg:
            inter_comp.plot_all_regressions(f'results/clv/{int_save_name}_reg.pdf')




    # --- comparison for each reviewer separately
    if by_rev_comp:
        methd_comp_by_rev = methd_comp.apply_to_df('query', f"Investigator!='Mean Investigator'", inplace=False)
        for rev in methd_comp_by_rev.df['Original_Investigator'].unique():
            if rev != test_arm and rev == rev:
                rev_row_filt = {"Original_Investigator": [rev, test_arm]}
                methd_comp_by_rev.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=vars_to_test,
                                                row_filters=rev_row_filt)
        methd_comp_by_rev.save_results(rf'results/clv/{save_name}_by_rev_reg.csv')
        # if plot_reg:
        #     methd_comp_by_rev.plot_all_regressions(f'results/clv/{save_name}_by_rev_reg.pdf')


    if sen_spec:
        bin_comp = methd_comp.apply_to_df('query', f"Investigator=='Mean Investigator' or Investigator=='{test_arm}'", inplace=False)
        bin_comp.batch_compare(levels_a=ref_arm, levels_b=test_arm,
                               variables=vars_to_test, comp_func='sen_spe', cis=False)
        bin_comp.batch_compare(levels_a=ref_arm, levels_b=test_arm,
                               variables=vars_to_test, comp_func='sen_spe', split_by='Site', cis=False)
        bin_comp.save_results(rf'results/clv/{save_name}_sen_spe.csv', result_type="sen_spe")


    if comp_with_cbm:
        if inv_names_in_export:
            methd_comp.df['Investigator'] = careful_map(methd_comp.df['Investigator'], inv_map)

        if exprt_mtrx:
            # csv with column for each investigator
            methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}_all_inv.csv',
                                                row_identifiers=["Site", "SampleID"],
                                                comparison_dims=("Variable", "Method", "Investigator"),
                                                needed_vals=vars_to_print,
                                                needed_grades=['ScanID'])

        # comparison of CBM will be only with MeanInvestigator
        methd_comp = methd_comp.apply_to_df('query', f"Investigator=='Mean Investigator' or Investigator=='{test_arm}'", inplace=False)

        if exprt_mtrx:
            # csv with column with only mean investigator and test arm
            methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}_mean_inv.csv',
                                                row_identifiers=["Site", "SampleID"],
                                                comparison_dims=("Variable", "Method"),
                                                needed_vals=vars_to_print,
                                                needed_grades=['ScanID'])

        methd_comp.batch_fit(ref_arm, test_arm, vars_to_test)
        methd_comp.batch_fit(ref_arm, test_arm, vars_to_test, site_filters=by_list)

        methd_comp.save_results(rf'results/clv/{save_name}_reg.csv')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/clv/{save_name}_reg.pdf')



