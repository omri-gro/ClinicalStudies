import os
import sys
import pandas as pd
# import numpy as np
# import seaborn as sns
# from scipy.stats import f_oneway
# from matplotlib.backends.backend_pdf import PdfPages
# from scipy.stats import spearmanr, kendalltau
# from statsmodels.miscmodels.ordinal_model import OrderedModel

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df, add_mean_investigator, add_grade_column, add_pos_column, create_derived_variables_long, write_df_to_file
from pipelines import medium_pipe, clv_pipe
from itertools import *




if __name__ == "__main__":
    inter = False
    comp_with_cbm = False

    min_inv = 1  # False or number
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

    save_name = f'clv_cbm_{crf_ssn}-ssn_mininv-{min_inv}_no_scrtch-{no_scrtch}_brdrmv-{rmv_brd}'

    intr_by_pair = True

    exprt_long = True
    exprt_mtrx = True
    plot_reg = False
    inv_names_in_export = False  # if False investigators will appear as Rev1 and Rev2 only
    by_rev_comp = True
    rbc_agg_params = True

    sites = ['BWH', 'LMU', 'TASMC']
    inv_map = {'Alina': 'Rev1', 'Alina KÃ¼pper': 'Rev1', 'Christine Lavoie': 'Rev1', 'Ebikebuna Rufus': 'Rev1', 'Sarah Pereira Rodrigues': 'Rev1',
               'Nikolic Sladana': 'Rev2', 'Sladana': 'Rev2', 'Christopher Wright': 'Rev2', 'Thu Tran': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Sladana Nikolic': 'Rev2',
               'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
    names_map = {'Alina': 'Alina', 'Alina KÃ¼pper': 'Alina', 'Christine Lavoie': 'Christine', 'Ebikebuna Rufus': 'Ebi',
                 'Sarah Pereira Rodrigues': 'Sarah', 'Nikolic Sladana': 'Sladana', 'Sladana Nikolic': 'Sladana',
                 'Sladana': 'Sladana', 'Christopher Wright': 'Chris', 'Thu Tran': 'Thu', 'YAEL SAYEGH': 'Yael',
                 'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
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
    vars_to_test = metadata.variable_groups.get('RBC morphology', []) + metadata.variable_groups.get('PLT morphology', []) + \
                   metadata.variable_groups.get('RBC combinations', [])


    vars_to_print = vars_to_test + ['TotalRBC'] + ['TotalPLT']
    id_vars_clv = ["SampleID", "Site", "Method", "FileName", 'Investigator']
    id_vars_cbm = ["SampleID", "Site", "Method", "FileName"]

    df_srcs_list = []
    for site in sites:
        # Currently cases with < min_inv investigators still appearing. Will be removed after mtrx_export.

        extra_calcs = True if crf_ssn not in ['pre', 'post'] else False  # in these cases mean investigator and derived variables will be added later
        df = clv_pipe(f'{site}_ClV.csv', site, metadata, dir=r'raw/cbm_method_comparison',
                      min_inv=min_inv, mean_inv=extra_calcs, drv_vars=extra_calcs, only_mean=False)
        df_srcs_list.append(df)

    # converting clv results into MethodComparator to use filter_by_df - will not be necessary once MethodComparator is split into multiple objects
    df_clv = pd.concat(df_srcs_list)
    mthd_comp_clv = MethodComparator(df_clv)
    if crf_ssn in ['pre', 'post']:
        pre_ssn_file = r'flt_lists/pre_session_reviews.csv'
        keep_pre = (crf_ssn == 'pre')
        methd_comp_flt = mthd_comp_clv.filter_by_df(pre_ssn_file, include_rows=keep_pre)
        df_clv_flt = methd_comp_flt.df.copy()

        # requires calculating mean_inv again (as mean_inv were dropped when only pre-session were kept)
        df = add_mean_investigator(df_clv_flt, mthd=ref_arm, min_inv=min_inv)
        df = add_grade_column(df, metadata)
        df = add_pos_column(df, metadata)
        df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
        df = create_derived_variables_long(df, metadata)

        df_srcs_list = []
        df_srcs_list.append(df)

    df = medium_pipe(f'6sites_CBM.csv', None, 'CBM', metadata, dir=r'raw/cbm_method_comparison',
                     id_vars=id_vars_cbm)
    df["Investigator"] = test_arm
    df_srcs_list.append(df)
    all_dfs = pd.concat(df_srcs_list)

    methd_comp = MethodComparator(all_dfs)

    rmv_file = r'flt_lists/low_quality.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    # rmv_file = r'flt_lists/wrong_analysis_area.csv'
    # rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    # methd_comp = methd_comp.filter_by_df(rmv_df)


    # remove for arbitration only the morphologies in the arbitration categories
    arb_categories = {morph: ctgr for ctgr in ['RBC inclusions', 'RBC shape', 'RBC color'] for morph in metadata.variable_groups[ctgr]}
    methd_comp.df['Arbitration Category'] = methd_comp.df['Variable'].map(arb_categories)
    rmv_file = r'flt_lists/for_arbitration.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    if rmv_brd:
        rmv_file = r'flt_lists/slides_to_remove_borderline.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if no_scrtch:
        rmv_file = r'flt_lists/scratched.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)


    if exprt_long:
        arb_vars = metadata.variable_groups['RBC inclusions'] + metadata.variable_groups['RBC shape'] + metadata.variable_groups['RBC color']
        df_long = methd_comp.df.query(
            f"Variable in @arb_vars and Method=='{ref_arm}' and Investigator=='Mean Investigator'")[
            ['SampleID', 'Site', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        df_long['Investigator'] = df_long['Investigator'].map(inv_map)
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long.csv')

    if inter:
        int_save_name = f'clv_inter_{crf_ssn}-ssn_no_scrtch-{no_scrtch}_arbrmv-{rmv_brd}_bypair-{intr_by_pair}'

        intr_df = methd_comp.only_when_cond(f"Investigator!='Mean Investigator' and Method=='{ref_arm}'").df.copy()
        if intr_by_pair:
            intr_df['Site'] = intr_df['Investigator'].map(pair_map)
        intr_df['Investigator'] = intr_df['Investigator'].map(inv_map)
        intr_df['Method'] = intr_df['Investigator']
        inter_comp = MethodComparator(intr_df)
        if exprt_mtrx:
            # csv with column for each investigator
            inter_comp.export_comparison_matrix(out_path=fr'comp_tables/{int_save_name}.csv',
                                                row_identifiers=["Site", "SampleID"],
                                                comparison_dims=("Variable", "Method"),
                                                needed_vals=vars_to_print)

        inter_comp.batch_fit('Rev1', 'Rev2', vars_to_test)
        by_list = pairs if intr_by_pair else sites
        inter_comp.batch_fit('Rev1', 'Rev2', vars_to_test, site_filters=by_list)

        inter_comp.save_results(rf'results/clv/{int_save_name}_reg.csv')
        if plot_reg:
            inter_comp.plot_all_regressions(f'results/clv/{int_save_name}_reg.pdf')



    # --- comparison for each reviewer separately
    if by_rev_comp:
        methd_comp_by_rev = methd_comp.apply_to_df('query', f"Investigator!='Mean Investigator'", inplace=False)
        methd_comp_by_rev.df['Investigator'] = methd_comp_by_rev.df['Investigator'].map(names_map)
        for rev in methd_comp_by_rev.df['Investigator'].unique():
            rev_row_filt = {"Investigator": [rev, test_arm]}
            methd_comp_by_rev.batch_compare(levels_a=ref_arm, levels_b=test_arm, variables=vars_to_test, row_filters=rev_row_filt)
        methd_comp_by_rev.save_results(rf'results/clv/{save_name}_by_rev_reg.csv')
        if plot_reg:
            methd_comp_by_rev.plot_all_regressions(f'results/clv/{save_name}_by_rev_reg.pdf')

    if comp_with_cbm:
        if inv_names_in_export:
            methd_comp.df['Investigator'] = methd_comp.df['Investigator'].map(inv_map)


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



