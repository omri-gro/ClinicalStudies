import os
import sys
import pandas as pd, matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr, kendalltau
from statsmodels.miscmodels.ordinal_model import OrderedModel

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df, add_mean_investigator, add_grade_column, add_pos_column, create_derived_variables_long
from pipelines import medium_pipe, clv_pipe
from itertools import *




if __name__ == "__main__":
    inter = False
    comp_with_cbm = True

    min_inv = 0  # False or number
    no_scrtch = False  # True to filter scratched slides out
    crf_ssn = 'all'  # 'all', 'pre' or 'post'
    rmv_arb = True

    intr_by_pair = False

    save_name = f'clv_cbm_{crf_ssn}-ssn_mininv-{min_inv}_no_scrtch-{no_scrtch}_arbrmv-{rmv_arb}'

    exprt_mtrx = True
    plot_reg = False
    inv_names_in_export = True  # if False investigators will appear as Rev1 and Rev2 only
    sites = ['BWH', 'LMU', 'TASMC']
    inv_map = {'Alina': 'Rev1', 'Alina KÃ¼pper': 'Rev1', 'Christine Lavoie': 'Rev1', 'Ebikebuna Rufus': 'Rev1', 'Sarah Pereira Rodrigues': 'Rev1',
               'Sladana': 'Rev2', 'Christopher Wright': 'Rev2', 'Thu Tran': 'Rev2', 'YAEL SAYEGH': 'Rev2',
               'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
    pair_map = {'Alina': 'Alina&Sladana', 'Alina KÃ¼pper': 'Alina&Sladana', 'Christine Lavoie': 'Christine&Chris', 'Ebikebuna Rufus': 'Ebi&Thu', 'Sarah Pereira Rodrigues': 'Sarah&Yael',
               'Sladana': 'Alina&Sladana', 'Christopher Wright': 'Christine&Chris', 'Thu Tran': 'Ebi&Thu', 'YAEL SAYEGH': 'Sarah&Yael',
               'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}
    pairs = ['Christine&Chris', 'Ebi&Thu', 'Sarah&Yael']
    by_list = pairs if intr_by_pair else sites
    test_arm = 'CBM'
    ref_arm = 'ClV'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = r'config.yaml'
    metadata = MetadataBundle(meta_path)

    vars_to_test = metadata.variable_groups['RBC morphology'] + metadata.variable_groups['PLT morphology'] + \
                   metadata.variable_groups['RBC combinations']
    vars_to_print = vars_to_test + ['TotalRBC'] + ['TotalPLT']
    id_vars_clv = ["SampleID", "Site", "Method", "FileName", 'Investigator']
    id_vars_cbm = ["SampleID", "Site", "Method", "FileName"]

    df_srcs_list = []
    for site in sites:
        # Currently cases with <min_inv investigators still appearing. Will be removed after mtrx_export.
        df = clv_pipe(f'{site}_ClV.csv', site, metadata, dir=r'raw/cbm_method_comparison',
                      min_inv=min_inv, mean_inv=True, drv_vars=False)
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

    if rmv_arb:
        rmv_file = r'flt_lists/slides_to_remove.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    if no_scrtch:
        rmv_file = r'flt_lists/scratched.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)


    if inter:
        int_save_name = f'clv_inter_{crf_ssn}-ssn_no_scrtch-{no_scrtch}_arbrmv-{rmv_arb}_bypair-{intr_by_pair}'

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



