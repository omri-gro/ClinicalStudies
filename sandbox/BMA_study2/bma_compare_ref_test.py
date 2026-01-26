import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle, raw_bma_to_df, _ensure_list
from sandbox import *
from pipelines import bma_prep_pipeline


# def bma_df_pipeline(paths: dict, metadata: sb.MetadataBundle, dir=None, measurement_col='Value',
#                     more_id_vars=None):


def removed_for_arbitration(df_raw, df_arb, arbitrator):
    arbitrators = _ensure_list(arbitrator)

    # check which samples went to arbitration
    df = pd.merge(df_raw, df_arb, on=['Site', 'SampleID', 'Method'], how='left', indicator=True)

    # keep reviews which were not sent to arbitration, or ones sent to arbitration and reviewed by arbitrator
    df = df[(df['_merge'] == 'left_only') | ((df['_merge'] == 'both') & (df['Investigator'].isin(arbitrators)))]

    return df


if __name__ == "__main__":
    save_name = 'two_sites'
    meta_path = r'config_BMA.yaml'
    sites = ["OHSU", "HUP"]
    arbitrators = ['Phil Raess', 'Olga Pozdnyakova']

    compare_methods = True
    inter = False
    raw_dss = False

    exprt_mtrx = False
    plot_reg = False
    keep_names = True  # use investigators' full names - creates very long 'all investigators' comparison matrix if True

    investigators_map = {'Todd Williams': 'Rev1', 'Wei Xie': 'Rev2', 'Phil Raess': 'Arbitrator',
                         'AB': 'Rev1', 'AS': 'Rev2', 'DL': 'Rev3', 'OP': 'Arbitrator',
                         'Adam Bagg': 'Rev1', 'Annapurna Saksena': 'Rev2', 'Dorottya Laczko': 'Rev3', 'Olga Pozdnyakova': 'Arbitrator',
                         'Mean Investigator': 'Mean Investigator'}
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    read_dir = os.path.join(cur_dir, 'raw')
    df_map = read_to_df(f'BMA_mapping.csv', file_dir=read_dir)

    metadata = MetadataBundle(meta_path)
    collect_dfs = []
    for site in sites:
        ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir)
        test_df = bma_prep_pipeline(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

        df_arb = read_to_df('to_arbitration.csv', ref_df, file_dir=read_dir)
        ref_df = removed_for_arbitration(ref_df, df_arb, arbitrators)
        test_df = removed_for_arbitration(test_df, df_arb, arbitrators)

        # need to find way add_mean_investigator's min_inv could handle arbitrator-only samples,
        # maybe split add_mean_investigator into 2-3 steps, or integrate with removed_for_arbitration
        # for now can do this by removing only 1 reviewer cases beforehand
        ref_df = add_mean_investigator(ref_df, mthd='REF', min_inv=2)
        test_df = add_mean_investigator(test_df, mthd='TEST', min_inv=2)

        # change all SampleIDs to the TEST Barcode based on site's mapping
        id_lookup = df_map.set_index('REF Barcode')['TEST Barcode']
        mapped_ids = ref_df['SampleID'].map(id_lookup)
        ref_df['SampleID'] = mapped_ids.fillna(ref_df['SampleID'])

        collect_dfs.append(pd.concat([ref_df, test_df]))

    if raw_dss:
        df_dss = read_to_df(f'raw_DSS.csv', file_dir=read_dir)
        df_dss = stnd_names(df_dss, metadata.alias_map)
        df_dss["Method"] = 'DSS'
        df_dss["FileName"] = os.path.basename(f'raw_DSS.csv')
        df_dss.columns = df_dss.columns.str.strip()
        df_dss = calc_diff(df_dss, metadata, diff_cells="NDC")
        id_vars = ["SampleID", "Site", "Method", "FileName", 'Investigator']
        df_dss = pivot_long(df_dss, id_vars=id_vars)
        df_dss = add_grade_column(df_dss, metadata)
        df_dss = df_dss.dropna(subset=["Value", "Grade"], how="all")
        df_dss = create_derived_variables_long(df_dss, metadata)
        df_dss = add_mean_investigator(df_dss, mthd='DSS', min_inv=0)

        # currently BWH only relevant when analysing DSS; add it to main for loop once BWH digital reviews arrive
        site = 'BWH'
        ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir)
        df_arb = read_to_df('to_arbitration.csv', ref_df, file_dir=read_dir)
        ref_df = removed_for_arbitration(ref_df, df_arb, arbitrators)
        ref_df = add_mean_investigator(ref_df, mthd='REF', min_inv=2)
        # id_lookup = df_map.set_index('REF Barcode')['TEST Barcode']
        # mapped_ids = ref_df['SampleID'].map(id_lookup)
        # ref_df['SampleID'] = mapped_ids.fillna(ref_df['SampleID'])

        collect_dfs.append(pd.concat([ref_df, df_dss]))


    all_dfs = pd.concat(collect_dfs)
    all_dfs = add_pos_column(all_dfs, metadata)

    if not keep_names:
        all_dfs['Investigator'] = all_dfs['Investigator'].map(investigators_map)

    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator'", inplace=False)


    ndc_vars_list = metadata.variable_groups['NDC'] + metadata.variable_groups['NDC lineage total']
    needed_rows = methd_comp.df[methd_comp.df['Variable'].isin(ndc_vars_list)]
    grade_vars_list = metadata.variable_groups['grade']

    if exprt_mtrx:
        comp_table = methd_comp.export_comparison_matrix(needed_vals=ndc_vars_list,
                                                         comparison_dims=("Variable", "Method", "Investigator"),
                                                         row_completeness="none")
        comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
        id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
        comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
        comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
        comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}.csv', index=True)

    comp_table = methd_comp_all_inv.export_comparison_matrix(needed_vals=ndc_vars_list, needed_grades=grade_vars_list,
                                                             comparison_dims=("Variable", "Method", "Investigator"),
                                                             row_completeness="none")
    comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
    id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
    comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
    comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
    comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}_all_investigators.csv', index=True)

    # main method comparison regressions + biases
    if compare_methods:
        methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list)
        methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list, site_filters=sites)
        methd_comp.calc_all_biases(metadata.crit_points)
        methd_comp.save_results(rf'results/{save_name}_bma_reg.csv')
        methd_comp.save_results(rf'results/{save_name}_bias.xlsx', result_type='bias')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/{save_name}_bma_reg.pdf')


    # inter-investigator
    if inter:
        if keep_names:
            levels_a = ['Todd Williams', 'AB', 'AS', 'Adam Bagg', 'Annapurna Saksena']
            levels_b = ['Wei Xie', 'AS', 'DL', 'Adam Bagg', 'Dorottya Laczko']
        else:
            levels_a = ['Rev1']
            levels_b = ['Rev2', 'Rev3']
        methd_comp_all_inv.batch_compare(levels_a=levels_a, levels_b=levels_b, variables=ndc_vars_list,
                                 dim_col='Investigator', split_by=['Method', 'Site'])
        methd_comp_all_inv.save_results(rf'results/{save_name}_inter.csv')
        if plot_reg:
            methd_comp_all_inv.plot_all_regressions(f'results/{save_name}_inter.pdf')


    # raw dss vs reference regression
    if raw_dss:
        # methd_comp.clean_calculations()
        methd_comp.batch_compare(levels_a=['REF', 'TEST'], levels_b='DSS', variables=ndc_vars_list)
        methd_comp.batch_compare(levels_a=['REF', 'TEST'], levels_b='DSS', variables=ndc_vars_list, split_by='Site')
        methd_comp.save_results(rf'results/{save_name}_raw_dss.csv')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/{save_name}_raw_dss.pdf')




