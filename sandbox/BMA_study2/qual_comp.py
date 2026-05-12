import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle, raw_bma_to_df, _ensure_list
from sandbox import *
from pipelines import bma_prep_pipeline

def removed_for_arbitration(df_raw, df_arb, arbitrator):
    arbitrators = _ensure_list(arbitrator)

    # check which samples went to arbitration
    df = pd.merge(df_raw, df_arb, on=['Site', 'SampleID', 'Method'], how='left', indicator=True)

    # keep reviews which were not sent to arbitration and reviewed by regular reviewer, or ones sent to arbitration and reviewed by arbitrator
    df = df[((df['_merge'] == 'left_only') & ~(df['Investigator'].isin(arbitrators))) | ((df['_merge'] == 'both') & (df['Investigator'].isin(arbitrators)))]

    return df

if __name__ == "__main__":
    suffix = '_additional_omr-based_arbitrations'
    save_name = f'sample_qual{suffix}'
    meta_path = r'config_BMA.yaml'
    sites = ["OHSU", "HUP", "BWH"]
    arbitrators = ['Phil Raess', 'Olga Pozdnyakova', 'Christopher Hergott', 'OP', 'Arbitrator']

    compare_methods = True
    raw_dss = False
    inter = False
    inter_to_include_arbitrated = False

    exprt_mtrx = True
    exprt_long = True
    plot_reg = True
    keep_names = False  # use investigators' full names - creates very wide 'all investigators' comparison matrix if True

    other_removed = True  # for filtering out samples for side analysis
    rmv_unclass = False  # re-calculate differential as if all unclassified were moved to dirt/other
    pooled_params = True  # analyze for pooled parameters like Erythroblast&BasophilicNormoblast - only when rmv_unclass False

    investigators_map = {'Todd Williams': 'Rev1', 'Wei Xie': 'Rev2', 'Phil Raess': 'Arbitrator',
                         'TW': 'Rev1', 'WX': 'Rev2', 'PR': 'Arbitrator',
                         'AB': 'Rev1', 'AS': 'Rev2', 'DL': 'Rev3', 'OP': 'Arbitrator',
                         'Adam Bagg': 'Rev1', 'Annapurna Saksena': 'Rev2', 'Dorottya Laczko': 'Rev3', 'Olga Pozdnyakova': 'Arbitrator',
                         'Elizabeth Morgan': 'Rev1', 'Habibe Kurt': 'Rev2', 'Robert Hasserjian': 'Rev3',
                         'Sam Sadigh': 'Rev4', "Megan Fitzpatrick": 'Rev5', "Vignesh Shanmugam": 'Rev6',
                         "Christopher Hergott-rev": 'Rev7', "Christopher Hergott": 'Arbitrator',
                         'Rev1': 'Rev1', 'Rev2': 'Rev2',
                         'Mean Investigator': 'Mean Investigator'}
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    read_dir = os.path.join(cur_dir, 'raw')
    df_map = read_to_df(f'BMA_mapping.csv', file_dir=read_dir)

    if rmv_unclass:
        meta_path = r'config_BMA_no_unclass.yaml'
    elif pooled_params:
        meta_path = r'config_BMA_pool.yaml'

    metadata = MetadataBundle(meta_path)
    collect_dfs = []
    for site in sites:
        min_inv_site = 2

        ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir, recalc_diff=rmv_unclass)
        test_df = bma_prep_pipeline(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

        ref_df = min_inv_filt(ref_df, 'REF', min_inv=min_inv_site)
        test_df = min_inv_filt(test_df, 'TEST', min_inv=min_inv_site)

        # df_arb = df_arb[~((df_arb['Site'] == 'HUP') & (df_arb['Method'] == 'TEST'))]

        """
        ref_df = removed_for_arbitration(ref_df, df_arb, arbitrators)
        test_df = removed_for_arbitration(test_df, df_arb, arbitrators)

        ref_df = add_mean_investigator(ref_df, mthd='REF', min_inv=0)
        test_df = add_mean_investigator(test_df, mthd='TEST', min_inv=0)
        """

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

        collect_dfs.append(df_dss)

    all_dfs = pd.concat(collect_dfs)
    if not keep_names:
        all_dfs['Investigator'] = all_dfs['Investigator'].map(investigators_map)
    all_dfs_all_inv = all_dfs

    df_arb = read_to_df('to_arbitration.csv', file_dir=read_dir)

    if not inter_to_include_arbitrated:
        all_dfs_all_inv = removed_for_arbitration(all_dfs_all_inv, df_arb, arbitrators)
    methd_comp_all_inv = MethodComparator(all_dfs_all_inv)

    all_dfs = removed_for_arbitration(all_dfs, df_arb, arbitrators)
    all_dfs = add_pos_column(all_dfs, metadata)

    quality_params = metadata.variable_groups['quality']
    qual_df = all_dfs.query("Variable in @quality_params")

    if exprt_long:
        qual_df.to_csv(fr'{cur_dir}/comp_tables/{save_name}_long.csv', index=False)

    methd_comp = MethodComparator(qual_df)

    methd_comp.batch_compare(levels_a='REF',
                             levels_b='TEST',
                             variables=quality_params,
                             comp_func='confusion_matrix',
                             id_cols=('SampleID', 'Investigator'))

    methd_comp.batch_compare(levels_a='REF',
                             levels_b='TEST',
                             variables=quality_params,
                             comp_func='confusion_matrix',
                             split_by='Site',
                             id_cols=('SampleID', 'Investigator'))

    methd_comp.save_results(fr'{cur_dir}/results/{save_name}_conf_mat.csv', result_type='confusion_matrix')
    methd_comp.save_results(fr'{cur_dir}/results/{save_name}_conf_mat_vis.xlsx', result_type='matrix_visual')

    methd_comp.clean_calculations()

    methd_comp.batch_compare(levels_a='REF',
                             levels_b='TEST',
                             variables=quality_params,
                             comp_func='confusion_matrix',
                             id_cols=('SampleID', 'Investigator'))
