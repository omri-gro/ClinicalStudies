import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle, raw_bma_to_df
from sandbox import *
from pipelines import bma_prep_pipeline


# def bma_df_pipeline(paths: dict, metadata: sb.MetadataBundle, dir=None, measurement_col='Value',
#                     more_id_vars=None):


def removed_for_arbitration(df_raw, df_arb, arbitrator):
    # check which samples went to arbitration
    df = pd.merge(df_raw, df_arb, on=['Site', 'SampleID', 'Method'], how='left', indicator=True)

    # keep reviews which were not sent to arbitration, or ones sent to arbitration and reviewed by arbitrator
    df = df[(df['_merge'] == 'left_only') | ((df['_merge'] == 'both') & (df['Investigator'] == arbitrator))]

    return df


if __name__ == "__main__":
    meta_path = r'config_BMA.yaml'
    site = 'OHSU'
    arbitrator = 'Phil Raess'
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    read_dir = os.path.join(cur_dir, 'raw')

    metadata = MetadataBundle(meta_path)
    srcs = {(site, 'REF'): f'{site}_CRF_REF.csv',
            (site, 'TEST'): f'{site}_CRF_TEST.csv'}

    ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir)
    test_df = bma_prep_pipeline(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

    df_arb = read_to_df('to_arbitration.csv', ref_df, file_dir=read_dir)
    ref_df_filtered = removed_for_arbitration(ref_df, df_arb, arbitrator)
    ref_df_filtered = removed_for_arbitration(test_df, df_arb, arbitrator)

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=read_dir, bma=True)
    ndc_vars_list = metadata.variable_groups['NDC']
    needed_rows = methd_comp.df[methd_comp.df['Variable'].isin(ndc_vars_list)]

    # calculate means for variable-identifier combinations that have exactly 2 rows
    variable_identifier = ['SampleID', 'Site', 'Method', 'Variable']
    m = needed_rows.groupby(variable_identifier)['Value'].transform('size').eq(2)
    out = needed_rows.loc[m].groupby(variable_identifier, as_index=False)['Value'].mean()
    methd_comp.df = out

    methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list)
    methd_comp.save_results(rf'results/{site}_bma_reg.csv')
    methd_comp.plot_all_regressions(f'results/{site}_bma_reg.pdf')

    print(methd_comp)
    print()


