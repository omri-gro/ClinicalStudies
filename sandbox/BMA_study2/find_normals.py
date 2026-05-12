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
    # version of this function that if a required arbitration was not performed yet takes whatever reviews we have
    arbitrators = _ensure_list(arbitrator)

    # Flag which Site-SampleID-Method combinations required arbitration
    df = df_raw.merge(df_arb.assign(is_arbitrated=True),
                      on=['Site', 'SampleID', 'Method'],
                      how='left').fillna({'is_arbitrated': False})

    # Identify which rows belong to an arbitrator
    df['is_arbitrator'] = df['Investigator'].isin(arbitrators)

    # Logic A: If not arbitrated, keep non-arbitrators
    cond_a = (~df['is_arbitrated']) & (~df['is_arbitrator'])

    # Logic B: If arbitrated and arbitrator result exists, keep arbitrator
    # Logic C: If arbitrated but no arbitrator result exists, keep the results (non-arbitrators)
    # To handle B and C, we check if an arbitrator actually performed work for that specific combo
    df['arb_performed'] = df.groupby(['Site', 'SampleID', 'Method'])['is_arbitrator'].transform('any')

    cond_b = df['is_arbitrated'] & df['arb_performed'] & df['is_arbitrator']
    cond_c = df['is_arbitrated'] & (~df['arb_performed'])

    # Combine and Filter
    df = df[cond_a | cond_b | cond_c].drop(columns=['is_arbitrated', 'is_arbitrator', 'arb_performed'])

    return df

if __name__ == "__main__":
    suffix = ''
    save_name = f'BMA_normals_search'
    meta_path = r'config_BMA.yaml'
    sites = ["OHSU", "HUP", "BWH"]
    arbitrators = ['Phil Raess', 'Olga Pozdnyakova', 'Christopher Hergott', 'OP', 'Arbitrator']

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    read_dir = os.path.join(cur_dir, 'raw')
    df_map = read_to_df(f'BMA_mapping.csv', file_dir=read_dir)

    metadata = MetadataBundle(meta_path)
    collect_dfs = []
    for site in sites:
        min_inv_site = 1

        ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir, recalc_diff=False)
        test_df = bma_prep_pipeline(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

        ref_df = min_inv_filt(ref_df, 'REF', min_inv=min_inv_site)
        test_df = min_inv_filt(test_df, 'TEST', min_inv=min_inv_site)

        # change all SampleIDs to the TEST Barcode based on site's mapping
        id_lookup = df_map.set_index('REF Barcode')['TEST Barcode']
        mapped_ids = ref_df['SampleID'].map(id_lookup)
        ref_df['SampleID'] = mapped_ids.fillna(ref_df['SampleID'])

        collect_dfs.append(pd.concat([ref_df, test_df]))

    all_dfs = pd.concat(collect_dfs)

    df_arb = read_to_df('to_arbitration.csv', file_dir=read_dir)
    all_dfs = removed_for_arbitration(all_dfs, df_arb, arbitrators)
    all_dfs = add_mean_investigator(all_dfs, mthd='REF', min_inv=0)
    all_dfs = add_mean_investigator(all_dfs, mthd='TEST', min_inv=0)

    all_dfs = add_pos_column(all_dfs, metadata)

    methd_comp = MethodComparator(all_dfs)
    methd_comp = methd_comp.apply_to_df('query', "Investigator=='Mean Investigator' and Method=='TEST'", inplace=False)

    rmv_file = 'flt_lists/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    ndc_vars_list = metadata.variable_groups['NDC'] + metadata.variable_groups['NDC lineage total']
    ndc_vars_list.remove('Unclassified')

    comp_table = methd_comp.export_comparison_matrix(needed_vars=ndc_vars_list,
                                                     value_col='Positive',
                                                     comparison_dims="Variable",
                                                     row_completeness="any")
    comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
    id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
    comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
    comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
    comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}_test_positivity.csv', index=True)

    comp_table = methd_comp.export_comparison_matrix(needed_vars=ndc_vars_list,
                                                     value_col='Value',
                                                     comparison_dims="Variable",
                                                     row_completeness="any")
    comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
    id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
    comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
    comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
    comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}_test_values.csv', index=True)

