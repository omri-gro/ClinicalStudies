import pandas as pd
import os
import sys
os.chdir('..')
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


def bma_prep_pipeline_alt(file_name, site, method, metadata, sheet_name='Sheet1', dir=None,
                          id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'], **kwargs):
    def stnd_bma_id(name):
        match = re.match(r"^([A-Za-z]*\d+)", str(name))
        return match.group(1) if match else name

    df = raw_bma_to_df(file_name, site, method, sheet_name, dir)
    # standardize column names
    df = stnd_names(df, metadata.alias_map)
    df['SampleID'] = df['SampleID'].apply(stnd_bma_id)

    if method == 'TEST':  # in future represent this as site rules
        df = calc_diff(df, metadata, diff_cells="NDC")
        # df = calc_diff(df, metadata, diff_cells="NDC", additional_cells="NDC-like")
        # df = calc_diff(df, metadata, diff_cells="NDC lineage")
    else:
        # print warning if WBCs in differential don't add up to ~100
        check_diff_sum(df, metadata, tolerance=5, diff_cells="NDC")
    df = pivot_long(df, id_vars=id_vars)

    # df = add_grade_column(df, metadata)
    # df = df.dropna(subset=["Value", "Grade"], how="all")

    df = df.dropna(subset="Value", how="all")
    return df


if __name__ == "__main__":
    save_name = 'two_sites'
    meta_path = r'config_BMA.yaml'
    sites = ["OHSU", "HUP"]
    arbitrators = ['Phil Raess', 'Olga Pozdnyakova']
    cur_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    read_dir = os.path.join(cur_dir, 'raw')
    df_map = read_to_df(f'BMA_mapping.csv', file_dir=read_dir)

    metadata = MetadataBundle(meta_path)
    metadata = MetadataBundle(meta_path)
    collect_dfs = []
    for site in sites:
        ref_df = bma_prep_pipeline_alt(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir)
        test_df = bma_prep_pipeline_alt(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

        df_arb = read_to_df('to_arbitration.csv', ref_df, file_dir=read_dir)
        ref_df = removed_for_arbitration(ref_df, df_arb, arbitrators)
        test_df = removed_for_arbitration(test_df, df_arb, arbitrators)

        # change all SampleIDs to the TEST Barcode based on site's mapping
        id_lookup = df_map.set_index('REF Barcode')['TEST Barcode']
        mapped_ids = ref_df['SampleID'].map(id_lookup)
        ref_df['SampleID'] = mapped_ids.fillna(ref_df['SampleID'])

        collect_dfs.append(pd.concat([ref_df, test_df]))

    all_dfs = pd.concat(collect_dfs)
    mk_df = all_dfs.query("Variable=='Megakaryocyte'")
    mk_df.to_csv(fr'{cur_dir}/comp_tables/mk_count.csv', index=False)

    all_dfs = pd.concat(collect_dfs)
    mk_df = all_dfs.query("Variable=='Megakaryocyte estimate'")
    mk_df.to_csv(fr'{cur_dir}/comp_tables/mk_estimate.csv', index=False)


