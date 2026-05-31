import pandas as pd
import os
import sys
import numpy as np
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools import safe_pivot
from clinstudtools.utils import write_df_to_file

# can later move to sandbox or clinstudtools
def assign_dynamic_roles(df, group_cols=['Site', 'SampleID'], inv_col='Investigator', prefix='Rev'):
    """
    Dynamically assigns relative reviewer roles per sample to avoid hardcoded mappings.
    Optionally preserves the original name in a new column for traceability.
    """
    df = df.copy()

    # Save the original name just in case you need it for debugging/arbitration
    if 'Original_Investigator' not in df.columns:
        df['Original_Investigator'] = df[inv_col]

    # Assign the dynamic role using a list comprehension to avoid NumPy string errors
    df[inv_col] = (
        df.groupby(group_cols)[inv_col]
        .transform(lambda x: [f"{prefix}{i + 1}" for i in pd.factorize(x)[0]])
    )

    return df


if __name__ == "__main__":
    # read long table containing cases reviewed by two reviewers
    long_file = r'mnl_allssn__mininv2_long_no_renaming.csv'
    long_df = read_to_df(long_file, file_dir=r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\comp_tables')
    long_df.sort_values(by=['Site', 'SampleID', 'Variable', 'Investigator'])

    mnl_df = long_df[long_df['Investigator'] != 'CBM']
    mnl_df['Method'] = 'Reference'
    mnl_df = assign_dynamic_roles(mnl_df, group_cols=['Site', 'SampleID'])
    mnl_df = min_inv_filt(mnl_df, 'Reference', min_inv=2, exact=True)  # keeps only (SampleID, Site, Variable) combinations with with exactly 2 Reference method rows

    #  comparing numeric variables
    meta_path = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\config.yaml'
    metadata = MetadataBundle(meta_path)
    diff_vars = metadata.variable_groups['WBC diff']

    num_df = mnl_df.dropna(subset=['Value'])
    num_df = num_df[num_df['Variable'].isin(diff_vars)]
    num_df = num_df[['SampleID', 'Site', 'Investigator', 'Variable', 'Value']]

    num_pairs_df = safe_pivot(df=num_df,
                              index=['Site', 'SampleID', 'Variable'],
                              columns=['Investigator'],
                              values='Value')  # regular pivots with warnings on duplicates

    num_pairs_df['Mean'] = num_pairs_df[['Rev1', 'Rev2']].mean(axis=1)
    num_pairs_df = num_pairs_df.reset_index()

    #  Add 'Mean Total WBC'
    mean_wbc = mnl_df[mnl_df['Variable'] == 'Total WBC'].groupby(['Site', 'SampleID'])['Value'].mean().reset_index()
    mean_wbc.rename(columns={'Value': 'Mean Total WBC'}, inplace=True)
    num_pairs_df = num_pairs_df.merge(mean_wbc, on=['Site', 'SampleID'], how='left')

    # check 99% CIs
    booschlo_df = pd.read_csv(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\special_config_files\booschlo_ranges.csv')
    booschlo_df = booschlo_df.sort_values('Mean')
    num_pairs_df['Row_Mean'] = np.ceil(num_pairs_df[['Rev1', 'Rev2']].mean(axis=1) * 2) / 2  # closest CI from csv
    num_pairs_df = num_pairs_df.sort_values('Row_Mean')

    num_pairs_df = pd.merge_asof(
        num_pairs_df,
        booschlo_df[['Mean', 'Low99', 'High99']],
        left_on='Row_Mean',
        right_on='Mean',
        direction='nearest'
    )

    num_pairs_df['BothIn99'] = (
            (num_pairs_df['Rev1'] >= num_pairs_df['Low99']) & (num_pairs_df['Rev1'] <= num_pairs_df['High99']) &
            (num_pairs_df['Rev2'] >= num_pairs_df['Low99']) & (num_pairs_df['Rev2'] <= num_pairs_df['High99'])
    )

    # also compare with CBM result if exists
    cbm_df = long_df[long_df['Investigator'] == 'CBM'][['Site', 'SampleID', 'Variable', 'Value']]
    cbm_df = cbm_df.rename(columns={'Value': 'CBM'})
    num_pairs_df = num_pairs_df.merge(cbm_df, on=['Site', 'SampleID', 'Variable'], how='left')
    num_pairs_df['CBM_in99'] = (
            (num_pairs_df['CBM'] >= num_pairs_df['Low99']) &
            (num_pairs_df['CBM'] <= num_pairs_df['High99'])
    )
    num_pairs_df['CBM_in99'] = num_pairs_df['CBM_in99'].where(num_pairs_df['CBM'].notna(), False)  # ensure CBM_in99 false if no CBM value known

    # dataframe only for disagreements
    num_dis_df = num_pairs_df[~num_pairs_df['BothIn99']].copy()

    # count disagreements
    disagreements_count = num_dis_df[~num_dis_df['CBM_in99']].groupby(['Site', 'SampleID']).size().reset_index(
        name='Disagreements count')
    num_dis_df = num_dis_df.merge(disagreements_count, on=['Site', 'SampleID'], how='left')
    # Fill NaNs with 0 for combinations that had no disagreements
    num_dis_df['Disagreements count'] = num_dis_df['Disagreements count'].fillna(0)

    # remove cases already checked for arbitration
    already_checked_for_arbitration_df = read_to_df('already_checked_for_arbitration.csv', file_dir=r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\flt_lists')
    checked_keys = already_checked_for_arbitration_df[['Site', 'SampleID']].drop_duplicates()
    checked_keys['Already checked'] = True
    num_dis_df = num_dis_df.merge(checked_keys, on=['Site', 'SampleID'], how='left')
    num_dis_df['Already checked'] = num_dis_df['Already checked'].fillna(False)

    # find potential arbitration cases
    ptnt_num_arb = num_dis_df[
        (~num_dis_df['Already checked']) &
        (~num_dis_df['CBM_in99']) &
        (num_dis_df['Rev1'].notna()) &
        (num_dis_df['Rev2'].notna())
    ].copy()

    ptnt_num_arb['CBM_exists'] = ptnt_num_arb['CBM'].notna()

    ptnt_num_arb = ptnt_num_arb.sort_values(
        by=['CBM_exists', 'Disagreements count', 'Site', 'SampleID', 'Variable'],
        ascending=[False, False, False, False, False]
    )

    ptnt_num_arb = ptnt_num_arb.drop(columns=['Mean_x', 'Mean_y', 'BothIn99', 'CBM_in99', 'Already checked', 'CBM_exists'])
    write_df_to_file(ptnt_num_arb, r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\results\side_games\potential_arbitration_mnl_3105.csv')



