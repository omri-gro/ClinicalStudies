import pandas as pd

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\clinstudtools')
from table_integrity import robust_dup

if __name__ == "__main__":
    mthd = 'TEST'  # for test, remember to change name of 'Total' raw to 'Total nucleated cells'
    site = 'HUP'

    investigators = {'Wei Xie': 'User1', 'Todd Williams': 'User2',
                     'Elizabeth Morgan': 'User1', 'Habibe Kurt': 'User2', 'Robert Hasserjian': 'User3', 'Sam Sadigh': 'User4',
                     'DL': 'User1', 'AB': 'User2', 'AS': 'User2',
                     "Phil Raess": 'Arbitrator', "Christopher Hergott": 'Arbitrator', 'OP': 'Arbitrator', "Olga Pozdnyakova": 'Arbitrator',
                     'ev1': 'User1', 'ev2': 'User2', 'arb': 'Arbitrator'}

    arbitrators = ["Phil Raess", "Christopher Hergott", "OP", "Olga Pozdnyakova"]
    # suffixes = ["rev1", "rev2", "arb"]

    # investigators = {'Elizabeth Morgan': 'User1', 'Habibe Kurt': 'User2', 'Robert Hasserjian': 'User1', 'Sam Sadigh': 'User2'}

    typed_file = rf'raw\{site}_CRF_{mthd}.csv'
    df = pd.read_csv(typed_file)

    # Identify empty specimen columns (entire column is empty below the header)
    empty_columns = [col for col in df.columns[1:] if df[col].iloc[2:].isna().all()]
    # Print message and remove empty columns if any exist
    if empty_columns:
        print(f"Removing empty specimen columns: {empty_columns}")
        df = df.drop(columns=empty_columns)

    # drop unnecessary rows
    users_row = df.iloc[0].map(investigators)
    df_cleaned = df.iloc[2:].reset_index(drop=True)
    df_cleaned = df_cleaned[df_cleaned['Case ID'] != 'Total']

    # Rename the first column to 'Parameter' for clarity
    df_cleaned.rename(columns={df_cleaned.columns[0]: 'Parameter'}, inplace=True)

    # Convert to long format (melted) for easier processing
    df_melted = df_cleaned.melt(id_vars=['Parameter'], var_name='Case ID', value_name='Value')
    df_melted['Sample'] = df_melted['Case ID'].str.extract(r'(\d+)').astype(int)
    # df_melted['Sample'] = df_melted['Case ID'].str[:-2]
    if mthd == 'TEST' and site == 'HUP':
        df_melted['User'] = df_melted['Case ID'].str[-3:].map(investigators)
    else:
        df_melted['User'] = df_melted['Case ID'].map(users_row)

    # Remove samples that exist for only one of the reviewers
    cases_counts = df_melted.groupby('Sample')['User'].nunique()
    single_reviewer_specimens = cases_counts[cases_counts < 2].index.tolist()
    if single_reviewer_specimens:
        print(f"Removing samples with only one reviewer: {single_reviewer_specimens}")
        df_melted = df_melted[~df_melted['Sample'].isin(single_reviewer_specimens)]

    # Pivot table to align values from both investigators
    robust_dup(df_melted, ['Sample', 'Parameter', 'User'])
    df_pivoted = df_melted.pivot(index=['Sample', 'Parameter'], columns='User', values='Value').reset_index()

    # arbitrator values do not need to go through additional arbitration
    if 'Arbitrator' in df_pivoted.columns:
        df_pivoted.drop(columns=['Arbitrator'], inplace=True)

    if len(df_pivoted) == 0:
        raise ValueError("Nothing to compare (no samples were reviewed by multiple reviewers yet)")

    # Rename columns for clarity
    df_pivoted.rename(columns={'User1': 'Value_1', 'User2': 'Value_2', 'User3': 'Value_3', 'User4': 'Value_4'}, inplace=True)

    # Convert numeric values where applicable

    for col in ['Value_1', 'Value_2', 'Value_3', 'Value_4']:  # might need to change script to fit number of reviewers per site
        if col in df_pivoted.columns:
            df_pivoted[col] = pd.to_numeric(df_pivoted[col], errors='ignore')

    df_pivoted.to_excel(rf'results\df_pivoted_{site}_{mthd}.xlsx', index=False)
    print(df_pivoted)



