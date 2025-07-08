import pandas as pd

if __name__ == "__main__":
    mthd = 'TEST'  # for test, remember to change name of 'Total' raw to 'Total nucleated cells'

    investigators = {'Wei Xie': 'User1', 'Todd Williams': 'User2'}
    typed_file = f'OHSU_CRF_{mthd}.csv'
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
    df_melted['User'] = df_melted['Case ID'].map(users_row)

    # Remove samples that exist for only one of the reviewers
    cases_counts = df_melted.groupby('Sample')['User'].nunique()
    single_reviewer_specimens = cases_counts[cases_counts < 2].index.tolist()
    if single_reviewer_specimens:
        print(f"Removing samples with only one reviewer: {single_reviewer_specimens}")
        df_melted = df_melted[~df_melted['Sample'].isin(single_reviewer_specimens)]

    # Pivot table to align values from both investigators
    df_pivoted = df_melted.pivot(index=['Sample', 'Parameter'], columns='User', values='Value').reset_index()

    # Rename columns for clarity
    df_pivoted.rename(columns={'User1': 'Value_1', 'User2': 'Value_2'}, inplace=True)

    # Convert numeric values where applicable
    for col in ['Value_1', 'Value_2']:
        df_pivoted[col] = pd.to_numeric(df_pivoted[col], errors='ignore')

    df_pivoted.to_excel(f'df_pivoted_{mthd}.xlsx', index=False)
    print(df_pivoted)



