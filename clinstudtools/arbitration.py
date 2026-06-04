import pandas as pd


def parse_arbitration_rules(rules_df: pd.DataFrame, metadata: "MetadataBundle") -> pd.DataFrame:
    """
    Converts a human-readable arbitration tracking dataframe into a strict long-format dataframe.
    Optionally looks for a 'Method' column to restrict overrides.
    """
    expanded_rules = []
    has_method = 'Method' in rules_df.columns

    for _, row in rules_df.iterrows():
        target_str = row.get('Target')
        if pd.isna(target_str):
            continue

        targets = [t.strip() for t in str(target_str).split(',')]

        for target in targets:
            vars_to_apply = metadata.variable_groups.get(target, [target])

            for v in vars_to_apply:
                rule = {
                    'SampleID': str(row['SampleID']).zfill(5),
                    'Site': row['Site'],
                    'Variable': v
                }
                # If a specific method was designated in the CSV, enforce it
                if has_method and pd.notna(row['Method']):
                    rule['Method'] = row['Method']

                expanded_rules.append(rule)

    if not expanded_rules:
        cols = ['SampleID', 'Site', 'Variable'] + (['Method'] if has_method else [])
        return pd.DataFrame(columns=cols)

    return pd.DataFrame(expanded_rules).drop_duplicates()


def apply_arbitration_override(df: pd.DataFrame, arb_df: pd.DataFrame, rules_df: pd.DataFrame,
                               metadata: "MetadataBundle") -> pd.DataFrame:
    """
    Selectively replaces 'Mean Investigator' values with Arbitrator values across all methods.
    Drops 'Mean Investigator' rows if arbitration is required but data is missing.
    """

    rules_keys = parse_arbitration_rules(rules_df, metadata)
    if rules_keys.empty:
        return df

    # Determine merge keys
    idx_cols = ['SampleID', 'Site', 'Variable', 'Investigator']
    if 'Method' in rules_keys.columns:
        idx_cols.append('Method')

    rules_keys['Investigator'] = 'Mean Investigator'
    arb_prep = arb_df.copy()
    arb_prep['Investigator'] = 'Mean Investigator'

    df_out = df.set_index(idx_cols)
    rules_idx = rules_keys.set_index(idx_cols).index
    arb_idx = arb_prep.set_index(idx_cols).index

    # Set Math: Fulfilled = Rules AND Arb. Pending = Rules MINUS Arb.
    fulfilled_idx = rules_idx.intersection(arb_idx)
    pending_idx = rules_idx.difference(arb_idx)

    # Drop pending cases (errors='ignore' prevents crashes if a rule doesn't exist in df)
    df_out = df_out.drop(index=pending_idx, errors='ignore')

    # Update fulfilled cases
    df_out.update(arb_prep.set_index(idx_cols).loc[fulfilled_idx, ['Value']])

    return df_out.reset_index()
