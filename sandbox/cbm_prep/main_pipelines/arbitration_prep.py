import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe


def keep_pairs_only(df):
    """
    Filters the dataframe to keep only SampleID-Site-Variable
    combinations that have exactly two unique investigators.
    """
    # Define the columns that make up a unique task combination
    group_cols = ['SampleID', 'Site', 'Variable']

    # Create a boolean mask: True if the group has exactly 2 unique investigators
    mask = df.groupby(group_cols)['Investigator'].transform('nunique') == 2

    # Apply the mask to the dataframe and return a new copy
    return df[mask].copy()


def rev_mapping(df):
    # can be added as future method

    grouped_investigators = df.groupby(['SampleID', 'Site', 'Variable'])['Investigator'].unique()
    reviewer_map = {}

    for investigators in grouped_investigators:
        if len(investigators) == 2:
            inv1, inv2 = investigators[0], investigators[1]

            # If neither investigator has been assigned a role yet
            if inv1 not in reviewer_map and inv2 not in reviewer_map:
                reviewer_map[inv1] = 'Rev1'
                reviewer_map[inv2] = 'Rev2'

            # If the first investigator already has a role, assign the opposite to the second
            elif inv1 in reviewer_map and inv2 not in reviewer_map:
                reviewer_map[inv2] = 'Rev2' if reviewer_map[inv1] == 'Rev1' else 'Rev1'

            # If the second investigator already has a role, assign the opposite to the first
            elif inv2 in reviewer_map and inv1 not in reviewer_map:
                reviewer_map[inv1] = 'Rev1' if reviewer_map[inv2] == 'Rev2' else 'Rev2'

            # If both are already in the map, we leave them as-is to preserve global consistency.

        elif len(investigators) == 1:
            inv1 = investigators[0]
            # Assign a default role if they work alone and haven't been seen yet
            if inv1 not in reviewer_map:
                reviewer_map[inv1] = 'Rev1'

    print(reviewer_map)
    df['Reviewer'] = df['Investigator'].map(reviewer_map)
    return df



if __name__ == "__main__":
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    ref_arm = 'manual'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    if ref_arm=='manual':
        df_srcs_list = []
        for site in sites:
            raw_df = raw_to_df(f'{site}_manual.csv', site, ref_arm, dir=r'raw/cbm_method_comparison')
            df = stnd_names(raw_df, metadata.alias_map)
            df = calc_diff(df, metadata, additional_cells="WBC-like")
            df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
            df_srcs_list.append(df)

        all_dfs = pd.concat(df_srcs_list)
        methd_comp = MethodComparator(all_dfs)

        # calculate mean investigator
        df = methd_comp.df

        df = df.query("Value!='--------'")

        binary_vars = metadata.variable_groups["PLT morphology"] + metadata.variable_groups["RBC arrangement"] + \
                      metadata.variable_groups["WBC morphology"]
        raw_grade_cond = lambda d: (
                d["Method"].isin([ref_arm])
                & d["Variable"].isin(binary_vars)
        )

        df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
        df = add_pos_column(df, metadata)
        df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value nor grade in row
        df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
        df = create_derived_variables_long(df, metadata)
        methd_comp = MethodComparator(df)


        # cases to always remove - waiting for arbitration, horrible scans, etc.
        rmv_file = 'flt_lists/low_quality.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

        rmv_file = 'flt_lists/for_arbitration.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

        vars_to_test = metadata.variable_groups['WBC&PLT compare']
        grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
            'PLT morphology'] + metadata.variable_groups['RBC arrangement']
        grades_to_print = grades_to_test + ['ScanID']
        morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
        print_also = ['Unclassified WBC', "Total WBC"]
        vals_to_print = vars_to_test + print_also
        morph_vals_to_print = morph_vals_to_test + print_also

    include_in_export = metadata.variable_groups['WBC diff'] + print_also

    df = methd_comp.df
    df = rev_mapping(df)
    df = keep_pairs_only(df)
    # df_long = methd_comp.df.query(f"Variable in @include_in_export")[
    #     ['SampleID', 'Site', 'Investigator', 'Reviewer', 'Variable', 'Value']]
    # write_df_to_file(df_long, rf'comp_tables/{ref_arm}_long.csv')

    df_long_grades = methd_comp.df.query(f"Variable in @grades_to_test")[
        ['SampleID', 'Site', 'Investigator', 'Reviewer', 'Variable', 'Value', "Grade", "Positive"]]
    write_df_to_file(df_long_grades, rf'comp_tables/{ref_arm}_grades_long.csv')


