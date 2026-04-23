import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    test_arm = 'CBM'
    ref_arm = 'manual'

    exprt_long = True

    min_inv = 2  # False or number  currently does not seem to make much of a difference
    aftr_2nd_ssn = False

    scnd_ssn_str = '_aftr2ndssn' if aftr_2nd_ssn else ''
    save_name = f'mnl_mininv{min_inv}{scnd_ssn_str}'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    investigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1', 'Ana Catarina Silva': 'Rev1',
                        'Sarah Pereira Rodrigues': 'Rev1', 'Maria Buen Viana De Perio': 'Rev1',
                        'Christine Lavoie': 'Rev1', 'Ebikebuna Rufus': 'Rev1', 'Donald': 'Rev1',
                        'Sladana': 'Rev2', 'Deborah Swearingen': 'Rev2', 'Tony Omigie': 'Rev2',
                        'Joy Arthur': 'Rev2', 'Tiffany I Highsmith': 'Rev2', 'Tiffany I. Highsmith': 'Rev2',
                        'Harsha Hirani': 'Rev2', 'Harsha HIrani': 'Rev2',
                        'YAEL ASYEGH': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Yael S': 'Rev2', 'Yael Sayegh': 'Rev2',
                        'Yael S ': 'Rev2',
                        'Christopher Wright': 'Rev2', 'Thu Tran': 'Rev2',
                        'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}

    df_srcs_list = []
    for site in sites:
        raw_df = raw_to_df(f'{site}_manual.csv', site, ref_arm, dir=r'raw/cbm_method_comparison')
        df = stnd_names(raw_df, metadata.alias_map)
        df = calc_diff(df, metadata, additional_cells="WBC-like")
        df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
        df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp = MethodComparator(all_dfs)

    if aftr_2nd_ssn:
        rmv_file = 'flt_lists/pre_2nd_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    # calculate mean investigator
    df = methd_comp.df

    df = df.query("Value!='--------'")

    df['Investigator'] = df['Investigator'].map(investigators_map)
    df = add_mean_investigator(df, mthd=ref_arm, min_inv=min_inv)

    binary_vars = metadata.variable_groups["PLT morphology"] + metadata.variable_groups["RBC arrangement"] + metadata.variable_groups["WBC morphology"]
    raw_grade_cond=lambda d: (
        d["Method"].isin([ref_arm])
        & d["Variable"].isin(binary_vars)
    )

    df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value nor grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)


    cbm_file_name = '6sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, test_arm, metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'

    all_dfs = pd.concat([df, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    # cases to always remove - waiting for arbitration, horrible scans, etc.
    rmv_file = 'flt_lists/low_quality.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    rmv_file = 'flt_lists/for_arbitration.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)


    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups[
        'PLT morphology']
    grades_to_print = grades_to_test + ['ScanID']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also

    # keep only mean investigator
    methd_comp = methd_comp.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                        inplace=False)

    if exprt_long:
        include_in_export = vals_to_print + grades_to_print
        df_long = methd_comp.df.query(f"Variable in @include_in_export")[['SampleID', 'Site', 'Investigator', 'Variable', 'Value', 'Grade', 'Positive']]
        write_df_to_file(df_long, rf'comp_tables/{save_name}_long.csv')

    methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test)
    methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test, site_filters=sites)
    methd_comp.save_results(rf'results/mnl/{save_name}_reg.csv')

