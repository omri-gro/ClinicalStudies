import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe


if __name__ == "__main__":
    sites = ['CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    save_dir = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\results\sandbox_results'

    crf_ssn = 'post'  # 'all', 'pre' or 'post'


    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    invstigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1',
                        'Sarah Pereira Rodrigues': 'Rev1', 'Maria Buen Viana De Perio': 'Rev1',
                        'Sladana': 'Rev2', 'Deborah Swearingen': 'Rev2', 'Tony Omigie': 'Rev2',
                        'Joy Arthur': 'Rev2', 'Tiffany I Highsmith': 'Rev2', 'Tiffany I. Highsmith': 'Rev2',
                        'YAEL ASYEGH': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Yael S': 'Rev2', 'Yael Sayegh': 'Rev2',
                        'Yael S ': 'Rev2',
                        'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}

    df_srcs_list = []
    for site in sites:
        raw_df = raw_to_df(f'{site}_manual.csv', site, "manual", dir=r'raw/cbm_method_comparison')
        df = stnd_names(raw_df, metadata.alias_map)
        df = calc_diff(df, metadata, additional_cells="WBC-like")
        df = pivot_long(df, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'])
        df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp = MethodComparator(all_dfs)

    if crf_ssn == "post":
        rmv_file = 'flt_lists/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    df = methd_comp.df

    df = df.query("Value!='--------'")

    df['Investigator'] = df['Investigator'].map(invstigators_map)


    binary_vars = metadata.variable_groups["binary"]
    raw_grade_cond=lambda d: (
        d["Method"].isin(['manual'])
        & d["Variable"].isin(binary_vars)
    )

    df = add_grade_column(df, metadata, raw_grade_cond=raw_grade_cond)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade", "Positive"], how='all')  # drop when neither value or grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)


    cbm_file_name = '6sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    cbm_df['Investigator'] = 'CBM'

    all_dfs = pd.concat([df, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    rmv_file = 'flt_lists/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2',
                             variables='Myelocyte', dim_col='Investigator',
                             row_filters={'Method': 'manual'}, split_by='Site')
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.batch_compare(levels_a='Rev1', levels_b='Rev2', comp_func='binary',
                             variables='Pelger Cell', dim_col='Investigator',
                             row_filters={'Method': 'manual'}, split_by='Site')

    methd_comp.save_results(fr'{save_dir}/myelo_manual_inter.csv')
    methd_comp.save_results(fr'{save_dir}/myelo_manual_inter_bias.csv', 'Bias')
    methd_comp.save_results(fr'{save_dir}/pelger_manual_inter.csv', 'Binary')

    methd_comp.plot_all_regressions(fr'{save_dir}/myelo_manual_inter.pdf')
