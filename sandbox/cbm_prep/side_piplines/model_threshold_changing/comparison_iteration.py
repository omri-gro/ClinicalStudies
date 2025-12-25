import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import mean_manual_pipe, medium_pipe

def comp_mnl(cbm_file_name, mnl_dir, vars_to_test, crf_ssn='post', diff500=False):
    # can break this into multiple functions later
    sites = ['CPG', 'LMU', 'SYN', 'TASMC']
    meta_path = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\config.yaml'
    flt_lists_dir = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\flt_lists'

    metadata = MetadataBundle(meta_path)

    invstigators_map = {'Alina': 'Rev1', 'Aubrey B Charlton': 'Rev1', 'Thomas Muddiman': 'Rev1',
                        'Sarah Pereira Rodrigues': 'Rev1',
                        'Sladana': 'Rev2', 'Deborah Swearingen': 'Rev2', 'Tony Omigie': 'Rev2',
                        'YAEL ASYEGH': 'Rev2', 'YAEL SAYEGH': 'Rev2', 'Yael S': 'Rev2', 'Yael Sayegh': 'Rev2',
                        'CBM': 'CBM', 'Mean Investigator': 'Mean Investigator'}

    df_srcs_list = []
    for site in sites:
        df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=mnl_dir, only_mean=False)
        df_srcs_list.append(df)

    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata,
                         dir=r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\side_piplines\model_threshold_changing\5sites_csvs')
    cbm_df['Investigator'] = 'CBM'
    df_srcs_list.append(cbm_df)

    all_dfs = pd.concat(df_srcs_list)
    all_dfs['Investigator'] = all_dfs['Investigator'].map(invstigators_map)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'",
                                                inplace=False)

    rmv_file = f'{flt_lists_dir}/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    if diff500:
        diff500_file = f'{flt_lists_dir}/500_WBC_mnl_cases.csv'
        diff500_df = read_to_df(diff500_file, file_dir=os.getcwd())
        diff500_comp_df = methd_comp.filter_by_df(diff500_df, include_rows=True).df.copy()

        cbm_df = methd_comp.only_when_cond("Method=='CBM'").df.copy()
        diff500_dfs = pd.concat([diff500_comp_df, cbm_df])
        methd_comp = MethodComparator(diff500_dfs)

    if crf_ssn == "post":
        rmv_file = f'{flt_lists_dir}/pre_session_reviews.csv'
        rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(rmv_df)

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=sites)
    return methd_comp
