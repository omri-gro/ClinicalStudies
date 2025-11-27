import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from itertools import *


if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = r'config.yaml'
    save_name = 'site_combinations_omr'

    metadata = MetadataBundle(meta_path)

    # build list of files to read from
    sites = ['BWH', 'CPG', 'LMU', 'SYN', 'TASMC']
    mthds = ['OMR', 'CBM']
    srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')

    methd_comp.df["SampleID"] = methd_comp.df["Site"] + methd_comp.df["SampleID"]


    vars_to_test = metadata.variable_groups['WBC&PLT compare']

    all_combinations = []
    for r in range(2, len(sites)):
        combinations_of_length_r = combinations(sites, r)
        all_combinations.extend(list(combinations_of_length_r))

    comb_methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    for comb in all_combinations:
        df_to_filter = methd_comp.df
        comb_methd_comp.df = df_to_filter.query("Site in @comb")
        comb_name = ','.join(comb)
        comb_methd_comp.df["Site"] = comb_name
        comb_methd_comp.batch_fit('OMR', 'CBM', vars_to_test, site_filters=comb_name)

    comb_methd_comp.save_results(rf'results/{save_name}.csv')


    print(all_combinations)


