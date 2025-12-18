import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle
from pipelines import medium_pipe

if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = r'config.yaml'
    metadata = MetadataBundle(meta_path)

    graded_vars = metadata.variable_groups['grade']
    rbc_vars = metadata.variable_groups['RBC morphology']
    rbc_grade_vars = list(set(graded_vars) & set(rbc_vars))


    sites = ['BWH', 'CPG', 'LMU', 'TASMC']
    srcs = {(site, 'OMR'): f'{site}_OMR.csv' for site in sites}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')

    df = methd_comp.df
    df_grades = df[df['Variable'].isin(rbc_grade_vars)]


    # variables that are not graded at all (or having no positives) in a specific sites to not be included in the analysis
    mask = (
        df_grades.groupby(['Site', 'Variable'])['Grade']
        .transform('nunique')
        .ne(1)
    )
    df_graded = df_grades[mask]
    df_graded.loc[:, "Value"] = df_graded["Grade"]

    cbm_file_name = '5sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')

    all_dfs = pd.concat([df_graded, cbm_df])
    methd_comp = MethodComparator(all_dfs)

    out_path = r'comp_tables/omr_cbm_rbc.csv'
    methd_comp.export_comparison_matrix(
        out_path=out_path,
        row_identifiers=["Site", "SampleID"],
        comparison_dims=("Variable", "Method"),
        needed_vals=rbc_grade_vars,
        needed_grades=(['ScanID'])
    )


