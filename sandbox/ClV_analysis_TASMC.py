import numpy as np
import pandas as pd
from test_comp import regression_comp, force_num_cols, plot_ver_reg, plot_ver_reg_dbl
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from sandbox import MetadataBundle


def rename_cols_base(df, rename_map, suffices):
    def rename_columns(col):
        # Check if column ends with a known suffix
        for suf in sorted(suffices, key=len, reverse=True):
            if col.endswith(" " + suf):  # suffix preceded by space
                var = col[:-(len(suf) + 1)]  # everything before suffix + space
                new_var = rename_map.get(var, var)
                return f"{new_var} {suf}"
        # If no suffix match, rename whole variable name if possible
        return rename_map.get(col, col)

    return df.rename(columns=rename_columns)

if __name__ == "__main__":
    inter_bars = True

    # clv_morphs = ['Polychromatic cells', 'Hypochromatic cells', 'Microcytes', 'Macrocytes', 'Target cells',
    #               'Schistocytes', 'Helmet cells', 'Sickle cells', 'Spherocytes', 'Elliptocytes', 'Ovalocytes',
    #               'Tear drop cells', 'Stomatocytes', 'Acanthocytes', 'Echinocytes', 'Bite cells', 'Blister cells',
    #               'Howell-Jolly', 'Pappenheimer', 'Basophilic stippling', 'Parasites',
    #               'Large platelets', 'Giant platelets', 'Agranular platelets']

    suffices = ['rev1', 'rev2', 'ClV mean', 'DSS']

    ID_var = 'Case ID'
    ref_var = 'ClV mean'
    test_var = 'DSS'

    meta_path = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\config.yaml'
    metadata = MetadataBundle(meta_path)
    rename_map = metadata.alias_map

    raw_data_name = "TASMC_rbc_side_by_side"
    pdf_path = fr'{raw_data_name}.pdf'
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    data_path = os.path.join(project_root, 'data', 'sandbox', rf"{raw_data_name}.csv")
    src_df = pd.read_csv(data_path)

    src_df = rename_cols_base(src_df, rename_map, suffices)

    clv_morphs = metadata.variable_groups['RBC morphology'] + metadata.variable_groups['PLT morphology']

    with PdfPages(pdf_path) as pdf:
        for cls in clv_morphs:
            try:
                cls_ref = f'{cls} {ref_var}'
                cls_test = f'{cls} {test_var}'
                cls_df = src_df

                fig, ax = plot_ver_reg_dbl(cls_df, cls_ref, cls_test, meas1_2=[f'{cls} rev1', f'{cls} rev2'])
                pdf.savefig(fig)
                plt.close()
            except Exception as e:
                print(f'Plot for {cls} not possible: {e}')
                continue

