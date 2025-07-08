import numpy as np
import pandas as pd
from test_comp import regression_comp, force_num_cols, plot_ver_reg, plot_ver_reg_dbl
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    inter_bars = True
    # clv_morphs = ['Acanthocytes', 'Basophilic stippling', 'Bite cells', 'Blister cells', 'Burr Cells', 'Elliptocytes',
    #               'Helmet cells', 'Howell-Jolly', 'Ovalocytes', 'Pappenheimer', 'Parasites', 'Poikilocytosis',
    #               'Schistocytes', 'Sickle cells', 'Spherocytes', 'Stomatocytes', 'Target cells', 'Tear drop cells']
    # clv_morphs += ['AcanBurr', 'EllipOval']
    # clv_morphs = ['Acanthocytes', 'Burr Cells', 'AcanBurr',
    #               'Elliptocytes', 'Ovalocytes', 'EllipOval']
    # clv_morph_pos_exmp = ['Burr Cells', 'Elliptocytes', 'Ovalocytes', 'Pappenheimer', 'Poikilocytosis', 'Schistocytes',
    #                       'Target cells', 'Tear drop cells', 'AcanBurr', 'EllipOval']
    # clv_morph_1_exmp = ['Acanthocytes', 'Helmet cells', 'Parasites', 'Sickle cells', 'Stomatocytes']
    # clv_morphs = clv_morph_pos_exmp + clv_morph_1_exmp
    clv_morphs = ['Acanthocytes', 'Basophilic stippling', 'Bite cells', 'Blister cells', 'Echinocytes', 'Elliptocytes',
                   'Helmet cells', 'Howell-Jolly', 'Ovalocytes', 'Pappenheimer', 'Parasites', 'Poikilocytosis',
                   'Schistocytes', 'Sickle cells', 'Spherocytes', 'Stomatocytes', 'Target cells', 'Tear drop cells',
                   'AcanEchin', 'EllipOval', 'BiteHelSchi', 'BiteHelSchiPoik']

    cols_order = ["class", "slope_str", "intercept_str", "correlation_coefficient",
                  "regression method", "CI method", "N", "iterations", "slope", "slope_ci_bottom", "slope_ci_top",
                  "intercept", "intercept_ci_bottom", "intercept_ci_top"]

    ID_var = 'Case ID'
    ref_var = 'ClV mean'
    test_var = 'DSS'
    group_var = 'OMR'

    raw_data_name = "BWH_ClV_293"
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    data_path = os.path.join(project_root, 'data', 'sandbox', rf"{raw_data_name}.csv")
    src_df = pd.read_csv(data_path)
    for suf in ['ClV A', 'ClV B', 'ClV mean']:
        src_df[f'Burr Cells {suf}'] = src_df[f'Echinocytes {suf}']
        src_df[f'Basophilic stippling {suf}'] = src_df[f'Basophilic stippling {suf}']
        src_df[f'AcanBurr {suf}'] = src_df[f'Acanthocytes {suf}'] + src_df[f'Echinocytes {suf}']
        src_df[f'EllipOval {suf}'] = src_df[f'Elliptocytes {suf}'] + src_df[f'Ovalocytes {suf}']
    src_df[f'AcanBurr DSS'] = src_df[f'Acanthocytes DSS'] + src_df[f'Echinocytes DSS']
    src_df[f'EllipOval DSS'] = src_df[f'Elliptocytes DSS'] + src_df[f'Ovalocytes DSS']
    final_rows = []

    for cls in clv_morphs:
        cls_ref = f'{cls} {ref_var}'
        cls_test = f'{cls} {test_var}'
        cls_group = f'{cls} {group_var}'
        cls_df = force_num_cols(df=src_df, num_cols=[cls_ref, cls_test])
        ref_num = cls_df[cls_ref]
        test_num = cls_df[cls_test]
        try:
            reg = regression_comp(ref_num, test_num, lambda_=1, n_bootstrap=1000)
            # reg = regression_comp(ref_num, test_num, lambda_=1, n_bootstrap=200, reg_method="passing")
            reg.update({"class": cls})
            print(reg)
            final_rows.append(reg)
            #  plot_ver_reg(cls_df, cls_ref, cls_test, reg_ser=reg, cls_var=cls_group)
            if inter_bars:
                plot_ver_reg_dbl(cls_df, cls_ref, cls_test, meas1_2=[f'{cls} ClV A', f'{cls} ClV B'], reg_ser=reg)
            else:
                plot_ver_reg(cls_df, cls_ref, cls_test, reg_ser=reg, cls_var=cls_group, reg_ci=True)
        except ValueError or IndexError:
            print(f'Regression for {cls} not possible')
            continue

    regs_df = pd.DataFrame(final_rows, columns=cols_order)
    regs_df.to_excel(fr"{raw_data_name}_deming.xlsx", index=False)


