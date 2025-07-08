import pandas as pd
from test_comp import regression_comp, force_num_cols, plot_ver_reg
import os

if __name__ == "__main__":
    # all_classes = ['Band Neut', 'Seg Neut', 'Lymphocyte', 'Atypical Lymph', 'LGL', 'Aberrant Lymph', 'Monocyte',
    #                'Eosinophil', 'Basophil', 'Plasma', 'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Blast', 'NRBC',
    #                'Smudge Cells', 'Unclassified', 'PLT', 'Var Lym', 'Total Lym', 'Neut', 'Immature']
    # all_classes = ['Lymphocyte', 'Monocyte', 'Eosinophil', 'Basophil', 'Plasma Cell', 'Blast', 'NRBC',
    #                'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Segmented Neutrophil', 'Band Neutrophil',
    #                'Aberrant Lymphocyte', 'Atypical Lymphocyte',
    #                'Total Var Lym', 'Total Lym', 'Total Neutrophil', 'Total Imm']
    # all_classes = ['Total Lym', 'Total Neutrophil', 'Monocyte', 'Eosinophil', 'Basophil',
    #                'Plasma Cell', 'Blast', 'NRBC', 'Total Imm',
    #                'Aberrant Lymphocyte', 'Atypical Lymphocyte', 'Large Granular Lymphocyte',
    #                'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Segmented Neutrophil', 'Band Neutrophil']
    all_classes = ['Acanthocytes', 'Basophilic stippling', 'Bite cells', 'Blister cells', 'Echinocytes', 'Elliptocytes',
                   'Helmet cells', 'Howell-Jolly', 'Ovalocytes', 'Pappenheimer', 'Parasites', 'Poikilocytosis',
                   'Schistocytes', 'Sickle cells', 'Spherocytes', 'Stomatocytes', 'Target cells', 'Tear drop cells',
                   'AcanEchin', 'EllipOval', 'BiteHelSchi', 'BiteHelSchiPoik']
    cols_order = ["site", "class", "slope_str", "intercept_str", "correlation_coefficient",
                  "regression method", "CI method", "N", "iterations", "slope", "slope_ci_bottom", "slope_ci_top",
                  "intercept", "intercept_ci_bottom", "intercept_ci_top"]
    ID_var = 'Barcode'
    ID_var = 'Case ID'
    # ref_var = 'MNL A'
    # test_var = 'MNL B'
    ref_var = 'ClV mean'
    test_var = 'DSS'
    # ref_var = 'OMR'
    # test_var = 'DSS'

    # raw_data_name = r"BWH_manual"
    # raw_data_name = r"LMU_training"
    raw_data_name = r"BWH_ClV_293"



    # Define the project root dynamically
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    df_path = os.path.join(project_root, rf"data/sandbox/{raw_data_name}.csv")
    # df_path = os.path.join(os.getcwd(), rf"{raw_data_name}.csv")
    raw_df = pd.read_csv(df_path)

    # # Filter cases with unclass > 5%
    # raw_df = raw_df[raw_df['Unclassified WBC DSS'] < 5]
    # raw_df = raw_df[raw_df['Unclassified WBC MNL mean'] < 5]

    # Add aggregate values

    for mthd in ['ClV mean', 'ClV A', 'ClV B', 'DSS']:
        raw_df[f'EllipOval {mthd}'] = raw_df[f'Elliptocytes {mthd}'] + raw_df[f'Ovalocytes {mthd}']
        raw_df[f'AcanEchin {mthd}'] = raw_df[f'Acanthocytes {mthd}'] + raw_df[f'Echinocytes {mthd}']
        raw_df[f'BiteHelSchi {mthd}'] = raw_df[f'Bite cells {mthd}'] + raw_df[f'Helmet cells {mthd}'] + raw_df[f'Schistocytes {mthd}']
        raw_df[f'BiteHelSchiPoik {mthd}'] = raw_df[f'BiteHelSchi {mthd}'] + raw_df[f'Poikilocytosis {mthd}']
    raw_df.to_excel(fr"{raw_data_name}_with_agg_vars.xlsx", index=False)


    # sites = list(raw_df['Site'].unique())
    # sites.append('All')
    # sites = ['All', 'BWH', 'LMU']
    sites = ['All']

    final_rows = []

    for site in sites:
        # Use only the specified site's results
        if site != 'All':
            site_df = raw_df.query("Site == @site")
        else:
            site_df = raw_df

        for cls in all_classes[:]:
            # prepare df for analysis (only rows with numbers)
            cls_ref = f'{cls} {ref_var}'  # name of column, in the format of 'LGL Ref'
            cls_test = f'{cls} {test_var}'
            if "Site" in site_df.columns:
                cls_df = site_df[["Site", ID_var, cls_ref, cls_test]]
            else:
                cls_df = site_df[[ID_var, cls_ref, cls_test]]
            try:
                cls_df = force_num_cols(df=cls_df, num_cols=[cls_ref, cls_test])
                ref_num = cls_df[cls_ref]
                test_num = cls_df[cls_test]

                # try calculating regression and save results
                # deming
                reg = regression_comp(ref_num, test_num, lambda_=1, n_bootstrap=1000)
                reg.update({"class": cls, "site": site})
                print(reg)
                final_rows.append(reg)

                # passing
                # reg = regression_comp(ref_num, test_num, reg_method="passing", n_jobs=8, n_bootstrap=200)
                # reg.update({"class": cls, "site": site})
                # print(reg)
                # final_rows.append(reg)

            except ValueError or IndexError:
                print(f'Regression for {cls} in {site} site not possible')
                continue

            # plot regression and data points
            plot_ver_reg(cls_df, cls_ref, cls_test, cls_var="Site", reg_ser=reg)


    regs_df = pd.DataFrame(final_rows, columns=cols_order)
    regs_df.to_excel(fr"{raw_data_name}_deming.xlsx", index=False)
    # regs_df.to_excel(fr"{raw_data_name}_deming_inter.xlsx", index=False)
