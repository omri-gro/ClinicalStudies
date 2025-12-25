import os
import pandas as pd
import numpy as np

from aggregation import results_df_from_jsons_dir_threshold_change
from comparison_iteration import comp_mnl

if __name__ == "__main__":
    orig_cell = 'Aberrant Lymphocyte'
    dst_cell = "Lymphocyte"

    vars_to_test = [orig_cell, dst_cell]
    print_also = ['Unclassified WBC', "Total WBC"]
    vals_to_print = vars_to_test + print_also

    cur_thrs = 0.55
    thresholds = np.arange(cur_thrs, 0.95, 0.05)

    jsons_path = r'result_dicts'
    dss_csvs_path = r'5sites_csvs'

    # for appropriate saving names
    cell_str = orig_cell.replace(' ', '_')

    scan_map = pd.read_csv('scans_in_study.csv')

    # for thrs in thresholds:
    #     thresh_str = f'{thrs:.2f}'.replace('.', '_')
    #     csv_file_name = rf'{dss_csvs_path}/{cell_str}_thresh{thresh_str}_cbm.csv'
    #     cbm_df = pd.read_csv(csv_file_name)
    #     cbm_df.rename(columns={'Scan UUID': 'ScanID'}, inplace=True)
    #     merged_df = pd.merge(scan_map, cbm_df, on='ScanID', how='inner')
    #     merged_df.to_csv(csv_file_name, index=False)




    for thrs in thresholds:
        thresh_str = f'{thrs:.2f}'.replace('.', '_')
        csv_file_name = f'{cell_str}_thresh{thresh_str}_cbm.csv'
        results_df = results_df_from_jsons_dir_threshold_change(jsons_path, orig_cell=orig_cell, dst_cell=dst_cell, thresh=thrs)
        results_df = pd.merge(scan_map, results_df, on='ScanID', how='inner')

        results_df.to_csv(rf'{dss_csvs_path}/{cell_str}_thresh{thresh_str}_cbm.csv')

        thresh_mthd_comp = comp_mnl(csv_file_name,
                                    r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\raw\cbm_method_comparison',
                                    [orig_cell, dst_cell])

        save_name = f'{cell_str}_thresh{thresh_str}_to_mnl'
        thresh_mthd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}.csv',
                                                  row_identifiers=["Site", "SampleID"],
                                                  comparison_dims=("Variable", "Method"),
                                                  needed_vals=vals_to_print,
                                                  needed_grades=['ScanID'])

        thresh_mthd_comp.save_results(rf'results/{save_name}.csv')
        thresh_mthd_comp.plot_all_regressions(rf'results/{save_name}.pdf')










