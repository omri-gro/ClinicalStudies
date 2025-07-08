import pandas as pd
from test_comp import sen_spe, create_conf_mtrx
import os
from scipy.stats import bootstrap


if __name__ == "__main__":

    
    norm_ranges = {
        'Basophil': [0, 2],
        'Blast': [0, 0.5],
        'Eosinophil': [0, 5],
        'Immature': [0, 1],
        'Lymphocyte': [9, 40],
        'Metamyelocyte': [0, 1],
        'Monocyte': [2, 10],
        'Myelocyte': [0, 1],
        'Neut': [47, 80],
        'NRBC': [0, 1],
        'Plasma': [0, 0],
        'PLT': [150, 400],
        'Promyelocyte': [0, 1],
        'Total Lym': [9, 45],
        'Var Lym': [0, 12]
        }

    ref_var = 'Ref'
    test_var = 'Test'

    raw_data_name = r"PBS_eval_2024_WBC_PLT"
    raw_df = pd.read_csv(f'{raw_data_name}.csv')

    sites = list(raw_df['Site'].unique())
    sites.append('All')

    final_rows = []

    for site in sites:
        # Use only the specified site's results
        if site != 'All':
            site_df = raw_df.query("Site == @site")
        else:
            site_df = raw_df

        for cls, rng in norm_ranges.items():
            try:
                cls_ref = site_df[f'{cls} {ref_var}']  # name of column, in the format of 'LGL Ref'
                cls_test = site_df[f'{cls} {test_var}']

                pred_vals = sen_spe(cls_ref, cls_test, rng)
                out_row = {'site': site, 'class': cls} | pred_vals
                print(out_row)
                final_rows.append(out_row)

            except ValueError or KeyError as err:
                print(f'Error during predictive values calculation: {err}')
                continue

    pv_df = pd.DataFrame(final_rows)
    pv_df.to_excel(fr"{raw_data_name}_pred_vals.xlsx", index=False)

"""
    # RBC contingency tables
    ref_var = 'ref'
    test_var = 'test'

    raw_data_name = r"PBS_eval_2024_WBC_PLT"
    raw_df = pd.read_csv(f'{raw_data_name}.csv')
    # Define the project root dynamically
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    events_json_path = os.path.join(project_root, rf"data/sandbox/{uuid}.json")
"""


