import pandas as pd
import os
from objects import MethodComparator
from sandbox import MetadataBundle
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    suffix = ''
    meta_path = r'../config.yaml'
    site = 'LMU'
    save_name = f'LMU_inter'
    analysis_name = "cbm_method_comparison"
    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'../raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    df_srcs_list = []

    # import and pre-process each method separately
    df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=r'../raw/cbm_method_comparison', only_mean=False)
    df['Method'] = df['Investigator']
    df = df[df['Investigator'].isin(['Alina', 'Sladana'])]

    methd_comp = MethodComparator(df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']

    methd_comp.batch_fit(['Alina'], ['Sladana'], vars_to_test)
    methd_comp.save_results(rf'results/{save_name}_manual_regs.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_manual_reg.pdf')







