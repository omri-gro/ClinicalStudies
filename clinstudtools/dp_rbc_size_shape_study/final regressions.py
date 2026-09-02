import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *


if __name__ == "__main__":
    variables = ['Spherocytes', 'Schistocytes', 'Macrocytes', 'Microcytes', 'Macroplatelets']

    df = read_to_df(f'DP_final_df.csv', file_dir=r'raw')
    methd_comp = MethodComparator(df)
    methd_comp.batch_fit('Investigators', 'CBM', variables)
    methd_comp.save_results(rf'results/dp_reg.csv')
    methd_comp.plot_all_regressions(f'results/dp_reg.pdf')




