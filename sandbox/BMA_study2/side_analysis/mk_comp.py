import pandas as pd
import os
import sys
os.chdir('..')
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator, RegressionResult
from sandbox import MetadataBundle, raw_bma_to_df, _ensure_list
from sandbox import *
from regressions import regression_comp


df = pd.read_csv(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\BMA_study2\raw\mk_raw.csv')
reg_dict = regression_comp(df['Mean Digital'], df['Mean DSS'], reg_method='deming')
mk_reg = RegressionResult.from_dict(reg_dict)
print(reg_dict)
biases = calc_bias_at_points(mk_reg, [2, 10])
print(biases)



