import os
import sys
import pandas as pd, matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr, kendalltau
from statsmodels.miscmodels.ordinal_model import OrderedModel

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import medium_pipe
from itertools import *




if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = r'config.yaml'
    save_name = 'clv_to_cbm'

    metadata = MetadataBundle(meta_path)

    rmv_file = 'slides_to_remove_long.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())


