import os
import sys
import pandas as pd

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df, add_mean_investigator, add_grade_column, add_pos_column, create_derived_variables_long
from pipelines import medium_pipe, clv_pipe


if __name__ == "__main__":

    # calculating anisocytosis,
