import pandas as pd
import numpy as np
import scipy.stats as stats
# from scipy.sparse import coo_matrix
from joblib import Parallel, delayed
from dataclasses import dataclass


@dataclass
class RegressionResult:
    slope: float
    intercept: float
    slope_ci: tuple[float, float]
    intercept_ci: tuple[float, float]
    r_squared: float
    method: str
    n_samples: int
    residuals: np.ndarray


"""
Functions for performing statistical analysis for arrays/Series
"""



