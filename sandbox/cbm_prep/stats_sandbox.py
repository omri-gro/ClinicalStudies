import pandas as pd
import numpy as np
import scipy.stats as stats
# from scipy.sparse import coo_matrix
from joblib import Parallel, delayed
from dataclasses import dataclass
from sklearn.utils import resample
from statsmodels.stats.proportion import proportion_confint


# Global cache dictionary to store previously calculated CIs
_ci_cache = {}

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


def binary_classification_metrics(y_true, y_pred):
    """
    Calculate sensitivity, specificity, and overall agreement.
    Parameters
    ----------
    y_true : array-like of bool
    y_pred : array-like of bool

    Returns
    -------
    dict with keys: sensitivity, specificity, agreement
    """
    y_true = np.asarray(y_true).astype(bool)
    y_pred = np.asarray(y_pred).astype(bool)

    tp = np.sum((y_true == True) & (y_pred == True))
    fn = np.sum((y_true == True) & (y_pred == False))
    tn = np.sum((y_true == False) & (y_pred == False))
    fp = np.sum((y_true == False) & (y_pred == True))

    sensitivity = tp / (tp + fn) if tp + fn > 0 else np.nan
    specificity = tn / (tn + fp) if tn + fp > 0 else np.nan
    agreement = (tp + tn) / len(y_true)

    return {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "agreement": agreement,
        "tp": tp,
        "fn": fn,
        "tn": tn,
        "fp": fp
    }


def bootstrap_ci(metric_func, y_true, y_pred, n_boot=1000, alpha=0.05, random_state=None,
                 stratify=None, data_df=None):
    """
    Generic bootstrap CI calculator for a metric based on paired y_true/y_pred.

    Parameters
    ----------
    metric_func : callable
        Should accept (y_true, y_pred) and return a float.
    y_true, y_pred : array-like
    n_boot : int
        Number of bootstrap samples.
    alpha : float
        Significance level (default 0.05 gives 95% CI).
    random_state : int or None
    If data_df and stratify are provided, it performs stratified resampling.

    Returns
    -------
    (lower_ci, upper_ci)
    """
    rng = np.random.default_rng(random_state)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    n = len(y_true)
    boot_vals = np.empty(n_boot)

    # Check if we have the full dataframe for stratification
    if stratify is not None and data_df is not None:
        stratify_data = data_df[stratify]
    else:
        stratify_data = None

    for i in range(n_boot):
        if stratify_data is not None:
            # Resample indices based on stratification
            boot_idx = resample(np.arange(len(y_true)), replace=True, stratify=stratify_data)
            t_boot, p_boot = y_true[boot_idx], y_pred[boot_idx]
        else:
            # Legacy simple resample
            t_boot, p_boot = resample(y_true, y_pred, replace=True)

        boot_vals[i] = metric_func(t_boot, p_boot)

    lower = np.nanpercentile(boot_vals, 100 * (alpha / 2))
    upper = np.nanpercentile(boot_vals, 100 * (1 - alpha / 2))

    return lower, upper


def binary_classification_metrics_bootstrap(
        y_true, y_pred, n_boot=1000, alpha=0.05, random_state=None,
        stratify_cols=None, data_df=None
):
    """
    Calculates sensitivity, specificity, agreement + bootstrap CIs.

    Returns
    -------
    dict with metric values and confidence intervals
        {"sensitivity": {"value": 0.8, "ci": (0.7, 0.9)},
        "specificity": ...,
        "tp": 83, ...}
    """
    base = binary_classification_metrics(y_true, y_pred)

    def make_metric_func(key):
        return lambda t, p: binary_classification_metrics(t, p)[key]

    results = {}
    for key in ["sensitivity", "specificity", "agreement"]:
        val = base[key]
        lwr, upr = bootstrap_ci(
            make_metric_func(key), y_true, y_pred,
            n_boot=n_boot, alpha=alpha, random_state=random_state,
            stratify=stratify_cols, data_df=data_df
        )
        results[key] = {
            "value": val,
            "ci": (lwr, upr)
        }

    for key in ["tp", "fn", "tn", "fp"]:
        results[key] = base[key]

    return results


# functions related to "Borderline / Equivalence Zone" Analysis (statistical artifacts on sensitivity/specificity based on Poisson noise - added due to FDA request
def check_borderline_equivalence(ref_pct, test_pct, total_cells=500, alpha=0.01):
    """
    Determines if the test percentage falls within the exact Clopper-Pearson CI of the reference.
    Useful for defining Equivocal Zones in cell counting studies.
    Uses memoization to avoid redundant, slow Beta distribution calculations.
    """
    if pd.isna(ref_pct) or pd.isna(test_pct):
        return False

    ref_count = round((ref_pct / 100) * total_cells)
    lower_prop, upper_prop = proportion_confint(ref_count, total_cells, alpha=alpha, method='beta')

    return (test_pct >= lower_prop * 100) and (test_pct <= upper_prop * 100)


def get_equivocal_zone_masks(y_true_cont, y_pred_cont, cutoff, total_cells=500, alpha=0.01):
    """
    Binarizes continuous arrays and identifies which discordant pairs fall within the statistical noise.
    Returns:
        y_true_bin: Strict binary reference array
        y_pred_bin: Strict binary test array
        exclude_mask: Boolean mask indicating which cases are equivocal and should be dropped
    """
    y_true_cont = np.asarray(y_true_cont)
    y_pred_cont = np.asarray(y_pred_cont)

    y_true_bin = (y_true_cont >= cutoff).astype(int)
    y_pred_bin = (y_pred_cont >= cutoff).astype(int)

    exclude_mask = np.zeros(len(y_true_cont), dtype=bool)

    for i in range(len(y_true_cont)):
        if y_true_bin[i] != y_pred_bin[i]:
            if check_borderline_equivalence(y_true_cont[i], y_pred_cont[i], total_cells, alpha):
                exclude_mask[i] = True

    return y_true_bin, y_pred_bin, exclude_mask


# might not be used anywhere, in which case can be deleted as replaced by get_equivocal_zone_masks
def calculate_adjusted_contingency(df, ref_col, test_col, cutoff, total_cells=500):
    """
    Takes a dataframe of continuous percentages and evaluates adjusted Sen/Spe.
    """
    df = df.copy()

    # Strict Binary Classification
    df['ref_binary'] = (df[ref_col] >= cutoff).astype(int)
    df['test_binary'] = (df[test_col] >= cutoff).astype(int)

    # Determine strict contingency status
    conditions = [
        (df['ref_binary'] == 1) & (df['test_binary'] == 1),
        (df['ref_binary'] == 0) & (df['test_binary'] == 0),
        (df['ref_binary'] == 0) & (df['test_binary'] == 1),
        (df['ref_binary'] == 1) & (df['test_binary'] == 0)
    ]
    choices = ['TP', 'TN', 'FP', 'FN']
    df['strict_status'] = np.select(conditions, choices, default='Unknown')

    # Identify Borderline Equivalence
    df['is_borderline'] = df.apply(
        lambda row: check_borderline_equivalence(row[ref_col], row[test_col], total_cells)[0],
        axis=1
    )

    # Adjust Status: If it's an FP or FN but within noise, it becomes Concordant (TN or TP)
    def adjust_status(row):
        if row['strict_status'] == 'FP' and row['is_borderline']:
            return 'TN (Adjusted)'
        elif row['strict_status'] == 'FN' and row['is_borderline']:
            return 'TP (Adjusted)'
        else:
            return row['strict_status']

    df['adjusted_status'] = df.apply(adjust_status, axis=1)

    # Map back to binary for adjusted bootstrap calculations
    df['adjusted_test_binary'] = df.apply(
        lambda row: row['ref_binary'] if 'Adjusted' in row['adjusted_status'] else row['test_binary'],
        axis=1
    )

    return df



