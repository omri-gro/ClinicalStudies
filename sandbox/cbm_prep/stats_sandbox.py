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


def bootstrap_ci(metric_func, y_true, y_pred, n_boot=1000, alpha=0.05, random_state=None):
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

    Returns
    -------
    (lower_ci, upper_ci)
    """
    rng = np.random.default_rng(random_state)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    n = len(y_true)
    boot_vals = np.empty(n_boot)

    for i in range(n_boot):
        idx = rng.integers(0, n, n)  # resample with replacement
        boot_vals[i] = metric_func(y_true[idx], y_pred[idx])

    lower = np.nanpercentile(boot_vals, 100 * (alpha / 2))
    upper = np.nanpercentile(boot_vals, 100 * (1 - alpha / 2))

    return lower, upper


def binary_classification_metrics_bootstrap(
        y_true, y_pred, n_boot=1000, alpha=0.05, random_state=None
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
            n_boot=n_boot, alpha=alpha, random_state=random_state
        )
        results[key] = {
            "value": val,
            "ci": (lwr, upr)
        }

    for key in ["tp", "fn", "tn", "fp"]:
        results[key] = base[key]

    return results
