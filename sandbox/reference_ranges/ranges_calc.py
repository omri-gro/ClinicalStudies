import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import binom
import statsmodels.api as sm


# ==========================================
# Helper Function for Zero-Heavy Parameters
# ==========================================
def calculate_onesided_upper_limit(data, ref_level=0.95, ci_level=0.90):
    """
    Calculates a one-sided upper reference limit for zero-heavy distributions.
    As per CLSI C28-A3c, the reference interval becomes 0 to X.
    """
    data_sorted = np.sort(data)
    n = len(data_sorted)

    # For a one-sided 95% limit, we typically look at the 95th percentile
    # (or 97.5th depending on strict clinical definitions). Here we use the upper bound
    # of the standard central interval.
    upper_p = ref_level + ((1.0 - ref_level) / 2.0)
    upper_limit = np.percentile(data_sorted, upper_p * 100)

    alpha = 1.0 - ci_level
    ul_rank_low = int(binom.ppf(alpha / 2, n, upper_p))
    ul_rank_high = int(binom.ppf(1 - (alpha / 2), n, upper_p))

    ul_idx_low = max(0, ul_rank_low - 1)
    ul_idx_high = min(n - 1, ul_rank_high)

    # Check if we actually have enough data to form a CI, otherwise return N/A for CIs
    try:
        ul_ci = (data_sorted[ul_idx_low], data_sorted[ul_idx_high])
    except IndexError:
        ul_ci = (np.nan, np.nan)

    return {
        "lower_limit": 0.0,
        "upper_limit": upper_limit,
        "lower_limit_ci": (0.0, 0.0),
        "upper_limit_ci": ul_ci,
        "method": "One-Sided Non-Parametric"
    }


# ==========================================
# 1. Non-Parametric Method (n >= 120)
# ==========================================
def calculate_nonparametric_ri(data, ref_level=0.95, ci_level=0.90):
    data_sorted = np.sort(data)
    n = len(data_sorted)

    lower_p = (1.0 - ref_level) / 2.0
    upper_p = 1.0 - lower_p

    lower_limit = np.percentile(data_sorted, lower_p * 100)
    upper_limit = np.percentile(data_sorted, upper_p * 100)

    alpha = 1.0 - ci_level

    ll_rank_low = int(binom.ppf(alpha / 2, n, lower_p))
    ll_rank_high = int(binom.ppf(1 - (alpha / 2), n, lower_p))
    ul_rank_low = int(binom.ppf(alpha / 2, n, upper_p))
    ul_rank_high = int(binom.ppf(1 - (alpha / 2), n, upper_p))

    ll_idx_low = max(0, ll_rank_low - 1)
    ll_idx_high = min(n - 1, ll_rank_high)
    ul_idx_low = max(0, ul_rank_low - 1)
    ul_idx_high = min(n - 1, ul_rank_high)

    return {
        "lower_limit": lower_limit,
        "upper_limit": upper_limit,
        "lower_limit_ci": (data_sorted[ll_idx_low], data_sorted[ll_idx_high]),
        "upper_limit_ci": (data_sorted[ul_idx_low], data_sorted[ul_idx_high]),
        "method": "Non-Parametric"
    }


# ==========================================
# 2. Robust Method (n < 120)
# ==========================================
def calculate_robust_reference_interval(data, ref_level=0.95, ci_level=0.90, n_bootstraps=1000):
    n = len(data)

    # Zero-heavy check: If MAD is 0, robust will fail with division by zero
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    if mad == 0 or np.any(data <= 0):
        return calculate_onesided_upper_limit(data, ref_level, ci_level)

    def get_robust_limits(sample_data):
        transformed, lmbda = stats.boxcox(sample_data)
        try:
            robust_mean, robust_std = sm.robust.scale.huber(transformed)
        except Exception:
            return None, None, None

        alpha_ref = 1.0 - ref_level
        z_score = stats.norm.ppf(1.0 - alpha_ref / 2.0)

        lower_trans = robust_mean - (z_score * robust_std)
        upper_trans = robust_mean + (z_score * robust_std)

        def inverse_boxcox(y, lmbda):
            if lmbda == 0:
                return np.exp(y)
            return (y * lmbda + 1.0) ** (1.0 / lmbda)

        return inverse_boxcox(lower_trans, lmbda), inverse_boxcox(upper_trans, lmbda), lmbda

    lower_limit, upper_limit, _ = get_robust_limits(data)

    # If the primary calculation fails, fallback to one-sided limit
    if lower_limit is None or np.isnan(lower_limit):
        return calculate_onesided_upper_limit(data, ref_level, ci_level)

    boot_lower_limits, boot_upper_limits = [], []
    for _ in range(n_bootstraps):
        boot_sample = np.random.choice(data, size=n, replace=True)
        ll, ul, _ = get_robust_limits(boot_sample)
        if ll is not None and ul is not None and not np.isnan(ll) and not np.isnan(ul):
            boot_lower_limits.append(ll)
            boot_upper_limits.append(ul)

    alpha_ci = 1.0 - ci_level
    ci_lower_percentile = (alpha_ci / 2.0) * 100
    ci_upper_percentile = (1.0 - alpha_ci / 2.0) * 100

    ll_ci = (
        np.percentile(boot_lower_limits, ci_lower_percentile) if boot_lower_limits else np.nan,
        np.percentile(boot_lower_limits, ci_upper_percentile) if boot_lower_limits else np.nan
    )
    ul_ci = (
        np.percentile(boot_upper_limits, ci_lower_percentile) if boot_upper_limits else np.nan,
        np.percentile(boot_upper_limits, ci_upper_percentile) if boot_upper_limits else np.nan
    )

    return {
        "lower_limit": lower_limit,
        "upper_limit": upper_limit,
        "lower_limit_ci": ll_ci,
        "upper_limit_ci": ul_ci,
        "method": "Robust"
    }


# ==========================================
# 3. Main Data Processing Pipeline
# ==========================================
def generate_reference_intervals(csv_path, parameters, output_excel, group_col='Group'):
    df = pd.read_csv(csv_path)
    groups = df[group_col].dropna().unique()
    results = []

    for param in parameters:
        row_data = {'Parameter': param}

        for group in groups:
            # Keep zeros this time to allow the one-sided logic to handle them
            group_data = df[df[group_col] == group][param].dropna().values
            n = len(group_data)

            ll, ul = np.nan, np.nan
            ll_ci_str, ul_ci_str = "N/A", "N/A"
            method_used = "Insufficient Data"

            if n > 0:
                try:
                    if n >= 120:
                        res = calculate_nonparametric_ri(group_data)
                    else:
                        res = calculate_robust_reference_interval(group_data)

                    ll = res['lower_limit']
                    ul = res['upper_limit']

                    # Format output cleanly
                    if res['method'] == "One-Sided Non-Parametric":
                        ll_ci_str = "N/A (0 bound)"
                    else:
                        ll_ci_str = f"{res['lower_limit_ci'][0]:.2f} - {res['lower_limit_ci'][1]:.2f}"

                    ul_ci_str = f"{res['upper_limit_ci'][0]:.2f} - {res['upper_limit_ci'][1]:.2f}"

                except Exception as e:
                    ll_ci_str, ul_ci_str = f"Error: {str(e)}", "Error"

            row_data[f"{group}_Lower_RI"] = ll
            row_data[f"{group}_Upper_RI"] = ul
            row_data[f"{group}_Lower_CI"] = ll_ci_str
            row_data[f"{group}_Upper_CI"] = ul_ci_str

        results.append(row_data)

    results_df = pd.DataFrame(results)
    results_df.set_index('Parameter', inplace=True)
    results_df.to_excel(output_excel)
    print(f"Reference intervals saved to {output_excel}")

if __name__ == "__main__":
    params = ['Segmented Neutrophil', 'Band Neutrophil', 'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Blast',
              'Monocyte', 'Lymphocyte', 'Large Granular Lymphocyte', 'Atypical Lymphocyte', 'Aberrant Lymphocyte',
              'Hairy Cell', 'Sezary Cell', 'Plasma Cell', 'Basophil', 'Eosinophil', 'Normoblast', 'Smudge Cell',
              'Platelet', 'Plt Clump', 'Giant Platelet', 'Large Platelet',
              'Hypochromatic', 'Macrocytes', 'Microcytes', 'Micro-organisms',
              'Poikilocytes', 'Polychromatic', 'Schistocytes+Helmet', 'Sickle', 'Spherocytes', 'Stomatocytes', 'Target', 'Tear Drop']
    generate_reference_intervals("raw/normals_CBM_raw.csv", params, "results/output.xlsx")
