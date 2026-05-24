import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
import statsmodels.formula.api as smf
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

warnings.simplefilter('ignore', ConvergenceWarning)

# Pass/Fail limit dictionaries
# Format: { 'Official Name': (Switch_Threshold, SD_Limit, CV_Limit, Evaluation_Mode) }
AC_CONFIG = {
    # Tier 1: Ultra-Rare
    'Basophil': (2.0, 0.5, 30.0, 'HYBRID'),
    'Mast Cell': (2.0, 0.5, 30.0, 'HYBRID'),

    # Tier 2: Minor / Diagnostic
    'Blast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Promyelocyte': (5.0, 2.0, 30.0, 'HYBRID'),
    'Plasma Cell': (5.0, 2.0, 30.0, 'HYBRID'),
    'Erythroblast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Basophilic Normoblast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Monocyte': (5.0, 2.0, 30.0, 'HYBRID'),
    'Eosinophil': (5.0, 2.0, 30.0, 'HYBRID'),

    # Tier 3: Intermediate
    'Myelocyte': (10.0, 3.0, 25.0, 'HYBRID'),
    'Metamyelocyte': (10.0, 3.0, 25.0, 'HYBRID'),
    'Band Neutrophil': (10.0, 3.0, 25.0, 'HYBRID'),
    'Normoblast': (10.0, 3.0, 25.0, 'HYBRID'),
    'Polychromatophilic Normoblast': (10.0, 3.0, 25.0, 'HYBRID'),

    # Tier 4: Major Populations
    'Lymphocyte': (20.0, 5.0, 25.0, 'HYBRID'),
    'Segmented Neutrophil': (20.0, 5.0, 25.0, 'HYBRID'),

    # Safety fallback (If a name doesn't match perfectly)
    'default': (5.0, 2.0, 30.0, 'HYBRID')
}

# --- HELPER FUNCTIONS ---
def fit_mixed_model(data, param, top_group_col, nested_group_col):
    """
    Fits the nested ANOVA mixed model and extracts the 3 variance components.
    This safely handles the statsmodels implementation, singular matrices, and empty data.

    Returns:
        tuple: (top_level_variance, nested_level_variance, residual_variance)
    """
    try:
        # Use Q('Column Name') to escape spaces
        formula = f"Q('{param}') ~ 1"

        # The 're_formula="1"' ensures the top group intercept is maintained
        model = smf.mixedlm(
            formula,
            data,
            groups=data[top_group_col],
            re_formula="1",
            vc_formula={"Nested": f"0 + C({nested_group_col})"}
        )
        fit = model.fit(method='powell')

        v_top = max(0, fit.cov_re.iloc[0, 0])
        v_nested = max(0, fit.vcomp[0])
        v_residual = max(0, fit.scale)

    except (ValueError, np.linalg.LinAlgError, KeyError, IndexError):
        # Fallback: if the model is mathematically singular/empty, variation is pushed to residual
        v_top, v_nested, v_residual = 0, 0, data[param].var()

    return v_top, v_nested, v_residual


def calculate_ci(sd_val, mean_val, df, conf_level=0.95):
    """
    Calculates the 95% Confidence Interval (Lower, Upper) for SD and %CV
    using the Chi-Square distribution per CLSI EP05-A3 Eq. 11.
    """
    if df <= 0 or sd_val == 0:
        return "N/A", "N/A"

    alpha = 1 - conf_level

    # Chi-Square critical values
    chi2_lower = stats.chi2.ppf(1 - alpha / 2, df)
    chi2_upper = stats.chi2.ppf(alpha / 2, df)

    # Calculate SD Confidence Limits
    sd_lower = sd_val * np.sqrt(df / chi2_lower)
    sd_upper = sd_val * np.sqrt(df / chi2_upper)

    # Translate to %CV Confidence Limits
    cv_lower = (sd_lower / mean_val) * 100 if mean_val != 0 else 0
    cv_upper = (sd_upper / mean_val) * 100 if mean_val != 0 else 0

    # Format as string for the output table
    sd_ci_str = f"[{sd_lower:.4f}, {sd_upper:.4f}]"
    cv_ci_str = f"[{cv_lower:.2f}%, {cv_upper:.2f}%]"

    return sd_ci_str, cv_ci_str


def calculate_dfs(data, top_col, nested_col, v_top, v_nested, v_residual):
    """
    Calculates the Repeatability DF and approximates the Total (Within-Lab/Repro) DF
    using the Satterthwaite formula per CLSI EP05-A3.
    """
    N = len(data)
    n_top_groups = data[top_col].nunique()
    n_nested_groups = data[nested_col].nunique()

    # 1. Repeatability DF
    df_residual = N - n_nested_groups

    # 2. Total DF (Satterthwaite Approximation)
    # Reconstruct Mean Squares (MS) and their DF from the variance components
    df_top = n_top_groups - 1
    df_nested = n_nested_groups - n_top_groups

    # Prevent division by zero if there's no nested grouping
    if df_nested <= 0 or df_top <= 0 or df_residual <= 0:
        return max(1, df_residual), max(1, df_residual)

    n_rep = N / n_nested_groups
    n_nested_per_top = n_nested_groups / n_top_groups

    ms_residual = v_residual
    ms_nested = (v_nested * n_rep) + ms_residual
    ms_top = (v_top * n_rep * n_nested_per_top) + ms_nested

    # Coefficients for expressing Total Variance as a linear combination of MS
    a_top = 1 / (n_rep * n_nested_per_top)
    a_nested = (1 / n_rep) - a_top
    a_residual = 1 - (1 / n_rep)

    # Satterthwaite Formula Numerator and Denominator
    numerator = (a_top * ms_top + a_nested * ms_nested + a_residual * ms_residual) ** 2
    denominator = ((a_top * ms_top) ** 2 / df_top) + \
                  ((a_nested * ms_nested) ** 2 / df_nested) + \
                  ((a_residual * ms_residual) ** 2 / df_residual)

    df_total = numerator / denominator if denominator > 0 else df_residual

    return df_residual, df_total


def evaluate_status(mean_val, sd_val, cv_val, param_name):
    # Get limits or use default
    thresh, sd_lim, cv_lim, mode = AC_CONFIG.get(param_name, AC_CONFIG['default'])

    if mode == 'SD_ONLY':
        status = 'Pass' if sd_val <= sd_lim else 'Fail'
        metric = f"SD (Limit {sd_lim})"
    else:
        # Hybrid Logic
        if mean_val < thresh:
            status = 'Pass' if sd_val <= sd_lim else 'Fail'
            metric = f"SD (Limit {sd_lim})"
        else:
            status = 'Pass' if cv_val <= cv_lim else 'Fail'
            metric = f"CV% (Limit {cv_lim}%)"

    return status, metric

def export_results(results_list, filename):
    """
    Handles formatting and saving the output.
    Extracted so we can easily switch to Excel or SQL exports in the future.
    """
    if results_list:
        df_out = pd.DataFrame(results_list)
        df_out.to_csv(filename, index=False)
        print(f"Generated: {filename}")


# --- MAIN PROCESSING PIPELINES ---

def process_repeatability(df, parameters):
    """
    Processes Single-Site data (Expected: Sample, Day, Run, Scan).
    Matches CLSI EP05-A3 Section 3.8.1
    """
    print("Processing Repeatability (Single-Site) Study...")

    # Create string columns for accurate categorical modeling
    df['Day_str'] = df['Day'].astype(str)
    df['Day_Run_str'] = df['Day'].astype(str) + "_" + df['Run'].astype(str)

    results = []

    for param in parameters:
        for sample, sample_data in df.groupby('Sample'):
            if sample_data[param].mean() == 0 or len(sample_data) < 5:
                continue

            mean_val = sample_data[param].mean()

            # Use the shared statistical engine
            v_day, v_run, v_scan = fit_mixed_model(sample_data, param, "Day_str", "Day_Run_str")
            v_wl = v_scan + v_run + v_day  # Total Within-Lab

            # Standard Deviations
            sd_scan, sd_run, sd_day, sd_wl = np.sqrt([v_scan, v_run, v_day, v_wl])

            # Calculate Degrees of Freedom
            df_rep, df_wl = calculate_dfs(sample_data, "Day_str", "Day_Run_str", v_day, v_run, v_scan)

            # Calculate Confidence Intervals
            ci_sd_rep, ci_cv_rep = calculate_ci(sd_scan, mean_val, df_rep)
            ci_sd_wl, ci_cv_wl = calculate_ci(sd_wl, mean_val, df_wl)

            # %CVs
            cv_scan = (sd_scan / mean_val) * 100
            cv_run = (sd_run / mean_val) * 100
            cv_day = (sd_day / mean_val) * 100
            cv_wl = (sd_wl / mean_val) * 100

            # Evaluate pass/fail status
            status, metric_used = evaluate_status(mean_val, sd_wl, cv_wl, param)

            results.append({
                'Parameter': param,
                'Sample': sample,
                'Mean': round(mean_val, 2),
                'Repeatability SD': round(sd_scan, 2),
                'Repeatability %CV': round(cv_scan, 2),
                'Between-Run SD': round(sd_run, 2),
                'Between-Run %CV': round(cv_run, 2),
                'Between-Day SD': round(sd_day, 2),
                'Between-Day %CV': round(cv_day, 2),
                'Within-Laboratory SD': round(sd_wl, 2),
                'Within-Laboratory %CV': round(cv_wl, 2),
                'Status': status,
                'Metric Evaluated': metric_used
            })

    export_results(results, f"results/Repeatability.csv")


def process_reproducibility(df, parameters):
    """
    Processes Multisite data (Expected: Sample, Machine, Day, Scan).
    Matches CLSI EP05-A3 Section 4.7
    """
    print("Processing Reproducibility (Multisite) Study...")

    # Create string columns
    df['Machine_str'] = df['Machine'].astype(str)
    df['Machine_Day_str'] = df['Machine'].astype(str) + "_" + df['Day'].astype(str)

    results = []

    for param in parameters:
        for sample, sample_data in df.groupby('Sample'):
            if sample_data[param].mean() == 0 or len(sample_data) < 5:
                continue

            mean_val = sample_data[param].mean()

            # Use the shared statistical engine
            v_machine, v_day, v_scan = fit_mixed_model(sample_data, param, "Machine_str", "Machine_Day_str")
            v_repro = v_scan + v_day + v_machine  # Total Reproducibility

            # Standard Deviations
            sd_scan, sd_day, sd_machine, sd_repro = np.sqrt([v_scan, v_day, v_machine, v_repro])

            # Calculate Degrees of Freedom
            df_rep, df_repro = calculate_dfs(sample_data, "Machine_str", "Machine_Day_str", v_machine, v_day, v_scan)

            # Calculate Confidence Intervals
            ci_sd_rep, ci_cv_rep = calculate_ci(sd_scan, mean_val, df_rep)
            ci_sd_repro, ci_cv_repro = calculate_ci(sd_repro, mean_val, df_repro)

            # %CVs
            cv_scan = (sd_scan / mean_val) * 100
            cv_day = (sd_day / mean_val) * 100
            cv_machine = (sd_machine / mean_val) * 100
            cv_repro = (sd_repro / mean_val) * 100

            # Evaluate pass/fail status
            status, metric_used = evaluate_status(mean_val, sd_repro, cv_repro, param)

            results.append({
                'Parameter': param,
                'Sample': sample,
                'Mean': round(mean_val, 4),
                'Repeatability SD': round(sd_scan, 4),
                'Repeatability %CV': round(cv_scan, 2),
                'Between-Day SD': round(sd_day, 4),
                'Between-Day %CV': round(cv_day, 2),
                'Between-Site SD': round(sd_machine, 4),
                'Between-Site %CV': round(cv_machine, 2),
                'Reproducibility SD': round(sd_repro, 4),
                'Reproducibility %CV': round(cv_repro, 2),
                'Status': status,
                'Metric Evaluated': metric_used
            })

        # export_results(results, f"results/Table_4.7_Reproducibility_{param}.csv")
    export_results(results, f"results/Reproducibility.csv")
