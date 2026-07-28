import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.odr import Model, Data, ODR
import os
import warnings

# Suppress runtime warnings from zero-variance bootstrap samples
warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURATION OPTIONS (Boolean toggles and parameters)
# ==============================================================================

# Data Parameters
FILE_PATH = r"results/RBC_Shape_results.csv"
REF_COL = "percent_in_roi_rev"
TEST_COL = "percent_in_scan_ai"  # Options: 'percent_in_scan_ai' or 'percent_in_roi_ai'

# Analysis Modifiers
AGGREGATE_BY_SAMPLE = True  # If True, calculates per-sample mean. If False, analyzes per-reviewer rows.
REQUIRE_TWO_REVIEWERS = True  # If True (and AGGREGATE is True), only keeps samples evaluated by exactly 2 reviewers.
GROUP_BY_PAIR_OR_REVIEWER = True  # If True, breaks down results by Reviewer Pair (if aggregated) or individual Reviewer (if not).

# Statistical Toggles
CALC_BOOTSTRAP_CI = True  # Toggle Bootstrap Confidence Intervals (takes a bit longer to run)
N_BOOTSTRAP = 2000  # Number of resampling iterations for CI calculation

# Positivity Thresholds (Morphology Name : Threshold Value)
# Defines what counts as a "Positive" for the "Number of Positives" metric based on the reference column.
THRESHOLDS = {
    'Spherocyte': 5.0,
    'Helmet&Schisto': 0.5,
    # Add other morphologies and their thresholds here. Default will be > 0 if not listed.
}
DEFAULT_THRESHOLD = 1.0

# Plotting Options
GENERATE_PLOTS = True  # Toggle scatter plot generation
EQUAL_AXES = True  # If True, x and y axes use the same scale and an x=y identity line is drawn
PLOT_DIR = "plots"  # Folder to save the generated plots

# Output Options
SAVE_CSV = True  # Toggle saving the summary table to a CSV
CSV_FILENAME = r"results/method_comparison_results.csv"


# ==============================================================================

def get_pair_name(reviewers):
    """Sorts and formats reviewer names into a standard pair string."""
    rev_list = sorted(list(set(reviewers.dropna())))
    return " & ".join(rev_list)


def linear_func(p, x):
    """Linear model for ODR (Deming). p[0] is slope, p[1] is intercept."""
    return p[0] * x + p[1]


def run_deming(x, y):
    """Runs standard Deming Regression."""
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan

    linear_model = Model(linear_func)
    odr_data = Data(x, y)
    odr_obj = ODR(odr_data, linear_model, beta0=[1.0, 0.0])
    try:
        res = odr_obj.run()
        return res.beta[0], res.beta[1]  # Slope, Intercept
    except:
        return np.nan, np.nan


def bootstrap_deming_ci(x, y, n_boot=2000, seed=42):
    """Calculates 95% CIs for slope and intercept using bootstrap resampling."""
    np.random.seed(seed)
    n = len(x)
    slopes, intercepts = [], []
    linear_model = Model(linear_func)

    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        x_b, y_b = x[idx], y[idx]

        if np.std(x_b) == 0 or np.std(y_b) == 0:
            continue

        odr_data = Data(x_b, y_b)
        odr_obj = ODR(odr_data, linear_model, beta0=[1.0, 0.0])
        try:
            res = odr_obj.run()
            slopes.append(res.beta[0])
            intercepts.append(res.beta[1])
        except:
            continue

    if len(slopes) < 10:  # Not enough successful bootstraps
        return np.nan, np.nan, np.nan, np.nan

    slope_ci_lo, slope_ci_hi = np.percentile(slopes, 2.5), np.percentile(slopes, 97.5)
    int_ci_lo, int_ci_hi = np.percentile(intercepts, 2.5), np.percentile(intercepts, 97.5)
    return slope_ci_lo, slope_ci_hi, int_ci_lo, int_ci_hi


def create_plot(df, morph, x_col, y_col, group_col, slope, intercept, title_suffix=""):
    """Generates and saves the method comparison scatter plot."""
    if not os.path.exists(PLOT_DIR):
        os.makedirs(PLOT_DIR)

    plt.figure(figsize=(8, 8) if EQUAL_AXES else (8, 6))

    x = df[x_col].values
    y = df[y_col].values

    sns.scatterplot(data=df, x=x_col, y=y_col, hue=group_col, palette='Set1', s=60)

    axis_min, axis_max = min(x), max(x)  # Defaults

    if EQUAL_AXES:
        min_val = min(np.min(x), np.min(y))
        max_val = max(np.max(x), np.max(y))
        padding = (max_val - min_val) * 0.05 if max_val != min_val else 0.1
        axis_min, axis_max = min(0, min_val - padding), max_val + padding

        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.plot([axis_min, axis_max], [axis_min, axis_max], color='gray', linestyle=':', linewidth=2,
                 label='x = y (Identity)')

    if not np.isnan(slope):
        x_vals = np.linspace(axis_min, axis_max, 100)
        y_vals = slope * x_vals + intercept
        plt.plot(x_vals, y_vals, color='black', linestyle='--', linewidth=2, label=f'Deming (slope={slope:.2f})')

    title = f'Method Comparison for {morph}'
    if EQUAL_AXES: title += '\n(Equal Axes)'
    if title_suffix: title += f'\n{title_suffix}'

    plt.title(title)
    plt.xlabel(f'Reference ({x_col})')
    plt.ylabel(f'AI Test ({y_col})')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    safe_morph = morph.replace("&", "_and_").replace(" ", "_")
    filename = f"{PLOT_DIR}/{safe_morph}_comparison.png"
    plt.savefig(filename)
    plt.close()


def main():
    print(f"Loading data from {FILE_PATH}...")
    df = pd.read_csv(FILE_PATH)

    # Identify reviewer column
    rev_col = 'Reviewer Name' if 'Reviewer Name' in df.columns else 'Reviewer'

    # PREPROCESSING
    if AGGREGATE_BY_SAMPLE:
        agg_df = df.groupby(['Internal Count', 'Morphology']).agg(
            **{REF_COL: (REF_COL, 'mean')},
            **{TEST_COL: (TEST_COL, 'mean')},
            Grouping=(rev_col, get_pair_name),
            Count=(REF_COL, 'count')
        ).reset_index()

        if REQUIRE_TWO_REVIEWERS:
            agg_df = agg_df[agg_df['Count'] == 2]
    else:
        agg_df = df.copy()
        agg_df['Grouping'] = agg_df[rev_col]
        agg_df = agg_df.dropna(subset=[REF_COL, TEST_COL])

    results = []
    morphologies = agg_df['Morphology'].unique()

    # ANALYSIS LOOP
    print("Starting analysis...")
    for morph in morphologies:
        morph_data = agg_df[agg_df['Morphology'] == morph]
        threshold = THRESHOLDS.get(morph, DEFAULT_THRESHOLD)

        # Determine subsets to iterate over (All Combined + Individual Groups if toggled)
        subsets = [('All Combined', morph_data)]
        if GROUP_BY_PAIR_OR_REVIEWER:
            for group_name in morph_data['Grouping'].unique():
                subsets.append((group_name, morph_data[morph_data['Grouping'] == group_name]))

        for group_name, subset_data in subsets:
            x = subset_data[REF_COL].values
            y = subset_data[TEST_COL].values
            n = len(x)

            num_pos = np.sum(x >= threshold)

            # Pearson r
            r = np.nan
            if n > 1 and np.std(x) > 0 and np.std(y) > 0:
                r, _ = stats.pearsonr(x, y)

            # Deming
            slope, intercept = run_deming(x, y)

            # Bootstrap CI
            slope_ci_lo, slope_ci_hi, int_ci_lo, int_ci_hi = np.nan, np.nan, np.nan, np.nan
            if CALC_BOOTSTRAP_CI and n > 2 and not np.isnan(slope):
                slope_ci_lo, slope_ci_hi, int_ci_lo, int_ci_hi = bootstrap_deming_ci(x, y, N_BOOTSTRAP)

            # Format results
            results.append({
                'Morphology': morph,
                'Group': group_name,
                'N': n,
                f'Positives (>= {threshold})': num_pos,
                'Pearson r': f"{r:.4f}" if not np.isnan(r) else "NaN",
                'Slope': f"{slope:.4f}" if not np.isnan(slope) else "NaN",
                'Slope 95% CI': f"[{slope_ci_lo:.4f}, {slope_ci_hi:.4f}]" if not np.isnan(slope_ci_lo) else "NaN",
                'Intercept': f"{intercept:.4f}" if not np.isnan(intercept) else "NaN",
                'Intercept 95% CI': f"[{int_ci_lo:.4f}, {int_ci_hi:.4f}]" if not np.isnan(int_ci_lo) else "NaN"
            })

        # Plotting (Done per morphology on the whole set)
        if GENERATE_PLOTS:
            create_plot(morph_data, morph, REF_COL, TEST_COL, 'Grouping', slope, intercept,
                        title_suffix=f"Test Arm: {TEST_COL}")

    # OUTPUT
    res_df = pd.DataFrame(results)

    print("\n--- Final Results ---\n")
    print(res_df.to_string(index=False))

    if SAVE_CSV:
        res_df.to_csv(CSV_FILENAME, index=False)
        print(f"\nResults successfully saved to '{CSV_FILENAME}'")

    if GENERATE_PLOTS:
        print(f"Plots successfully saved to the '{PLOT_DIR}/' directory.")


if __name__ == "__main__":
    main()
