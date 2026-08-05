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
STUDY_MODE = "Spherocytes"  # RBC_Shape, Spherocytes, RBC_SIZE or PLT_SIZE
FILE_PATH = rf"results/{STUDY_MODE}_results.csv"
REF_COL = "percent_in_roi_rev"
TEST_COL = "percent_test_arm_ai"  # Options: 'percent_test_arm_ai', 'percent_in_scan_ai' or 'percent_in_roi_ai'

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
    'RBC macrocyte': 5.0,
    'RBC microcyte': 5.0,
    # Add other morphologies and their thresholds here. Default will be > 0 if not listed.
}
DEFAULT_THRESHOLD = 1.0

# Plotting Options
GENERATE_PLOTS = True  # Toggle scatter plot generation
GENERATE_ERROR_BAR_PLOTS = True  # Toggle scatter plots with inter-reviewer disagreement error bars
EQUAL_AXES = True  # If True, x and y axes use the same scale and an x=y identity line is drawn
PLOT_DIR = "plots"  # Folder to save the generated plots

# Output Options
SAVE_CSV = True  # Toggle saving the summary table to a CSV
CSV_FILENAME = rf"results/{STUDY_MODE}_regression_results.csv"


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


def format_plot_axes(x, y, slope, intercept, title, x_col, y_col):
    """Helper function to format axes, limits, identity line, and Deming line."""
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

    plt.title(title)
    plt.xlabel(f'Reference ({x_col})')
    plt.ylabel(f'AI Test ({y_col})')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()


def create_plot(df, morph, x_col, y_col, group_col, slope, intercept, title_suffix=""):
    """Generates and saves the standard method comparison scatter plot."""
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

    plt.figure(figsize=(8, 8) if EQUAL_AXES else (8, 6))
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=group_col, palette='Set1', s=60, zorder=3)

    title = f'Method Comparison for {morph}'
    if EQUAL_AXES: title += '\n(Equal Axes)'
    if title_suffix: title += f'\n{title_suffix}'

    format_plot_axes(df[x_col].values, df[y_col].values, slope, intercept, title, x_col, y_col)

    safe_morph = morph.replace("&", "_and_").replace(" ", "_")
    plt.savefig(f"{PLOT_DIR}/{safe_morph}_comparison.png")
    plt.close()


def create_error_bar_plot(df, morph, x_col, y_col, group_col, slope, intercept, title_suffix=""):
    """Generates and saves a scatter plot including + shaped error bars indicating individual reviewer spread."""
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

    plt.figure(figsize=(8, 8) if EQUAL_AXES else (8, 6))

    # Define color palette matching Seaborn's Set1
    groups = df[group_col].unique()
    palette = sns.color_palette('Set1', n_colors=len(groups))
    color_map = dict(zip(groups, palette))

    # Plot error bars behind the scatter points
    for group_name in groups:
        group_data = df[df[group_col] == group_name]

        plt.errorbar(
            x=group_data[x_col],
            y=group_data[y_col],
            xerr=group_data['ref_err'],
            yerr=group_data['test_err'],
            fmt='none',
            ecolor=color_map[group_name],
            elinewidth=1.5,
            alpha=0.4,  # Keep bars semi-transparent so they don't clutter the view
            zorder=1
        )

    # Overlay the scatter points
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=group_col, palette='Set1', s=60, zorder=3)

    title = f'Method Comparison (with inter-reviewer variance) for {morph}'
    if EQUAL_AXES: title += '\n(Equal Axes)'
    if title_suffix: title += f'\n{title_suffix}'

    format_plot_axes(df[x_col].values, df[y_col].values, slope, intercept, title, x_col, y_col)

    safe_morph = morph.replace("&", "_and_").replace(" ", "_")
    plt.savefig(f"{PLOT_DIR}/{safe_morph}_comparison_errorbars.png")
    plt.close()


def main():
    print(f"Loading data from {FILE_PATH}...")
    try:
        df = pd.read_csv(FILE_PATH)
    except FileNotFoundError:
        print(f"Error: Could not find the file '{FILE_PATH}'.")
        return

    if TEST_COL == "percent_test_arm_ai":
        if 'test_arm' not in df.columns:
            print("Error: The column 'test_arm' was not found in the CSV. Please add it to use 'percent_test_arm_ai'.")
            return

        conditions = [
            df['test_arm'].astype(str).str.strip().str.lower() == 'scan',
            df['test_arm'].astype(str).str.strip().str.lower() == 'roi'
        ]
        choices = [
            df['percent_in_scan_ai'],
            df['percent_in_roi_ai']
        ]

        df[TEST_COL] = np.select(conditions, choices, default=np.nan)

    # Identify reviewer column
    rev_col = 'Reviewer Name' if 'Reviewer Name' in df.columns else 'Reviewer'

    # PREPROCESSING
    if AGGREGATE_BY_SAMPLE:  # might need to remove 'Site' from list of columns when RBC_shape
        agg_df = df.groupby(['Internal Count', 'Morphology']).agg(
            ref_mean=(REF_COL, 'mean'),
            ref_min=(REF_COL, 'min'),
            ref_max=(REF_COL, 'max'),
            test_mean=(TEST_COL, 'mean'),
            test_min=(TEST_COL, 'min'),
            test_max=(TEST_COL, 'max'),
            Grouping=(rev_col, get_pair_name),
            Count=(REF_COL, 'count')
        ).reset_index()

        # Calculate standard symmetric errors for plotting (max - mean is identical to mean - min when n=2)
        agg_df['ref_err'] = agg_df['ref_max'] - agg_df['ref_mean']
        agg_df['test_err'] = agg_df['test_max'] - agg_df['test_mean']

        # Rename back to maintain compatibility with existing analysis blocks
        agg_df.rename(columns={'ref_mean': REF_COL, 'test_mean': TEST_COL}, inplace=True)

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

        combined_slope = np.nan
        combined_intercept = np.nan

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

            if group_name == 'All Combined':
                combined_slope = slope
                combined_intercept = intercept

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
            create_plot(morph_data, morph, REF_COL, TEST_COL, 'Grouping', combined_slope, combined_intercept,
                        title_suffix=f"Test Arm: {TEST_COL}")

        if GENERATE_ERROR_BAR_PLOTS and AGGREGATE_BY_SAMPLE:
            create_error_bar_plot(morph_data, morph, REF_COL, TEST_COL, 'Grouping', combined_slope, combined_intercept,
                                  title_suffix=f"Test Arm: {TEST_COL}")

    # OUTPUT
    res_df = pd.DataFrame(results)

    print("\n--- Final Results ---\n")
    print(res_df.to_string(index=False))

    if SAVE_CSV:
        res_df.to_csv(CSV_FILENAME, index=False)
        print(f"\nResults successfully saved to '{CSV_FILENAME}'")

    if GENERATE_PLOTS or GENERATE_ERROR_BAR_PLOTS:
        print(f"Plots successfully saved to the '{PLOT_DIR}/' directory.")


if __name__ == "__main__":
    main()
