import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.odr import Model, Data, ODR
import matplotlib.pyplot as plt


def get_deming(x_val, y_val):
    linear_model = Model(lambda p, x: p[0] * x + p[1])
    odr_res = ODR(Data(x_val, y_val), linear_model, beta0=[1., 0.]).run()
    return odr_res.beta[0], odr_res.beta[1]


def bootstrap_deming(x_val, y_val, n_bootstraps=1000):
    slopes = []
    intercepts = []
    n = len(x_val)
    x_val = np.array(x_val)
    y_val = np.array(y_val)

    np.random.seed(42)
    for _ in range(n_bootstraps):
        indices = np.random.choice(n, size=n, replace=True)
        x_boot = x_val[indices]
        y_boot = y_val[indices]
        if len(np.unique(x_boot)) > 1:
            m, b = get_deming(x_boot, y_boot)
            slopes.append(m)
            intercepts.append(b)

    slope_ci = np.percentile(slopes, [2.5, 97.5])
    intercept_ci = np.percentile(intercepts, [2.5, 97.5])
    return slope_ci, intercept_ci


def analyze_and_plot(df, morphology, x, y, scatter_filename, csv_filename):
    col_cbm = f'{morphology}|CBM'
    col_clv = f'{morphology}|ClV'
    poikilo_col = 'Poikilocytes|CBM'

    # --- 1. Calculate the conversion logic ---
    # Determine the amount of poikilocytes to convert for each row
    condition = df[col_cbm] > x
    converted_amount = np.where(condition, y * df[poikilo_col], 0)

    # Apply logic to create new values
    modified_cbm = df[col_cbm] + converted_amount
    modified_poikilo = df[poikilo_col] - converted_amount

    # --- 2. Save modified raw data to CSV ---
    out_df = pd.DataFrame({
        'Site': df['Site'],
        'SampleID': df['SampleID'],
        f'Original {morphology}|CBM': df[col_cbm],
        f'New {morphology}|CBM': modified_cbm,
        f'{morphology}|ClV': df[col_clv],
        f'Original {poikilo_col}': df[poikilo_col],
        f'New {poikilo_col}': modified_poikilo
    })
    out_df.to_csv(csv_filename, index=False)
    print(f"Saved modified data to {csv_filename}")

    # --- 3. Calculate Stats ---
    # Filter NaNs for plotting and math
    mask = ~np.isnan(df[col_clv]) & ~np.isnan(modified_cbm)
    x_data = df[col_clv][mask].values
    y_data = modified_cbm[mask]

    r, _ = pearsonr(x_data, y_data)
    slope, intercept = get_deming(x_data, y_data)
    slope_ci, intercept_ci = bootstrap_deming(x_data, y_data, n_bootstraps=500)

    # --- 4. Setup Plot ---
    plt.figure(figsize=(8, 8))
    plt.scatter(x_data, y_data, alpha=0.5, label='Data points', color='blue', edgecolor='k')

    max_val = max(np.max(x_data), np.max(y_data))
    min_val = min(np.min(x_data), np.min(y_data))
    padding = (max_val - min_val) * 0.05
    axis_max = max_val + padding
    axis_min = max_val * -0.05 if min_val == 0 else min_val - padding

    plt.xlim(axis_min, axis_max)
    plt.ylim(axis_min, axis_max)
    plt.gca().set_aspect('equal', adjustable='box')

    line_x = np.array([axis_min, axis_max])
    plt.plot(line_x, slope * line_x + intercept, color='red', label=f'Deming Regression (m={slope:.2f})')
    plt.plot(line_x, line_x, color='black', linestyle='--', label='y = x')

    plt.title(f'Method Comparison: {morphology}\n(x={x}, y={y})')
    plt.xlabel('Reference (ClV)')
    plt.ylabel('Test (Modified CBM)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(scatter_filename)
    plt.close()


# --- Execution ---
df = pd.read_csv(r"../../comp_tables/clv_cbm_all-ssn_mininv-2_no_scrtch-False_brdrmv-False_v325_NoTechFlgs_mean_inv.csv")

# Generate outputs for Sickle Cells
analyze_and_plot(df, "Sickle cells", 0, 0, r'results/sickle_scatter.png', r'results/sickle_modified_data.csv')

# Generate outputs for Tear Drop Cells
analyze_and_plot(df, "Tear drop cells", 0, 0, r'results/teardrop_scatter.png', r'results/teardrop_modified_data.csv')

