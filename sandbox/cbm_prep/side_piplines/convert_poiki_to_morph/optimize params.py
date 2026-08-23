import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.odr import Model, Data, ODR
import matplotlib.pyplot as plt
import seaborn as sns


def get_deming(x_val, y_val):
    """Calculates Deming Regression slope and intercept."""
    linear_model = Model(lambda p, x: p[0] * x + p[1])
    odr_res = ODR(Data(x_val, y_val), linear_model, beta0=[1., 0.]).run()
    return odr_res.beta[0], odr_res.beta[1]


# --- 1. Load Data ---
df = pd.read_csv(r"../../comp_tables/clv_cbm_all-ssn_mininv-2_no_scrtch-False_brdrmv-False_v325_NoTechFlgs_mean_inv.csv")

# --- 2. Define Parameters ---
morphologies = ["Sickle cells", "Tear drop cells"]
poikilo_col = 'Poikilocytes|CBM'

x_range = np.arange(0.5, 5.5, 0.5)
y_range = np.arange(0.0, 0.55, 0.05)

all_results = []

# --- 3. Iterate & Calculate ---
for morph in morphologies:
    col_cbm = f'{morph}|CBM'
    col_clv = f'{morph}|ClV'

    # Check column existence
    if col_cbm not in df.columns or col_clv not in df.columns:
        print(f"Skipping {morph}: Columns not found in dataset.")
        continue

    for x in x_range:
        for y in y_range:
            modified_cbm = np.where(df[col_cbm] > x,
                                    df[col_cbm] + (y * df[poikilo_col]),
                                    df[col_cbm])

            mask = ~np.isnan(df[col_clv]) & ~np.isnan(modified_cbm)
            if mask.sum() > 2:
                r, _ = pearsonr(df[col_clv][mask], modified_cbm[mask])
                slope, intercept = get_deming(df[col_clv][mask], modified_cbm[mask])

                # Check criteria & score distance from ideal (r=1, slope=1)
                meets_criteria = (r > 0.7) and (0.7 <= slope <= 1.3)
                score = (1 - r) ** 2 + (1 - slope) ** 2

                all_results.append({
                    'Morphology': morph,
                    'x': x,
                    'y': round(y, 2),
                    'Pearson_r': r,
                    'Slope': slope,
                    'Intercept': intercept,
                    'Score': score,
                    'Meets_Criteria': meets_criteria
                })

res_df = pd.DataFrame(all_results)

# Print Top Results
success_df = res_df[res_df['Meets_Criteria']]
for morph in morphologies:
    print(f"\n--- Top 3 Parameters for {morph} ---")
    top_morph = success_df[success_df['Morphology'] == morph].sort_values('Score').head(3)
    if top_morph.empty:
        print("No parameters met criteria.")
    else:
        print(top_morph[['x', 'y', 'Pearson_r', 'Slope', 'Intercept']].to_string(index=False))

# --- 4. Generate Heatmaps ---
# Dynamically size figure based on number of morphologies analyzed
fig, axes = plt.subplots(len(morphologies), 2, figsize=(14, 5 * len(morphologies)))
fig.suptitle('Parameter Optimization Heatmaps', fontsize=16)

# Handle edge case where there is only one morphology (axes is 1D)
if len(morphologies) == 1:
    axes = np.array([axes])

for i, morph in enumerate(morphologies):
    sub_df = res_df[res_df['Morphology'] == morph]
    if sub_df.empty:
        continue

    # Pivot for Pearson r
    pivot_r = sub_df.pivot(index='x', columns='y', values='Pearson_r')
    sns.heatmap(pivot_r, ax=axes[i, 0], cmap='viridis', annot=False, cbar=True, vmin=0.6, vmax=1.0)
    axes[i, 0].set_title(f'{morph} - Pearson r')
    axes[i, 0].set_ylabel('x (Threshold)')
    axes[i, 0].set_xlabel('y (Fraction of Poikilocytes)')
    axes[i, 0].invert_yaxis()  # Put lowest x value at the bottom

    # Pivot for Deming Slope
    pivot_slope = sub_df.pivot(index='x', columns='y', values='Slope')
    sns.heatmap(pivot_slope, ax=axes[i, 1], cmap='coolwarm', center=1.0, annot=False, cbar=True, vmin=0.5, vmax=1.5)
    axes[i, 1].set_title(f'{morph} - Deming Slope')
    axes[i, 1].set_ylabel('x (Threshold)')
    axes[i, 1].set_xlabel('y (Fraction of Poikilocytes)')
    axes[i, 1].invert_yaxis()

plt.tight_layout()
plt.subplots_adjust(top=0.92)  # Ensure title isn't clipped
plt.savefig(r'results/optimization_heatmaps.png')
plt.show()
