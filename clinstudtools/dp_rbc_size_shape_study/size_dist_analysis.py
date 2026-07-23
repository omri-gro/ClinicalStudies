import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. Configuration
# ==========================================
STUDY_MODE = "PLT_SIZE"  # "PLT_SIZE" or "RBC_SIZE"
DATA_FILE = rf"./side_results/behavioral_cell_data_{STUDY_MODE}.csv"
OUTPUT_DIR = "./visualizations"

# Set visualization style
sns.set_theme(style="whitegrid", palette="muted")


# ==========================================
# 2. Plotting Functions
# ==========================================
def plot_aggregate_distribution(df, size_col, title, output_filename, valid_actions):
    """Plots the overall KDE for all reviewers combined."""
    plt.figure(figsize=(10, 6))

    sns.kdeplot(
        data=df,
        x=size_col,
        hue="Reviewer_Action",
        hue_order=valid_actions,  # <-- LOCKS THE COLOR ORDER
        fill=True,
        common_norm=False,
        alpha=0.4,
        linewidth=2
    )

    plt.title(title, fontsize=14, pad=15)
    plt.xlabel("AI Calculated Size (μm)", fontsize=12)
    plt.ylabel("Density", fontsize=12)

    plt.xlim(df[size_col].quantile(0.001), df[size_col].quantile(0.999))
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()


def plot_per_reviewer_distribution(df, size_col, title, output_filename, valid_actions):
    """Plots a FacetGrid showing individual KDEs for each reviewer."""

    # Move hue and hue_order to the FacetGrid level so it enforces it globally
    g = sns.FacetGrid(
        df,
        col="Reviewer_Name",
        col_wrap=3,
        height=4,
        aspect=1.2,
        sharey=False,
        hue="Reviewer_Action",
        hue_order=valid_actions,  # <-- LOCKS THE COLOR ORDER
        palette="muted"  # <-- ENFORCES THE SAME PALETTE
    )

    g.map_dataframe(
        sns.kdeplot,
        x=size_col,
        fill=True,
        common_norm=False,
        alpha=0.4,
        linewidth=1.5
    )

    g.add_legend(title="Action")
    g.set_axis_labels("AI Calculated Size (μm)", "Density")
    g.set_titles(col_template="{col_name}")
    g.fig.suptitle(title, fontsize=16, y=1.05)

    g.set(xlim=(df[size_col].quantile(0.001), df[size_col].quantile(0.999)))

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================
# 3. Execution Wrapper
# ==========================================
if __name__ == "__main__":
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Loading data from {DATA_FILE}...")
    try:
        df = pd.read_csv(DATA_FILE)
    except FileNotFoundError:
        print(f"Error: Could not find {DATA_FILE}. Run behavioral_extraction.py first.")
        exit()

    # Drop rows where the AI bounding box failed to calculate a size
    df = df.dropna(subset=['Cell_Size_Max_um', 'Cell_Size_Mean_um'])

    # Determine which actions to plot based on study mode to keep legends clean
    if STUDY_MODE == "PLT_SIZE":
        valid_actions = ["Ignored", "Marked as PLT large", "Marked as PLT giant"]
    elif STUDY_MODE == "RBC_SIZE":
        valid_actions = ["Ignored", "Marked as RBC macrocyte", "Marked as RBC microcyte"]
    else:
        valid_actions = df['Reviewer_Action'].unique().tolist()

    df_filtered = df[df['Reviewer_Action'].isin(valid_actions)]

    print("Generating aggregate plots...")
    # 1. Aggregate Plots
    plot_aggregate_distribution(
        df_filtered,
        size_col="Cell_Size_Max_um",
        title=f"Aggregate Cell Size Distribution (Max Width/Height) - {STUDY_MODE}",
        output_filename=os.path.join(OUTPUT_DIR, f"{STUDY_MODE}_aggregate_max.png"),
        valid_actions=valid_actions
    )

    plot_aggregate_distribution(
        df_filtered,
        size_col="Cell_Size_Mean_um",
        title=f"Aggregate Cell Size Distribution (Mean Width/Height) - {STUDY_MODE}",
        output_filename=os.path.join(OUTPUT_DIR, f"{STUDY_MODE}_aggregate_mean.png"),
        valid_actions=valid_actions
    )

    print("Generating per-reviewer plots...")
    # 2. Per-Reviewer Plots
    plot_per_reviewer_distribution(
        df_filtered,
        size_col="Cell_Size_Max_um",
        title=f"Reviewer Bias: Cell Size Distribution (Max Width/Height) - {STUDY_MODE}",
        output_filename=os.path.join(OUTPUT_DIR, f"{STUDY_MODE}_per_reviewer_max.png"),
        valid_actions=valid_actions
    )

    plot_per_reviewer_distribution(
        df_filtered,
        size_col="Cell_Size_Mean_um",
        title=f"Reviewer Bias: Cell Size Distribution (Mean Width/Height) - {STUDY_MODE}",
        output_filename=os.path.join(OUTPUT_DIR, f"{STUDY_MODE}_per_reviewer_mean.png"),
        valid_actions=valid_actions
    )

    print(f"Complete! Visualizations saved to the '{OUTPUT_DIR}' directory.")