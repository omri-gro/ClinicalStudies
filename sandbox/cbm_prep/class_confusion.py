import pandas as pd
import numpy as np
from itertools import combinations
from scipy.stats import pearsonr, t
import matplotlib.pyplot as plt
import seaborn as sns
from sandbox import MetadataBundle

if __name__ == "__main__":
    site = 'TASMC'
    triangle_matrix = False
    data_path = r'raw/side_analysis/5site_WBC.csv'
    group_tested = 'WBC diff'
    save_name = f'{group_tested.replace(" ", "_")}_{site}'

    ref_arm = 'OMR'
    test_arm = 'CBM'

    meta_path = r'config.yaml'
    metadata = MetadataBundle(meta_path)
    vars_to_test = set(metadata.variable_groups[group_tested])
    vars_dont_use = set(metadata.variable_groups['derived'])
    vars_to_test = vars_to_test - vars_dont_use - {"Poikilocytosis", "Poikilocytes", "Other WBC"}

    raw = pd.read_csv(data_path)

    if site != 'All':
        raw = raw.query("Site == @site")

    # find cell classes in combined dataframe
    ref_cols = [col for col in raw.columns if f'|{ref_arm}' in col]
    ref_classes = set([col.replace(f'|{ref_arm}', '') for col in ref_cols])
    test_cols = [col for col in raw.columns if f'|{test_arm}' in col]
    test_classes = set([col.replace(f'|{test_arm}', '') for col in test_cols])
    common_classes = test_classes.intersection(ref_classes)

    common_classes = sorted(list(vars_to_test.intersection(common_classes)))

    # calculate residuals for each combination
    residuals = pd.DataFrame(index=raw.index)
    for cls in common_classes:
        # treat non-numeric values like nan
        test_num = pd.to_numeric(raw[f'{cls}|{test_arm}'], errors='coerce')
        ref_num = pd.to_numeric(raw[f'{cls}|{ref_arm}'], errors='coerce')
        residuals[cls] = test_num - ref_num

    # calculate pairwise correlation of residuals
    cls_pairs = list(combinations(common_classes, 2))
    correlations = {}
    p_vals = {}
    mean_residual_matrix = pd.DataFrame(np.nan, index=common_classes, columns=common_classes)

    for cls1, cls2 in cls_pairs:
        valid_rows = residuals[[cls1, cls2]].dropna()
        if len(valid_rows) > 2:
            corr, p_val = pearsonr(valid_rows[cls1], valid_rows[cls2])
            correlations[(cls1, cls2)] = corr
            p_vals[(cls1, cls2)] = p_val
        else:
            correlations[(cls1, cls2)] = np.nan
            p_vals[(cls1, cls2)] = np.nan
        # Also calculate mean residual difference
        if len(valid_rows) > 0:
            mean_residual_diff = valid_rows[cls1].mean() - valid_rows[cls2].mean()
            mean_residual_matrix.loc[cls1, cls2] = mean_residual_diff
            mean_residual_matrix.loc[cls2, cls1] = -mean_residual_diff

    # organize as df
    corrs_df = pd.DataFrame.from_dict(correlations, orient='index', columns=['Correlations'])
    p_df = pd.DataFrame.from_dict(p_vals, orient='index', columns=['Correlations'])

    # calculate t-values according to degrees of freedom
    t_vals = {}
    for pair, corr in correlations.items():
        valid_rows = residuals[list(pair)].dropna()
        df = len(valid_rows) - 2
        if df > 0 and not np.isnan(corr):
            t_val = corr * np.sqrt(df / (1 - corr ** 2))
            t_vals[pair] = (t_val, df)
        else:
            t_vals[pair] = (np.nan, df)
    t_vals_df = pd.DataFrame.from_dict(t_vals, orient='index', columns=['T-Value', 'Degrees of Freedom'])

    # Visualize correlations as heatmap
    corr_mat = residuals.corr()  # calculate correlations again, this time as matrix
    corr_mat = corr_mat.reindex(index=common_classes, columns=common_classes)
    plt.figure(figsize=(10, 8))  # consider using matplotlib instead of seaborn
    if triangle_matrix:
        mask = np.triu(np.ones_like(corr_mat, dtype=bool))
        ax = sns.heatmap(corr_mat, annot=True, mask=mask, cmap='coolwarm', vmin=-1, vmax=1, fmt='.2f',
                         annot_kws={"fontsize": 6}, cbar_kws={"label": "Correlation"})
        mat_classes = corr_mat.index.tolist()
        ax.set_xticks(np.arange(len(mat_classes) - 1) + 0.5)
        ax.set_yticks(np.arange(len(mat_classes) - 1) + 1.5)
        ax.set_xticklabels(mat_classes[:-1], fontsize=7, rotation=50)  # Exclude last label
        ax.set_yticklabels(mat_classes[1:], fontsize=7)  # Exclude first label
    else:
        sns.heatmap(corr_mat, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt='.2f',
                    annot_kws={"fontsize": 6}, cbar_kws={"label": "Correlation"})
        plt.xticks(fontsize=7, rotation=50)
        plt.yticks(fontsize=7)
    plt.title(f"Correlation of Classes\n{save_name}")
    plt.savefig(f'{save_name}.jpg')
    plt.show()


