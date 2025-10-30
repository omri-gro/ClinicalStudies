import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from itertools import product
from scipy.stats import spearmanr, kendalltau
from statsmodels.miscmodels.ordinal_model import OrderedModel


if __name__ == "__main__":
    meta_path = r'config.yaml'
    save_name = 'all_sites_combined'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    # build list of files to read from
    sites = ['BWH', 'CPG', 'LMU', 'SYN']
    mthds = ['OMR', 'CBM']
    srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')

    graded_params = metadata.variable_groups['grade']
    is_in_list = methd_comp.df['Variable'].isin(graded_params)
    grades_df = methd_comp.df[is_in_list]

    # get only samples that have both OMR and CBM results
    df = grades_df.groupby(['SampleID', "Site", "Variable"]).filter(lambda g: g["Method"].nunique() == 2)

    # drop rows without numeric grade
    df["Grade"] = pd.to_numeric(df["Grade"], errors="coerce")
    df = df.dropna(subset=["Grade"])

    # for values that were already supplied as grades (so all RBC values from OMR), Grade is equal is value
    df = df.reset_index(drop=True)
    df.loc[df["Method"] == "OMR", "Grade"] = df["Value"]

    # convert BWH CBM values from a 0-1 to a 0-100 scale
    df.loc[(df["Site"] == "BWH") & (df["Method"] == "CBM"), "Value"] = 100 * df["Value"]

    # pivot both Value and Grade by Method
    wide = (
        df.pivot_table(
            index=['SampleID', "Site", "Variable"],
            columns='Method',
            values=["Value", "Grade"],
            aggfunc="first"
        )
        .reset_index()
    )

    # flatten MultiIndex columns
    wide.columns = [
        "Sample", "Site", "Variable",
        *[f"{v}_{m}" for v, m in wide.columns.tolist()[3:]]
    ]


    # get info on Value_CBM for each Site-Variable-Grade_OMR combination
    wide["Value_CBM"] = pd.to_numeric(wide["Value_CBM"], errors="coerce")
    summary = (
        wide.groupby(["Site", "Variable", "Grade_OMR"])["Value_CBM"]
        .agg(["size", "mean", "median", "min", "max"])
        .reset_index()
        .sort_values(["Variable", "Site", "Grade_OMR"])
        .round({col: 2 for col in ["mean", "median", "min", "max"]})
    )

    pd.DataFrame(summary).to_csv('omr_grades_summary.csv', index=False)


    # plot and calculate regression for each Site-Variable combination (where multiple grades exist in OMR)
    correlations = []

    wide[["Value_CBM", "Grade_OMR"]] = wide[["Value_CBM", "Grade_OMR"]].apply(pd.to_numeric, errors='coerce')
    grouped_by = wide.groupby(['Variable', 'Site'])
    for (variable, site), df in grouped_by:
        var_site = f'{variable} in {site}'
        # df["Value_CBM"] = pd.to_numeric(df["Value_CBM"], errors="coerce")

        # verify regression can be calculated
        max_percent = df['Value_CBM'].max()
        if max_percent == 0:
            print(f'No positive measurements for {var_site}')
            continue
        if df["Grade_OMR"].nunique() == 1:
            print(f'No OMR grades assigned for {var_site}')
            continue

        # correlation measures
        rho, p_rho = spearmanr(df['Value_CBM'], df['Grade_OMR'])
        tau, p_tau = kendalltau(df['Value_CBM'], df['Grade_OMR'])

        print(f'{var_site}')
        print(f"Spearman’s ρ = {rho:.3f}, p = {p_rho:.3g}")
        print(f"Kendall’s τ = {tau:.3f}, p = {p_tau:.3g}")


        model = OrderedModel(df['Grade_OMR'], df[['Value_CBM']], distr='logit')
        res = model.fit(method='bfgs')


        
        # Plot regression (smooth curves)
        x_grid = np.linspace(0, max_percent + 1, 200)
        pred_grid = res.model.predict(res.params, exog=pd.DataFrame({'Value_CBM': x_grid}))

        plt.figure(figsize=(8, 5))
        for k in range(pred_grid.shape[1]):
            plt.plot(x_grid, pred_grid[:, k], label=f'P(grade={k})')


        # Jittered observed points
        pred_at_obs = res.model.predict(res.params, exog=df[['Value_CBM']])
        p_obs_grade = pred_at_obs[np.arange(len(df)), df['Grade_OMR'].to_numpy()]
        y_jitter = p_obs_grade + 0.01 * np.random.randn(len(df))
        plt.scatter(df['Value_CBM'], y_jitter, s=15, alpha=0.6, c=df['Grade_OMR'], cmap='viridis')

        plt.xlabel('CBM')
        plt.ylabel('Probability')
        plt.legend(loc='best')
        plt.title(var_site)
        plt.tight_layout()
        plt.show()
        




        # save correlations and number of positive examples
        grd_count_sum = summary.query(f"Site=='{site}' and Variable=='{variable}'")
        cors_dict = dict(zip(grd_count_sum['Grade_OMR'], grd_count_sum['size']))
        cors_dict.update({'Site': site, 'Variable': variable, 'Rho': rho, 'Tau': tau})

        correlations.append(cors_dict)

    cors = pd.DataFrame(correlations)
    cors.to_csv('omr_grades_correlations.csv', index=False)


    print(wide)
