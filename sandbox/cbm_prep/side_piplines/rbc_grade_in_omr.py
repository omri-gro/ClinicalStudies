import os
import sys
import pandas as pd, matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr, kendalltau
from statsmodels.miscmodels.ordinal_model import OrderedModel

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import medium_pipe
from itertools import *


def grades_by_rprt_frmts(df_long: pd.DataFrame, df_frmt: pd.DataFrame, meta: "MetadataBundle" = None, other: str = 'stay'):
    """
    Expected columns for df_frmt: Variable, Site, Reporting Format (nice to have to do: add option for Method?)
    Adds columns for grades and positivity based on reporting methods:
        'graded' - value is already grade - copy value to grade column and delete from value column
        'binary' - leave grade and value as Na but set positive/negative
        'no' - leave grade and value as Na
        'stay' - leave as is (will become more relevant once frmts by method are added)
    stay - instructions if any
    to do: consider adding verification of grades being possible for variable based on meta
    """
    # to do: add verification that formats in df_frmt are correct
    possible_frmts = ['graded', 'binary', 'no', 'stay']

    rpfrm = 'Reporting Format'
    val = 'Value'
    grade = 'Grade'
    positivity = 'Positive'

    df = df_long.merge(df_frmt, how='left', on=['Site', 'Variable'])
    df[rpfrm] = df[rpfrm].fillna(other)

    for col_name in [grade, positivity]:
        if col_name not in df.columns:
            df[col_name] = np.nan

    # if 'graded', grade is equal to provided value
    df[grade] = np.where(df[rpfrm] == 'graded', df[val], df[grade])

    # add positivity by positive value options if 'binary' or 'graded'
    pos_frmts = ['graded', 'binary']
    normal_grades = [0, "0", "Normal", "Negative", "normal", "negative"]
    df[positivity] = np.where(df[rpfrm].isin(pos_frmts),
                              ~df[val].isin(normal_grades),
                              df[positivity])
    # cases where value is any type of NaN should also have NaN in grade and positivity
    df[grade] = np.where(df[val].isna(), np.nan, df[grade])
    df[positivity] = np.where(df[val].isna(), np.nan, df[positivity])

    # remove values/grades/positives that are known not to be real
    no_val_frmts = ['graded', 'binary', 'no']
    no_grade_frmts = ['binary', 'no']
    df.loc[df[rpfrm].isin(no_val_frmts), val] = np.nan
    df.loc[df[rpfrm].isin(no_grade_frmts), grade] = np.nan
    df.loc[df[rpfrm] == 'no', positivity] = np.nan

    return df


def perc_to_grade_bplot(df, variable):
    fig = plt.figure(figsize=(8, 7))
    xname = f'{variable} OMR Grade'
    yname = f'{variable} CBM %'
    box = sns.boxplot(x='Grade_OMR', y='Value_CBM', data=df, width=0.4, showmeans=True,
                fill=False, fliersize=0, linecolor='grey',
                meanprops={"marker": "o", "markerfacecolor": "red", "markeredgecolor": "black"})
    for patch in box.artists:  # The boxes
        patch.set_alpha(0.5)  # Transparency level
    for line in box.lines:  # The whiskers, caps, and medians
        line.set_alpha(0.5)  # Transparency level
    sns.stripplot(x='Grade_OMR', y='Value_CBM', data=df,
                  color='black', size=3, alpha=0.9, dodge=True, jitter=True)
    plt.ylim(0, None)
    plt.xlabel(xname, fontweight='bold')
    plt.ylabel(yname, fontweight='bold')
    plt.title(f'{variable} in {site}')
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.7)
    return fig

if __name__ == "__main__":
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    meta_path = r'config.yaml'
    save_name = 'site_combinations_omr'
    pdf_path = f'{save_name}_bplots.pdf'

    metadata = MetadataBundle(meta_path)

    # all cbm numbers are already in single file
    cbm_file_name = '6sites_CBM.csv'
    cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    cbm_df = cbm_df[cbm_df['Value'].notna()]
    cbm_df_vals = cbm_df.drop(['Grade', 'Positive', 'Grade_from'], axis=1)

    # build df from omr files
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'TASMC']  # for now not dealing with SYN for simplicity
    srcs = {(site, 'OMR'): f'{site}_OMR.csv' for site in sites}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    omr_df = methd_comp.df

    # redo grade and positivity column using omr reporting formats file
    frmt_file = 'omr_rbc_reporting_formats.csv'
    frmt_df = read_to_df(frmt_file, file_dir=os.getcwd())
    frmt_df['Reporting Format'] = frmt_df['Reporting Format'].replace(metadata.alias_map)

    # just for now, treat the binary as grade
    frmt_df.loc[frmt_df['Reporting Format'] == 'binary', 'Reporting Format'] = 'graded'

    # omr_df = grades_by_rprt_frmts(omr_df, frmt_df)
    omr_df = omr_df[omr_df['Grade'].notna()]

    # create dataframe from both methods
    df = pd.concat([cbm_df_vals, omr_df])




    # filter out inappropriate scans
    # rmv_file = 'slides_to_remove.csv'
    # rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())

    graded_params = frmt_df['Variable']
    is_in_list = df['Variable'].isin(graded_params)
    grades_df = df[is_in_list]

    # get only samples that have both OMR and CBM results
    df = grades_df.groupby(['SampleID', "Site", "Variable"]).filter(lambda g: g["Method"].nunique() == 2)

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
    """
    pd.DataFrame(summary).to_csv('omr_grades_summary.csv', index=False)
    """
    val = "Value_CBM"
    grd = "Grade_OMR"


    # plot and calculate regression for each Site-Variable combination (where multiple grades exist in OMR)
    correlations = []

    wide[["Value_CBM", "Grade_OMR"]] = wide[["Value_CBM", "Grade_OMR"]].apply(pd.to_numeric, errors='coerce')


    # box plots
    grouped_by = wide.groupby(['Variable', 'Site'])
    with PdfPages(pdf_path) as pdf:
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

            # calculate statistical difference between grades
            groups = [group[val].values for name, group in df.groupby(grd)]
            f_stat, p_value = f_oneway(*groups)

            # calculate eta squared
            grand_mean = df[val].mean()
            ss_total = ((df[val] - grand_mean) ** 2).sum()
            ss_between = sum([
                len(group) * (group[val].mean() - grand_mean) ** 2
                for name, group in df.groupby(grd)
            ])
            eta_squared = ss_between / ss_total


            # save correlations and number of positive examples
            grd_count_sum = summary.query(f"Site=='{site}' and Variable=='{variable}'")
            cors_dict = dict(zip(grd_count_sum['Grade_OMR'], grd_count_sum['size']))
            cors_dict.update({'Site': site, 'Variable': variable, 'f': f_stat, 'p': p_value, 'eta^2': eta_squared})

            correlations.append(cors_dict)

            fig = perc_to_grade_bplot(df, variable)
            pdf.savefig(fig)
            plt.close()



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

        max_grade_index = pred_at_obs.shape[1] - 1
        grade_idx = np.clip(df['Grade_OMR'].astype(int), 0, max_grade_index)
        p_obs_grade = pred_at_obs[np.arange(len(df)), grade_idx]

        # p_obs_grade = pred_at_obs[np.arange(len(df)), df['Grade_OMR'].astype(int)]
        y_jitter = p_obs_grade + 0.01 * np.random.randn(len(df))
        plt.scatter(df['Value_CBM'], y_jitter, s=15, alpha=0.6, c=df['Grade_OMR'], cmap='viridis')

        plt.xlabel('CBM')
        plt.ylabel('Probability')
        plt.legend(loc='best')
        plt.title(var_site)
        plt.tight_layout()
        plt.show()


    cors = pd.DataFrame(correlations)
    cors.to_csv('omr_grades_correlations.csv', index=False)






