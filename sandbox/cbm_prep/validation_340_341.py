import pandas as pd
import sandbox as sb
import os
import sys
from objects import MethodComparator
import plotting
sys.path.append(r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification')
import reg_types as reg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import close

def short_pipe(df, metadata):
    # standardize column names
    df = sb.curate_df(df, metadata)

    if df.Site[0] == 'CPG' and df.Method[0] == 'OMR':
        df = sb.fill_nans(df, metadata, 'WBC diff', fill_value=0)
        df = sb.fill_nans(df, metadata, 'WBC-like', fill_value=0)

    df = sb.pivot_long(df)
    df = sb.add_grade_column(df, metadata)  # in future, include that step after concatenation
    df = df.dropna(subset=["Value", "Grade"], how='all')

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = sb.create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    return df


def medium_pipe(file_name, site, method, metadata, sheet_name='Sheet1', dir=None):
    df = sb.raw_to_df(file_name, site, method, sheet_name, dir)
    return short_pipe(df, metadata)


if __name__ == "__main__":
    site = 'Scopio-new-test-set'
    analysis_name = "August_validations"
    ref_mthd = '3.4.0'
    test_mthd = '3.4.1'
    save_name = f'{analysis_name}_{site}'

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    ref_df = medium_pipe(f'sb1438-3.4.0.csv', site, ref_mthd, metadata, dir=raw_dir)
    test_df = medium_pipe(f'sb3201-3.4.1.csv', site, test_mthd, metadata, dir=raw_dir)

    site_df = pd.concat([ref_df, test_df])
    methd_comp = MethodComparator(site_df)

    methd_comp.save_results(rf'results/{save_name}_reg.csv')
    methd_comp.save_results(rf'results/{save_name}_bias.xlsx', result_type='bias')

    vars_to_test = metadata.variable_groups['percent'] + metadata.variable_groups['percent-like'] + metadata.variable_groups['count']
    methd_comp.batch_fit(ref_mthd, test_mthd, vars_to_test)
    methd_comp.calc_all_biases(metadata.crit_points)

    figures = []
    for data, rslt in methd_comp.results.items():
        xname, yname, varname, _ = data
        style = {'title': varname,
                 'xlabel': xname,
                 'ylabel': yname,
                 'equal_limits': True,
                 'ci': True,
                 'ci_color': 'r',
                 'ci_mode': 'shade'}
        fig, ax = plotting.plot_scatter_basic(rslt, style=style)
        # to do: when creating plotting tools, add optional descriptor text for regression next to plot
        plotting.overlay_regression_line(fig=fig, ax=ax, result=rslt['reg'], style=style)

    with PdfPages(fr"results/{save_name}.pdf") as pdf:
        for data, rslt in methd_comp.results.items():
            if rslt['reg'].r is not None:
                xname, yname, varname, _ = data
                style = {'title': varname,
                         'xlabel': xname,
                         'ylabel': yname,
                         'equal_limits': True,
                         'ci': True,
                         'ci_color': 'r',
                         'ci_mode': 'shade'}
                fig, ax = plotting.plot_scatter_basic(rslt, style=style)
                plotting.overlay_regression_line(fig=fig, ax=ax, result=rslt['reg'], style=style)
                pdf.savefig(fig)
                close()
