import pandas as pd
import sandbox as sb
import os
import sys
from objects import MethodComparator
sys.path.append(r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification')
import reg_types as reg
from matplotlib.backends.backend_pdf import PdfPages


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
    site = 'CPG'
    analysis_name = "new_sites"
    save_name = "256_scans"
    ref_arm = 'OMR'
    test_arm = 'CBM'

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    omr_df = medium_pipe(f'{site}_OMRs.csv', site, ref_arm, metadata, dir=raw_dir)
    mnl_df = medium_pipe(f'{site}_CBM.csv', site, test_arm, metadata, dir=raw_dir)

    site_df = pd.concat([omr_df, mnl_df])
    methd_comp = MethodComparator(site_df)

    vars_to_test = metadata.variable_groups['percent'] + metadata.variable_groups['percent-like'] + metadata.variable_groups['count']

    methd_comp.batch_fit([ref_arm], [test_arm], vars_to_test)
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.save_results(rf'results/{site}_{save_name}_reg.csv')
    methd_comp.save_results(rf'results/{site}_{save_name}_bias.xlsx', result_type='bias')
    methd_comp.plot_all_regressions(f'results/{site}_{save_name}_reg.pdf')

