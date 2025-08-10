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
    df = sb.stnd_names(df, metadata.alias_map)

    # print warning if WBCs in differential don't add up to ~100
    sb.check_diff_sum(df, metadata, tolerance=5)

    if df.Site[0] == 'CPG' and df.Method[0] == 'OMR':
        df = sb.fill_nans(df, metadata, 'WBC diff', fill_value=0)
        df = sb.fill_nans(df, metadata, 'WBC-like', fill_value=0)

    df = sb.pivot_long(df)
    df = df.dropna(subset=["Value"])

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = sb.create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    return df


def medium_pipe(file_name, site, method, metadata, sheet_name='Sheet1', dir=None):
    df = sb.raw_to_df(file_name, site, method, sheet_name, dir)
    return short_pipe(df, metadata)


if __name__ == "__main__":
    site = 'CPG'
    analysis_name = "new sites"

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    omr_df = medium_pipe(f'{site}_OMRs.csv', site, 'OMR', metadata, dir=raw_dir)
    mnl_df = medium_pipe(f'{site}_DSS.xlsx', site, 'DSS', metadata, dir=raw_dir)

    site_df = pd.concat([omr_df, mnl_df])
    site_df = site_df[site_df["Value"].notna()]  # remember to include in larger pipeline later
    methd_comp = MethodComparator(site_df)
    methd_comp.batch_fit(['OMR'], ['DSS'], metadata.variable_groups['percent'])
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.save_results(rf'results/{site}_reg.csv')
    methd_comp.save_results(rf'results/{site}_bias.xlsx', result_type='bias')
