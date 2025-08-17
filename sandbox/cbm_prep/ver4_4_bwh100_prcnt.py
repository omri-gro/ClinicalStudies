import pandas as pd
import sandbox as sb
import os
import sys
sys.path.append(r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification')
import reg_types as reg
from matplotlib.backends.backend_pdf import PdfPages



def short_pipe(df):
    # standardize column names
    df = sb.stnd_names(df, metadata.alias_map)

    # print warning if WBCs in differential don't add up to ~100
    sb.check_diff_sum(df, metadata, tolerance=5)

    df = sb.pivot_long(df)
    df = df.dropna(subset=["Value"])

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = sb.create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    return df


def get_matched_arrays(df,
                       ref_method,
                       test_method,
                       variable,
                       site_filter=None,
                       measurement_col='Value'):
    """
    Filters for a given variable and two review methods across one or more sites.
    Matches and aligns samples by (Site, SampleID), ensuring within-site comparison.

    Parameters:
        df: DataFrame in long format
        ref_method: name of reference review method (e.g., 'OMR')
        test_method: name of test review method (e.g., 'DSS')
        variable: variable to filter on (e.g., 'Total Neutrophil')
        site_filter: list or set of sites to include (optional)
        measurement_col: column name holding the measurement values

    Returns:
        x: array of reference method measurements
        y: array of test method measurements
        sample_ids: list of (Site, SampleID) tuples corresponding to each x/y pair
    """
    # note this function actually does same thing as _prepare_arrays method from MethodComparator

    # Step 1: Filter by variable
    var_df = df[df['Variable'] == variable]

    # Step 2: Filter by site if requested
    if site_filter is not None:
        var_df = var_df[var_df['Site'].isin(site_filter)]

    # Step 3: Filter to the two methods
    method_df = var_df[var_df['Method'].isin([ref_method, test_method])]

    # Step 4: Create a MultiIndex for matching on both Site and SampleID
    method_df = method_df.set_index(['Site', 'SampleID'])

    # Step 5: Pivot so each (Site, SampleID) has one column per method
    pivot_df = (
        method_df
        .pivot(columns='Method', values=measurement_col)
        .dropna(subset=[ref_method, test_method])  # keep only matched pairs
    )

    # Step 6: Extract aligned arrays and keys
    x = pivot_df[ref_method].values
    y = pivot_df[test_method].values
    sample_ids = list(pivot_df.index)  # list of (Site, SampleID) tuples

    return x, y, sample_ids


if __name__ == "__main__":
    analysis_name = 'bwh100_investigator_training'
    suffix = '_old_version'
    # suffix = ''
    save_name = f"{analysis_name}{suffix}"
    site = 'BWH'

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    clv_df_raw = sb.raw_to_df('BWH_ClV%.csv', site, 'ClV', dir=raw_dir)
    mnl_df_raw = sb.raw_to_df('BWH_manual%.csv', site, 'Manual', dir=raw_dir)
    cbm_df_raw = sb.raw_to_df(f'BWH_CBM%{suffix}.csv', site, 'CBM', dir=raw_dir)

    clv_df = short_pipe(clv_df_raw)
    mnl_df = short_pipe(mnl_df_raw)
    cbm_df = short_pipe(cbm_df_raw)

    bwh_df = pd.concat([clv_df, mnl_df, cbm_df])

    bwh_df = bwh_df[bwh_df["Value"].notna()]  # remember to include in larger pipeline later

    reg_results = []
    figs = []

    cell_to_check = metadata.variable_groups["WBC diff"] + metadata.variable_groups["WBC-like"]
    for cell in cell_to_check:
        try:
            x, y, _ = get_matched_arrays(bwh_df, 'Manual', 'CBM', cell)
            dem_dict = reg.regression_comp(x, y)
            dem_dict["Parameter"] = cell
            reg_results.append(dem_dict)
            print(dem_dict)
            if cell in ['Segmented Neutrophil', 'Band Neutrophil', 'Metamyelocyte', 'Myelocyte', 'Promyelocyte', 'Blast', 'Lymphocyte (not regular)', 'LGL', 'Atypical Lymphocyte', 'Aberrant Lymphocyte', 'Monocyte', 'Basophil', 'Eosinophil', 'Total Neutrophil', 'Variant Lymphocyte', 'Total Lymphocyte', 'Immature', 'Immature&Blast', 'NRBC', 'Smudge Cell']:
                fig, ax = sb.plot_ver_reg(x, y, reg_ser=dem_dict)
                fig.show()
                sb.set_equal_limits_and_scale(ax=ax)
                ax.set_title(cell)
                ax.set_xlabel('Manual')
                ax.set_ylabel('Scopio')
                figs.append(fig)
        except KeyError as e:
            print(e)

    cell_to_check = metadata.variable_groups["RBC morphology"]
    for cell in cell_to_check:
        try:
            x, y, _ = get_matched_arrays(bwh_df, 'ClV', 'CBM', cell)
            if suffix == '':
                y = y * 100
            x = x * 100
            dem_dict = reg.regression_comp(x, y)
            dem_dict["Parameter"] = cell
            reg_results.append(dem_dict)
            print(dem_dict)
            if cell in ['Burr Cells', 'Elliptocytes', 'Ovalocytes', 'Pappenheimer', 'Parasites', 'Schistocytes', 'Sickle cells', 'Spherocytes', 'Stomatocytes', 'Target cells', 'Tear drop cells', 'EllipOval', 'AcanoEchino']:
                fig, ax = sb.plot_ver_reg(x, y, reg_ser=dem_dict)
                fig.show()
                sb.set_equal_limits_and_scale(ax=ax)
                ax.set_title(cell)
                ax.set_xlabel('CellaVision')
                ax.set_ylabel('Scopio')
                figs.append(fig)
        except KeyError as e:
            print(e)

    final_df = pd.DataFrame(reg_results)
    final_df.to_excel(fr"results/{save_name}.xlsx", index=False)

    with PdfPages(fr"results/{save_name}.pdf") as pdf:
        for fig in figs:
            pdf.savefig(fig)

    print('')


