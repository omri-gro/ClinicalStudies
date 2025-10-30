"""
Calculate performance of two different CBM versions when comparing to OMR, with both comparisons using the same samples.
Will also act as buildup for future function to calculate multiple regressions/stats only on samples that exist in multiple methods.
"""
import pandas as pd
import os
import sys
import sandbox as sb
from objects import MethodComparator

if __name__ == "__main__":
    site = 'BWH'
    # site = 'TASMC'
    analysis_name = "cbm_method_comparison"
    output_name = f"{site}_cbm_clv_rbc_comparison"
    output_name = f"{site}_cbm_clv_rbc_comparison"

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)
    metadata = sb.MetadataBundle('config.yaml')

    srcs = {(site, 'ClV'): f'{site}_ClV%_means_decimal.csv',
            (site, 'CBM'): f'{site}_CBM.csv'}

    # srcs = {(site, 'ClV'): f'BWH_ClV%_means_decimal.csv',
    #         # (site, 'ClV'): f'BWH_ClV%_means.csv',
    #         # (site, 'DSS'): f'BWH_DSS_old_dec.csv',
    #         # (site, 'DSS'): f'BWH_DSS_old.csv',
    #         (site, 'CBM'): f'BWH_CBM.csv'}

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=raw_dir, bma=False)
    vars_to_test = metadata.variable_groups['RBC morphology'] + metadata.variable_groups['PLT morphology']

    methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{output_name}.csv',
                                        row_identifiers=["SampleID"],
                                        comparison_dims=("Variable", "Method"),
                                        needed_vars=vars_to_test)

    methd_comp.batch_fit(['ClV'], ['CBM'], vars_to_test, site_filters=site)
    methd_comp.save_results(rf'results/{output_name}_all_scans_reg.csv')
    methd_comp.plot_all_regressions(f'results/{output_name}_all_scans_reg.pdf')
