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
    analysis_name = "bwh_version_compare"

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)
    metadata = sb.MetadataBundle('config.yaml')

    srcs = {(site, 'OMR'): f'OMR.csv',
            (site, 'old CBM'): f'old_CBM.csv',
            (site, 'new CBM'): f'new_CBM.csv'}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=raw_dir, bma=False)
    vars_to_test = metadata.variable_groups['OMR compare']

    methd_comp.export_comparison_matrix(out_path=fr'{analysis_name}/mthd_comp.csv',
                                        row_identifiers=("SampleID", "Variable"))

    methd_comp.batch_fit(['OMR'], ['old CBM', 'new CBM'], vars_to_test)
    methd_comp.save_results(rf'results/{analysis_name}_reg.csv')
    methd_comp.plot_all_regressions(f'results/{analysis_name}_reg.pdf')
