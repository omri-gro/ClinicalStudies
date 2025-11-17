import pandas as pd
import os
import sys
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from itertools import product


if __name__ == "__main__":
    meta_path = r'config.yaml'
    save_name = 'all_site_omr'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    # build list of files to read from
    sites = ['BWH', 'CPG', 'LMU', 'SYN', 'TASMC']
    mthds = ['OMR', 'CBM']
    srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    vars_to_write = metadata.variable_groups['WBC&PLT compare'] + ['Total WBC', 'PLT', 'WBC', 'Parasites', 'Unclassified WBC']


    comp_mtrx = methd_comp.export_comparison_matrix(out_path=f'comp_tables/{save_name}.csv',
                                                    comparison_dims=("Variable", "Method"),
                                                    needed_vals=vars_to_write,
                                                    needed_grades=["scan_id"])

    vars_to_test = metadata.variable_groups['WBC&PLT compare']


    methd_comp.batch_fit(['OMR'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.batch_fit(['OMR'], ['CBM'], vars_to_test)
    methd_comp.calc_all_biases(metadata.crit_points)

    methd_comp.save_results(rf'results/{save_name}_reg.csv')
    # methd_comp.save_results(rf'results/{save_name}_bias.xlsx', result_type='bias')
    methd_comp.plot_all_regressions(f'results/{save_name}_reg.pdf')


