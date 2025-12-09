import pandas as pd
import os
import sys
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df
from pipelines import medium_pipe
from itertools import product


if __name__ == "__main__":
    meta_path = r'config.yaml'
    save_name = 'all_site_omr'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    sptr_cbm_csv = False
    cases_to_filter = r'slides_to_remove_long.csv'

    metadata = MetadataBundle(meta_path)

    # build list of files to read from
    sites = ['BWH', 'CPG', 'LMU', 'SYN', 'TASMC']
    mthds = ['OMR', 'CBM']
    if sptr_cbm_csv:
        srcs = {(site, mthd): f'{site}_{mthd}.csv' for site, mthd in product(sites, mthds)}
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    else:
        # most cbm numbers are already in single file
        cbm_file_name = '5sites_CBM.csv'
        cbm_df = medium_pipe(cbm_file_name, None, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
        # syn_cbm_df = medium_pipe('SYN_CBM.csv', 'SYN', 'CBM', metadata, dir=r'raw/cbm_method_comparison')
        # cbm_df = pd.concat([cbm_df, syn_cbm_df])
        # gather omr as usual
        srcs = {(site, 'OMR'): f'{site}_OMR.csv' for site in sites}
        methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
        omr_df = methd_comp.df
        # combine
        df = pd.concat([cbm_df, omr_df])
        methd_comp = MethodComparator(df)

    vars_to_write = metadata.variable_groups['WBC&PLT compare'] + ['Total WBC', 'PLT', 'WBC', 'Parasites', 'Unclassified WBC']

    comp_mtrx = methd_comp.export_comparison_matrix(out_path=f'comp_tables/{save_name}.csv',
                                                    comparison_dims=("Variable", "Method"),
                                                    needed_vals=vars_to_write,
                                                    needed_grades=["scan_id"])

    filt_df = pd.read_csv(cases_to_filter)
    methd_comp.filter_by_df(filt_df)

    vars_to_test = metadata.variable_groups['WBC&PLT compare']


    methd_comp.batch_fit(['OMR'], ['CBM'], vars_to_test)
    methd_comp.batch_fit(['OMR'], ['CBM'], vars_to_test, site_filters=sites)
    methd_comp.calc_all_biases(metadata.crit_points)

    methd_comp.save_results(rf'results/{save_name}_reg.csv')
    methd_comp.save_results(rf'results/{save_name}_bias.xlsx', result_type='bias')
    methd_comp.plot_all_regressions(f'results/{save_name}_reg.pdf')


