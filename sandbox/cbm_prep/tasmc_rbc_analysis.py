from objects import MethodComparator
from sandbox import MetadataBundle
import sandbox as sb
import os

if __name__ == "__main__":
    site = 'TASMC'
    meta_path = r'config.yaml'
    # save_name = f'{site}_rbc'
    save_name = f'both_rbc'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)

    # srcs = {(site, 'CBM'): f'{site}_CBM.csv',
    #         (site, 'ClV'): f'{site}_ClV%.csv'}

    srcs = {('TASMC', 'CBM'): f'TASMC_CBM.csv',
            ('TASMC', 'ClV'): f'TASMC_ClV%.csv',
            ('BWH', 'CBM'): f'BWH_CBM.csv',
            ('BWH', 'ClV'): f'BWH_ClV%_means_decimal.csv'}

    # srcs = {(site, 'CBM'): f'{site}_CBM.csv'}

    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=r'raw/cbm_method_comparison')
    vars_to_test = metadata.variable_groups['RBC morphology']
    # methd_comp.batch_fit(['ClV'], ['CBM'], vars_to_test)
    # methd_comp.batch_fit(['ClV'], ['CBM'], vars_to_test, site_filters=['TASMC', 'BWH'])
    #
    # methd_comp.save_results(rf'results/{save_name}_reg.csv')
    # methd_comp.plot_all_regressions(f'results/{save_name}_reg.pdf')

    # comp_mtrx = methd_comp.export_comparison_matrix(out_path=f'comp_tables/{save_name}.csv',
    #                                                 comparison_dims=("Variable", "Method"),
    #                                                 needed_vars=vars_to_test)

    comp_mtrx = methd_comp.export_comparison_matrix(out_path=f'comp_tables/{save_name}.csv',
                                                    comparison_dims=("Variable", "Method"),
                                                    needed_vals=vars_to_test,
                                                    needed_grades=['scan_id'])

