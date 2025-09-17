import pandas as pd
import os
from objects import MethodComparator
from sandbox import MetadataBundle

if __name__ == "__main__":
    meta_path = r'config.yaml'
    site = 'LMU'
    save_name = 'LMU_Alina'
    analysis_name = "cbm_method_comparison"
    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    srcs = {(site, 'manual'): f'{site}_manual_Alina.csv',
            (site, 'CBM'): f'{site}_CBM.csv'}


    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=raw_dir, bma=False)
    vars_to_test = metadata.variable_groups['WBC&PLT compare']

    methd_comp.export_comparison_matrix(out_path=fr'comp_tables/omr_{save_name}.csv',
                                        row_identifiers=["SampleID"],
                                        comparison_dims=("Variable", "Method"),
                                        needed_vars=vars_to_test)

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=site)
    methd_comp.save_results(rf'results/{save_name}_all_scans_reg.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_all_scans_reg.pdf')


