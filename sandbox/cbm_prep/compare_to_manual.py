import pandas as pd
import os
from objects import MethodComparator
from sandbox import MetadataBundle
from pipelines import mean_manual_pipe, medium_pipe

if __name__ == "__main__":
    suffix = ''
    meta_path = r'config.yaml'
    site = 'TASMC'
    save_name = f'{site}_manual_both'
    analysis_name = "cbm_method_comparison"
    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    srcs = {(site, 'manual'): f'{site}_manual.csv',
            (site, 'CBM'): f'{site}_CBM.csv'}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=raw_dir, bma=False)

    # sites_investigators = {'LMU': ['Alina', 'Sladana'], 'CPG': ['Aubrey B Charlton', 'Deborah Swearingen']}
    # investigators = sites_investigators[site]

    df_srcs_list = []

    # import and pre-process each method separately
    df = mean_manual_pipe(f'{site}_manual.csv', site, metadata, dir=r'raw/cbm_method_comparison', only_mean=False)
    df_srcs_list.append(df)
    df = medium_pipe(f'{site}_CBM.csv', site, 'CBM', metadata, dir=r'raw/cbm_method_comparison')
    df["Investigator"] = "CBM"
    df_srcs_list.append(df)

    all_dfs = pd.concat(df_srcs_list)
    methd_comp_all_inv = MethodComparator(all_dfs)
    methd_comp = methd_comp_all_inv.apply_to_df('query', "Investigator=='Mean Investigator' or Investigator=='CBM'", inplace=False)



    vars_to_test = metadata.variable_groups['WBC&PLT compare']
    grades_to_test = ['scan_id'] + metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']
    morph_vals_to_test = metadata.variable_groups['WBC morphology'] + metadata.variable_groups['PLT morphology']

    methd_comp.export_comparison_matrix(
        out_path=fr'comp_tables/{save_name}_vals.csv',
        row_identifiers=["SampleID"],
        comparison_dims=("Variable", "Method", "Investigator"),
        needed_vals=vars_to_test,
        needed_grades=['scan_id'])
    #
    #
    # methd_comp_all_inv.export_comparison_matrix(
    #     out_path=fr'comp_tables/{save_name}_vals_all_inv.csv',
    #     row_identifiers=["SampleID"],
    #     comparison_dims=("Variable", "Method", "Investigator"),
    #     needed_grades=vars_to_test)

    methd_comp_all_inv.export_comparison_matrix(
        out_path=fr'comp_tables/{save_name}_grades_all_inv.csv',
        row_identifiers=["SampleID"],
        comparison_dims=("Variable", "Method", "Investigator"),
        needed_vals=[],
        needed_grades=grades_to_test)

    # methd_comp_all_inv.export_comparison_matrix(
    #     out_path=fr'comp_tables/{save_name}_morphs_all_inv.csv',
    #     row_identifiers=["SampleID"],
    #     comparison_dims=("Variable", "Method", "Investigator"),
    #     needed_vars=morph_vals_to_test)

    # methd_comp.export_comparison_matrix(out_path=fr'comp_tables/{save_name}.csv',
    #                                     row_identifiers=["SampleID"],
    #                                     comparison_dims=("Variable", "Method"),
    #                                     needed_vars=vars_to_test)

    methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=site)
    methd_comp.save_results(rf'results/{save_name}_all_scans_reg.csv')
    methd_comp.plot_all_regressions(f'results/{save_name}_all_scans_reg.pdf')

    # for inv in investigators:
    #     inv_methd_comp = methd_comp_all_inv.apply_to_df(f'query', f"Investigator=='{inv}' or Investigator=='CBM'", inplace=False)
    #     inv_methd_comp.batch_fit(['manual'], ['CBM'], vars_to_test, site_filters=site)
    #     inv_methd_comp.save_results(rf'results/{save_name}_all_scans_reg_{inv}.csv')
    #     inv_methd_comp.plot_all_regressions(f'results/{save_name}_all_scans_reg_{inv}.pdf')


