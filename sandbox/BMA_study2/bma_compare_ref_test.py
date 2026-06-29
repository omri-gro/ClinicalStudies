import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle, raw_bma_to_df, _ensure_list
from sandbox import *
from pipelines import bma_prep_pipeline
# trgt_dict = os.path.abspath(r'../../clinstudtools')
# sys.path.append(trgt_dict)
# from table_integrity import robust_dup
from bma_specific_functs import removed_for_arbitration,generate_fda_equivocal_report


if __name__ == "__main__":
    """ currently added change where arbitrator is counted as just another investigator for calculation of mean 
        return commented section (remove for arbitration + add mean investigator) in loop and remove it from for all_df """

    suffix = ''
    save_name = f'BMA_study_results{suffix}'
    meta_path = r'config_BMA.yaml'
    sites = ["OHSU", "HUP", "BWH"]
    arbitrators = ['Phil Raess', 'Olga Pozdnyakova', 'Christopher Hergott', 'OP', 'Arbitrator']

    compare_methods = True
    raw_dss = False
    inter = False
    inter_to_include_arbitrated = True  # also means that reviews which were replaced by arbitration will not appear in comp_mtrx

    # comp_mk_rois = True

    exprt_mtrx = False
    exprt_long = False
    plot_reg = False
    keep_names = False  # use investigators' full names - creates very wide 'all investigators' comparison matrix if True
    min_inv_site = 2

    other_removed = True   # for filtering out samples for side analysis
    rmv_unclass = False    # re-calculate differential as if all unclassified were moved to dirt/other
    pooled_params = True   # analyze for pooled parameters like Erythroblast&BasophilicNormoblast - only when rmv_unclass False
    add_filtering = True   # analysis of promyelocytes without APLs, analysis of erythroid stages without HUP
    only_merged = False


    investigators_map = {'Todd Williams': 'Rev1', 'Wei Xie': 'Rev2', 'Phil Raess': 'Arbitrator',
                         'TW': 'Rev1', 'WX': 'Rev2', 'PR': 'Arbitrator',
                         'AB': 'Rev1', 'AS': 'Rev2', 'DL': 'Rev3', 'OP': 'Arbitrator',
                         'Adam Bagg': 'Rev1', 'Annapurna Saksena': 'Rev2', 'Dorottya Laczko': 'Rev3', 'Olga Pozdnyakova': 'Arbitrator',
                         'Elizabeth Morgan': 'Rev1', 'Habibe Kurt': 'Rev2', 'Robert Hasserjian': 'Rev3',
                         'Sam Sadigh': 'Rev4', "Megan Fitzpatrick": 'Rev5', "Vignesh Shanmugam": 'Rev6',
                         "Christopher Hergott-rev": 'Rev7', "Christopher Hergott": 'Arbitrator',
                         'Rev1': 'Rev1', 'Rev2': 'Rev2',
                         'Mean Investigator': 'Mean Investigator'}
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    read_dir = os.path.join(cur_dir, 'raw')
    df_map = read_to_df(f'BMA_mapping.csv', file_dir=read_dir)

    if rmv_unclass:
        meta_path = r'config_BMA_no_unclass.yaml'
    elif pooled_params:
        meta_path = r'config_BMA_pool.yaml'

    metadata = MetadataBundle(meta_path)
    collect_dfs = []
    for site in sites:
        ref_df = bma_prep_pipeline(f'{site}_CRF_REF.csv', site, 'REF', metadata, dir=read_dir, recalc_diff=rmv_unclass)
        test_df = bma_prep_pipeline(f'{site}_CRF_TEST.csv', site, 'TEST', metadata, dir=read_dir)

        ref_df = min_inv_filt(ref_df, 'REF', min_inv=min_inv_site)
        test_df = min_inv_filt(test_df, 'TEST', min_inv=min_inv_site)

        # df_arb = df_arb[~((df_arb['Site'] == 'HUP') & (df_arb['Method'] == 'TEST'))]

        """
        ref_df = removed_for_arbitration(ref_df, df_arb, arbitrators)
        test_df = removed_for_arbitration(test_df, df_arb, arbitrators)

        ref_df = add_mean_investigator(ref_df, mthd='REF', min_inv=0)
        test_df = add_mean_investigator(test_df, mthd='TEST', min_inv=0)
        """

        # change all SampleIDs to the TEST Barcode based on site's mapping
        id_lookup = df_map.set_index('REF Barcode')['TEST Barcode']
        mapped_ids = ref_df['SampleID'].map(id_lookup)
        ref_df['SampleID'] = mapped_ids.fillna(ref_df['SampleID'])

        collect_dfs.append(pd.concat([ref_df, test_df]))

    all_dfs = pd.concat(collect_dfs)
    if not keep_names:
        all_dfs['Investigator'] = all_dfs['Investigator'].map(investigators_map)
    all_dfs_all_inv = all_dfs
    # all_dfs_all_inv = add_mean_investigator(all_dfs_all_inv, mthd='REF', min_inv=0)
    # all_dfs_all_inv = add_mean_investigator(all_dfs_all_inv, mthd='TEST', min_inv=0)
    # if raw_dss:
    #     all_dfs_all_inv = add_mean_investigator(all_dfs_all_inv, mthd='DSS', min_inv=0)
    all_dfs_all_inv = add_pos_column(all_dfs_all_inv, metadata)

    df_arb = read_to_df('to_arbitration.csv', file_dir=read_dir)

    if not inter_to_include_arbitrated:
        all_dfs_all_inv = removed_for_arbitration(all_dfs_all_inv, df_arb, arbitrators)
    methd_comp_all_inv = MethodComparator(all_dfs_all_inv)

    all_dfs = removed_for_arbitration(all_dfs, df_arb, arbitrators)
    # all_dfs = all_dfs[~(all_dfs['Investigator'].isin(arbitrators))]  # when not using arbitration at all

    all_dfs = add_mean_investigator(all_dfs, mthd='REF', min_inv=0)
    all_dfs = add_mean_investigator(all_dfs, mthd='TEST', min_inv=0)
    if raw_dss:
        df_dss = read_to_df(f'raw_DSS.csv', file_dir=read_dir)
        df_dss = stnd_names(df_dss, metadata.alias_map)
        df_dss["Method"] = 'DSS'
        df_dss["FileName"] = os.path.basename(f'raw_DSS.csv')
        df_dss.columns = df_dss.columns.str.strip()

        # for raw dss we are not including unclassified in differential calculation
        df_dss = df_dss.drop('Unclassified', axis=1)

        df_dss = calc_diff(df_dss, metadata, diff_cells="NDC")
        id_vars = ["SampleID", "Site", "Method", "FileName", 'Investigator']
        df_dss = pivot_long(df_dss, id_vars=id_vars)
        df_dss = add_grade_column(df_dss, metadata)
        df_dss = df_dss.dropna(subset=["Value", "Grade"], how="all")
        df_dss = create_derived_variables_long(df_dss, metadata)

        # uncomment this section to only show mean of all scans in the wide comparison matrix
        df_dss = add_mean_investigator(df_dss, mthd='DSS', min_inv=0)
        df_dss = df_dss.query("Investigator=='Mean Investigator'", inplace=False)
        all_dfs = pd.concat([all_dfs, df_dss])

    all_dfs = add_pos_column(all_dfs, metadata)

    methd_comp = MethodComparator(all_dfs)
    methd_comp = methd_comp.apply_to_df('query', "Investigator=='Mean Investigator'", inplace=False)

    rmv_file = 'flt_lists/slides_to_remove.csv'
    rmv_df = read_to_df(rmv_file, file_dir=os.getcwd())
    methd_comp = methd_comp.filter_by_df(rmv_df)

    if only_merged:
        merged_cases_path = 'flt_lists/merged_cases.csv'
        merged_cases_df = read_to_df(merged_cases_path, file_dir=os.getcwd())
        methd_comp = methd_comp.filter_by_df(merged_cases_df, include_rows=True)
        save_name = f'{save_name}_merged_cases'

    if exprt_long:
        long_df_to_exprt = all_dfs   # add here query if needed
        long_df_to_exprt.to_csv(fr'{cur_dir}/comp_tables/{save_name}_long.csv', index=False)

    ndc_vars_list = metadata.variable_groups['NDC'] + metadata.variable_groups['NDC lineage total']
    ndc_vars_list_to_print = ndc_vars_list + ['Total Nucleated']
    if not rmv_unclass:
        ndc_vars_list.remove('Unclassified')
    needed_rows = methd_comp.df[methd_comp.df['Variable'].isin(ndc_vars_list)]
    grade_vars_list = metadata.variable_groups['grade']
    ids_list = ['scan_id', 'case_id']


    if exprt_mtrx:
        comp_table = methd_comp.export_comparison_matrix(needed_vars=ndc_vars_list_to_print,
                                                         # needed_grades=ids_list,
                                                         comparison_dims=("Variable", "Method"),
                                                         row_completeness="none")
        comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
        id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
        comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
        comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
        comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}.csv', index=True)

        comp_table = methd_comp_all_inv.export_comparison_matrix(needed_vals=ndc_vars_list_to_print,
                                                                 needed_grades=grade_vars_list + ids_list,
                                                                 comparison_dims=("Variable", "Method", "Investigator"),
                                                                 row_completeness="none")
        comp_table.rename(columns={'SampleID': 'TEST Barcode'}, inplace=True)
        id_lookup = df_map.set_index('TEST Barcode')['REF Barcode']
        comp_table['REF Barcode'] = comp_table['TEST Barcode'].map(id_lookup)
        comp_table = comp_table.set_index(['TEST Barcode', 'REF Barcode', 'Site'])
        comp_table.to_csv(fr'{cur_dir}/comp_tables/{save_name}_all_investigators.csv', index=True)


    if add_filtering:
        # analysis of promyelocytes is without APL slides
        apls_file = 'flt_lists/APLs.csv'
        apls_df = read_to_df(apls_file, file_dir=os.getcwd())
        apls_df['Variable'] = 'Promyelocyte'
        methd_comp = methd_comp.filter_by_df(apls_df)

        # the multi-site analysis for erythroid stages does not cover HUP
        erythroid_stages = ['Erythroblast', 'Basophilic normoblast', 'Polychromatophilic normoblast', 'Normoblast', 'Myelocyte', 'Metamyelocyte']
        methd_comp = methd_comp.apply_to_df('query', "Variable not in @erythroid_stages or Site != 'HUP'", inplace=False)


    # main method comparison regressions + biases
    if compare_methods:
        methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list)
        # if not only_merged:
        #     methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list, site_filters=sites)
        methd_comp.calc_all_biases(metadata.crit_points)
        methd_comp.save_results(rf'results/{save_name}_bma_reg.csv')
        methd_comp.save_results(rf'results/{save_name}_bias.xlsx', result_type='bias')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/{save_name}_bma_reg.pdf')

        # compare using sensitivity/specificity
        methd_comp.batch_compare(levels_a='REF', levels_b='TEST', variables=ndc_vars_list, comp_func='binary')
        methd_comp.batch_compare(levels_a='REF', levels_b='TEST', variables=ndc_vars_list, split_by='Site', comp_func='binary')
        methd_comp.save_results(rf'results/{save_name}_bin.csv', result_type="binary")

        # FDA Medical Decision Level (MDL) Equivocal Report
        # Generates the strict and equivocal-adjusted metrics using the
        # critical points already defined in your metadata config.
        print("Generating FDA Equivocal Zone MDL Report...")
        fda_report_df = generate_fda_equivocal_report(
            methd_comp=methd_comp,
            decision_dict=metadata.crit_points,
            level_a='REF',
            level_b='TEST',
            dim_col='Method'
        )
        fda_report_path = rf'results/{save_name}_FDA_MDL_Equivocal_Report.csv'
        fda_report_df.to_csv(fda_report_path, index=False)
        print(f"Exported FDA MDL Equivocal Report to: {fda_report_path}")

    # inter-investigator
    if inter:
        if keep_names:
            levels_a = ['Todd Williams', 'AB', 'AS', 'Adam Bagg', 'Annapurna Saksena',
                        'Elizabeth Morgan', 'Habibe Kurt', 'Robert Hasserjian']
            levels_b = ['Wei Xie', 'AS', 'DL', 'Adam Bagg', 'Dorottya Laczko',
                        'Habibe Kurt', 'Robert Hasserjian', 'Sam Sadigh']
        else:
            levels_a = ['Rev1']
            levels_b = ['Rev2', 'Rev3']
        methd_comp_all_inv.batch_compare(levels_a=levels_a, levels_b=levels_b, variables=ndc_vars_list,
                                 dim_col='Investigator', split_by=['Method', 'Site'])
        methd_comp_all_inv.save_results(rf'results/{save_name}_inter.csv')
        if plot_reg:
            methd_comp_all_inv.plot_all_regressions(f'results/{save_name}_inter.pdf')


    # raw dss vs reference regression
    if raw_dss:
        # methd_comp.clean_calculations()
        methd_comp.batch_compare(levels_a=['REF', 'TEST'], levels_b='DSS', variables=ndc_vars_list)
        methd_comp.batch_compare(levels_a=['REF', 'TEST'], levels_b='DSS', variables=ndc_vars_list, split_by='Site')
        methd_comp.save_results(rf'results/{save_name}_raw_dss.csv')
        if plot_reg:
            methd_comp.plot_all_regressions(f'results/{save_name}_raw_dss.pdf')

    # ---------------------------------------------------------
    # Generate MS Word Appendix
    # ---------------------------------------------------------
    from appendix_generator import create_word_appendix

    appendix_order = [
        'Total Myeloid',
        'Total Erythroid',
        'Blast',
        'Promyelocyte',
        'MyeloMetamyelo',
        'Myelocyte',
        'Metamyelocyte',
        'Neutrophil',
        'Band neutrophil',
        'Segmented neutrophil',
        'Plasma cell',
        'Lymphocyte',
        'Erythroblast&BasophilicNormoblast',
        'BasNorm&PolychromNorm',
        'PolychromNorm&Normoblast',
        'Erythroblast',
        'Basophilic normoblast',
        'Polychromatophilic normoblast',
        'Normoblast',
        'Monocyte',
        'Eosinophil',
        'Basophil',
        'Mast cell'
    ]

    fig_titles = {
        'Total Myeloid': 'Total Myeloid Lineage',
        'Total Erythroid': 'Total Erythroid Lineage',
        'Blast': 'Blast',
        'Promyelocyte': 'Promyelocyte',
        'MyeloMetamyelo': 'Intermediate Myeloid Precursors',
        'Myelocyte': 'Myelocyte',
        'Metamyelocyte': 'Metamyelocyte',
        'Neutrophil': 'Total neutrophils',
        'Band neutrophil': 'Band neutrophil',
        'Segmented neutrophil': 'Segmented neutrophil',
        'Plasma cell': 'Plasma cell',
        'Lymphocyte': 'Lymphocyte',
        'Erythroblast&BasophilicNormoblast': 'Early Erythroid Precursors',
        'BasNorm&PolychromNorm': 'Intermediate Erythroid Precursors',
        'PolychromNorm&Normoblast': 'Late Erythroid Precursors',
        'Erythroblast': 'Erythroblast',
        'Basophilic normoblast': 'Basophilic normoblast',
        'Polychromatophilic normoblast': 'Polychromatophilic normoblast',
        'Normoblast': 'Normoblast',
        'Monocyte': 'Monocyte',
        'Eosinophil': 'Eosinophil',
        'Basophil': 'Basophil',
        'Mast cell': 'Mast cell',
        'Megakaryocytes': 'Megakaryocytes Count'
    }

    doc_title = {
        'Total Myeloid': 'Total Myeloid Lineage',
        'Total Erythroid': 'Total Erythroid Lineage',
        'Blast': 'Blast',
        'Promyelocyte': 'Promyelocyte',
        'MyeloMetamyelo': 'Intermediate Myeloid Precursors (Myelocyte and Metamyelocyte)',
        'Myelocyte': 'Myelocyte',
        'Metamyelocyte': 'Metamyelocyte',
        'Neutrophil': 'Total neutrophils (Band and Segmented neutrophils)',
        'Band neutrophil': 'Band neutrophil',
        'Segmented neutrophil': 'Segmented neutrophil',
        'Plasma cell': 'Plasma cell',
        'Lymphocyte': 'Lymphocyte',
        'Erythroblast&BasophilicNormoblast': 'Early Erythroid Precursors (Erythroblast and Basophilic normoblast)',
        'BasNorm&PolychromNorm': 'Intermediate Erythroid Precursors (Basophilic and Polychromatophilic normoblast)',
        'PolychromNorm&Normoblast': 'Late Erythroid Precursors (Polychromatophilic normoblast and Normoblast)',
        'Erythroblast': 'Erythroblast',
        'Basophilic normoblast': 'Basophilic normoblast',
        'Polychromatophilic normoblast': 'Polychromatophilic normoblast',
        'Normoblast': 'Normoblast',
        'Monocyte': 'Monocyte',
        'Eosinophil': 'Eosinophil',
        'Basophil': 'Basophil',
        'Mast cell': 'Mast cell',
        'Megakaryocytes': 'Megakaryocytes Count [MK / 10X FOV]'
    }

    # Define the path to your MK raw file
    # Adjust this path based on your script's execution directory
    mk_data_path = os.path.join(read_dir, 'mk_raw.csv')

    create_word_appendix(
        methd_comp=methd_comp,
        mk_csv_path=mk_data_path,
        output_filename=f'results/{save_name}_Appendix.docx',
        ordered_variables=appendix_order,
        fig_title_mapping=fig_titles,
        doc_title_mapping=doc_title
    )
