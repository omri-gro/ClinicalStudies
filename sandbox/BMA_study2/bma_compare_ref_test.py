import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle


# def bma_df_pipeline(paths: dict, metadata: sb.MetadataBundle, dir=None, measurement_col='Value',
#                     more_id_vars=None):



if __name__ == "__main__":
    meta_path = r'config_BMA.yaml'
    site = 'OHSU'
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    metadata = MetadataBundle(meta_path)
    srcs = {(site, 'REF'): f'{site}_CRF_REF.csv',
            (site, 'TEST'): f'{site}_CRF_TEST.csv'}
    methd_comp = MethodComparator.from_paths_dict(srcs, metadata, dir=cur_dir, bma=True)
    ndc_vars_list = metadata.variable_groups['NDC']
    needed_rows = methd_comp.df[methd_comp.df['Variable'].isin(ndc_vars_list)]

    # calculate means for variable-identifier combinations that have exactly 2 rows
    variable_identifier = ['SampleID', 'Site', 'Method', 'Variable']
    m = needed_rows.groupby(variable_identifier)['Value'].transform('size').eq(2)
    out = needed_rows.loc[m].groupby(variable_identifier, as_index=False)['Value'].mean()
    methd_comp.df = out

    methd_comp.batch_fit(['REF'], ['TEST'], ndc_vars_list)
    methd_comp.save_results(rf'results/{site}_bma_reg.csv')
    methd_comp.plot_all_regressions(f'results/{site}_bma_reg.pdf')

    print(methd_comp)
    print()


