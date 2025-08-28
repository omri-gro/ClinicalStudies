import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from objects import MethodComparator
from sandbox import MetadataBundle, read_to_df


if __name__ == "__main__":
    mthd = 'REF'
    site = 'OHSU'

    cell_classes = ['Blast', 'Promyelocyte', 'Myelocyte', 'Metamyelocyte', 'Segmented neutrophil', 'Band neutrophil', 'Plasma cell', 'Lymphocyte', 'Erythroblast', 'Basophilic normoblast', 'Polychromatophilic normoblast', 'Normoblast', 'Monocyte', 'Eosinophil', 'Basophil', 'Mast cell']

    data_file_name = f'df_pivoted_{mthd}.xlsx'
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    df = read_to_df(data_file_name, file_dir=cur_dir)
    df2 = pd.wide_to_long(df, "Value", i=["Sample", "Parameter"], j="Method", sep="_").reset_index()
    df2 = df2[df2['Parameter'].isin(cell_classes)]
    df2["Site"] = "OHSU"
    df2.rename(columns={"Sample": "SampleID", "Parameter": "Variable"}, inplace=True)
    df2['Method'] = df2['Method'].astype(str)
    df2['Value'] = df2['Value'].astype(float)
    mthd_comp = MethodComparator(df=df2)
    mthd_comp.batch_fit(["1"], ["2"], cell_classes)
    mthd_comp.save_results(rf'{site}_{mthd}_inter_rev_agreement.xlsx')
    mthd_comp.plot_all_regressions(rf'{site}_{mthd}_inter_rev_agreement.pdf')


