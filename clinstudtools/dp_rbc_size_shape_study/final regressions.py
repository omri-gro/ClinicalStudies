import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *


if __name__ == "__main__":
    variables = ['Spherocytes', 'Schistocytes', 'Macrocytes', 'Microcytes', 'Macroplatelets']

    meta_path = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\config.yaml'
    metadata = MetadataBundle(meta_path)


    df = read_to_df(f'DP_final_df.csv', file_dir=r'raw')
    methd_comp = MethodComparator(df)
    methd_comp.batch_fit('Investigators', 'CBM', variables)
    methd_comp.batch_compare(levels_a='Investigators', levels_b='CBM', variables=variables, split_by='Reviewing Site')
    methd_comp.batch_compare(levels_a='Investigators', levels_b='CBM', variables=variables, split_by='Site')
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.save_results(rf'results/dp_reg.csv')
    methd_comp.save_results(rf'results/dp_reg_bias.xlsx', result_type='bias')
    methd_comp.plot_all_regressions(f'results/dp_reg.pdf')


    # ---------------------------------------------------------
    # Generate MS Word Appendix
    # ---------------------------------------------------------
    sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\BMA_study2')
    from appendix_generator import create_word_appendix

    appendix_order = variables

    fig_titles = {
        'Spherocytes': 'Spherocytes',
        'Schistocytes': 'Schistocytes',
        'Macrocytes': 'Macrocytes',
        'Microcytes': 'Microcytes',
        'Macroplatelets': 'Macroplatelets',
    }

    doc_title = {
        'Spherocytes': 'Spherocytes',
        'Schistocytes': 'Schistocytes',
        'Macrocytes': 'Macrocytes',
        'Microcytes': 'Microcytes',
        'Macroplatelets': 'Large Platelets (Macroplatelets)',
    }

    create_word_appendix(
        methd_comp=methd_comp,
        mk_csv_path=None,
        output_filename=f'results/regressions_Appendix.docx',
        ordered_variables=appendix_order,
        fig_title_mapping=fig_titles,
        doc_title_mapping=doc_title,
        ref_arm_name='Digital Viewer Reviews [%]',
        test_arm_name='CBM Analyzer [%]',
    )




