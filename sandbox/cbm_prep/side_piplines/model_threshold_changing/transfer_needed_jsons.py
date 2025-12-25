import pandas as pd
import os
from shutil import copy2

import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from sandbox import read_to_df


if __name__ == "__main__":
    # study_uuids_dir = r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep\raw\cbm_method_comparison'
    study_uuids_dir = r'C:\Users\omrig\Downloads'
    study_uuids_name = r'scans_in_study.csv'
    uuids_df = read_to_df(study_uuids_name, file_dir=study_uuids_dir)
    uuids = uuids_df['ScanID']

    src_dir = r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\all6'
    dest_dir = r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\in_study'
    for uuid in uuids:
        src_path = os.path.join(src_dir, f'pbs\{uuid}.json')
        dst_path = os.path.join(dest_dir, f'pbs\{uuid}.json')
        if os.path.exists(src_path):
            copy2(src_path, dst_path)
        else:
            print(f'{uuid}.json not found')

        src_path = os.path.join(src_dir, f'logs\{uuid}.log')
        dst_path = os.path.join(dest_dir, f'logs\{uuid}.log')
        if os.path.exists(src_path):
            copy2(src_path, dst_path)
        else:
            print(f'{uuid}.log not found')


