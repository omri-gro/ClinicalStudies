import pandas as pd
import os
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\cbm_prep')
from objects import MethodComparator
from sandbox import *
from pipelines import mean_manual_pipe, medium_pipe

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools import careful_map, apply_arbitration_override
from clinstudtools.transforms import filter_by_reference, filter_by_condition, filter_samples_by_condition
from clinstudtools.utils import read_to_df


if __name__ == "__main__":
    suffix = ''
    sites = ['BWH', 'CPG', 'HUP', 'LMU', 'SYN', 'TASMC']
    analysis_name = "cbm_method_comparison"
    meta_path = r'config.yaml'
    test_arm = 'CBM'
    cbm_version = 'v325'

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    os.chdir(os.path.join(cur_dir, ".."))
    raw_dir = os.path.join(cur_dir, r'raw', analysis_name)

    metadata = MetadataBundle(meta_path)

    investigators_map = {
        # Standardize typos/variations
        'Christopher Wright': 'Chris',
        'Christine Lavoie': 'Christine',
        'Ebikebuna Rufus': 'Ebi',
        'Ebikebuna Rufus F.': 'Ebi',
        'Ebikebuna Rufus F': 'Ebi',
        'Thu Tran': 'Thu',
        'THU TRAN': 'Thu',
        'thu tran': 'Thu',
        'Aubrey B Charlton': 'Aubrey',
        'Deborah Swearingen': 'Deborah',
        'Maria Buen Viana De Perio': 'Buen',
        'Joy Arthur': 'Joy',
        'Madison Brooks': 'Madison',
        'Michelle Huynh': 'Michelle',
        'Michelle huynh': 'Michelle',
        'Tiffany I Highsmith': 'Tiffany',
        'Tiffany I. Highsmith': 'Tiffany',
        'Alina KÃƒÂ¼pper': 'Alina',
        'Alina': 'Alina',
        'Sladana Nikolic': 'Sladana',
        'Nikolic Sladana': 'Sladana',
        'Sladana': 'Sladana',
        'Ana Catarina Silva': 'Ana',
        'Harsha Hirani': 'Harsha',
        'Harsha HIrani': 'Harsha',
        'Thomas Muddiman': 'Thomas',
        'Tony Omigie': 'Tony',
        'Sarah Pereira Rodrigues': 'Sarah',
        'YAEL SAYEGH': 'Yael',
        'Yael Sayegh': 'Yael',
        'Yael S': 'Yael',
        'YAEL ASYEGH': 'Yael',

        # Explicitly tag Arbitrators
        'Jared Block': 'Arbitrator',
        'Jennifer Egan': 'Arbitrator',
        'Dr. med. Weigand, Michael': 'Arbitrator',
        'Dr. med. Michael Weigand': 'Arbitrator',
        'Dr Guy Hannah': 'Arbitrator',
        'Andrew Doyle': 'Arbitrator',
        'Ben-Zion Katz': 'Arbitrator',
        'Dan BENISTY': 'Arbitrator',
        'Olga Pozdnyakova': 'Arbitrator',
        'Christopher Hergott': 'Arbitrator',
        'Robert P Hasserjian': 'Arbitrator',

        # Preserve system/automated roles
        'CBM': 'CBM',
        'Mean Investigator': 'Mean Investigator'}






