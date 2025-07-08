import json

import pandas as pd
from pandas import DataFrame
import os

# parameters that are not manually graded by user, affect maturation DSS
lineage_vars = ['promyelocyte_grade', 'myelocyte_grade', 'metamyelocyte_grade', 'band_neutrophil_grade',
                'segmented_neutrophil_grade', 'erythroblast_grade', 'basophilic_normoblast_grade',
                'polychromatophilic_normoblast_grade', 'normoblast_grade']


def clinical_sign_off(rev, percent=True, with_tot=False):
    """
    Convert the standard BMA CASE_SIGNED_OFF review payload to dictionary (currently only Romanowsky)
    Args:
        rev (dict): content of event['payload']['payload']['review']
        percent (bool): for differential count classes, include percentages instead of absolute count
        with_tot (bool): whether to include lineage aggregates in results
    Returns:
        clin (dict): contains all case review clinical parameters
    """
    na_val = 'NA'  # what should be returned when grade not assigned
    diff_var = 'percentage' if percent else 'count'
    clin = {}
    obs = rev['observations']
    for ob in obs:
        param = ob['id']
        ob_type = ob['obs_type']
        if ob['group'] == 'prussian':  # remove if ever decide to cover Prussian in function
            continue
        na_str = f"\033[91mGrade not assigned for {param}.\033[0m"  # what should be printed when grade not assigned

        # parameters with only count (and total nucleated)
        if ob_type == 'additional' or param == 'total_differential':
            clin[param] = ob['count']

        # differential count parameters
        elif ob_type == 'differential':
            # exclude aggregates unless told otherwise
            if with_tot or not ob['isAggregation']:
                clin[param] = ob[diff_var]

        # grading parameters
        if ob_type in ['select', 'select_percent', 'select_density', 'select_ratio']:
            if param in lineage_vars:
                grade = ob['grade']
            elif 'selected' not in ob.keys():  # check for missing grade
                clin[param] = na_val
                print(na_str)
                continue
            else:
                grade = ob['selected']

            # convert numeric grades to string
            if isinstance(grade, int):
                grade_options = {item['id']: item['name'] for item in ob['options']['ranges']}
                clin[param] = grade_options[grade]
            else:
                clin[param] = grade

        # quantitative parameters grade is based on
        if ob_type == 'select_density':
            clin['megakaryocytes_density'] = ob['result']['count_per_10_hpf']
        elif ob_type == 'select_ratio':
            if ob['result']['denominator'] == 0:
                clin['me_ratio'] = 1000
            else:
                clin['me_ratio'] = ob['result']['nominator'] / ob['result']['denominator']
    # include reviewer's comments
    if 'comments' in rev.keys():
        comments = [com['message'] for com in rev['comments'] if isinstance(com, str)]
        clin['comments'] = '; '.join(comments)
    else:
        clin['comments'] = ''

    return clin


def parse_events(src_path):
    """
    Retrieve all needed data from a case's events file
    Args:
        src_path (str): path containing json file
    Returns:
        case_dict (dict):
    """
    case_dict = []

    evjson = json.load(open(src_path))
    for event in evjson[::-1]:
        if event['event_type'] == 'CASE_SIGNED_OFF':
            pay = event['payload']['payload']
            # ignore sign-offs performed by Scopio users
            reviewer = pay['session_info']['reviewer']['given_name']
            if reviewer == 'Scopio':
                continue

            # get clinical results
            clin_dict = clinical_sign_off(pay['review'])

            # add review metadata (consider adding scan IDs, age, sex and report number (from attachments link))
            clin_dict['reviewer'] = reviewer
            clin_dict['case_id'] = pay['session_info']['case_data']['case_id']
            clin_dict['case_name'] = pay['session_info']['case_data']['reference']
            clin_dict['sign_time'] = event['registered_at'][:16]  # only by minutes, not seconds

    case_dict = clin_dict
    return case_dict

if __name__ == "__main__":
    uuid = '407096ea-7124-4243-974d-b69f98911748'

    # Define the project root dynamically
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    events_json_path = os.path.join(project_root, rf"data/sandbox/{uuid}.json")
    needed_info = parse_events(events_json_path)
    df = pd.DataFrame([needed_info])
    df.to_csv('example.json', index=False)
    print(needed_info)
os.getcwd()
