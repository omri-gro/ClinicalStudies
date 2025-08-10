import pandas as pd
import os
import json
from collections import Counter

def read_json(path):
    data = json.load(open(path))
    return data


def diff_from_cells_json(data, result_type='percentage', dist_id='all_selected'):
    # result_type is 'count' or 'percentage'
    for dist in data:
        if dist['distribution_id'] == dist_id:
            dist_results_dict = {}
            for part in dist['observations']:
                if 'result' in part.keys():
                    if result_type in part['result'].keys():
                        dist_results_dict[part['observation_id']] = part['result'][result_type]
                else:
                    dist_results_dict[part['observation_id']] = part['grade']
    return dist_results_dict


def first_last_counts(data):
    # return 2 dicts by inserting data from a labels json
    orig_obs = [lbl['history'][0]['class_id'] for lbl in data['labels']]
    orig_obs_dict = dict(Counter(orig_obs))
    fnl_obs = [lbl['history'][-1]['class_id'] for lbl in data['labels']]
    fnl_obs_dict =  dict(Counter(fnl_obs))
    return orig_obs_dict, fnl_obs_dict


def save_dict_as_series(data_dict, save_path, ext='xlsx'):
    df = pd.Series(data_dict)
    if ext == 'csv':
        df.to_csv(f'{save_path}.{ext}')
    elif ext in ['xls', 'xlsx', 'xlsm']:
        df.to_excel(f'{save_path}.{ext}')
    else:
        print(f'extension {ext} unknown')




multi_uuid = '9471e325-ada8-4eab-b82d-59cfdbf81707'
scan_uuid = '47ade2d0-94e5-495f-8e56-076fb599f7e9'

multi_events = read_json(f'{multi_uuid}_events.json')
scan_events = read_json(f'{scan_uuid}_events.json')
scan_cells = read_json(f'{scan_uuid}_cells.json')
multi_labels = read_json(f'{multi_uuid}_labels.json')

final_diff = diff_from_cells_json(scan_cells)
final_counts = diff_from_cells_json(scan_cells, result_type='count')
not_only_selected_counts = diff_from_cells_json(scan_cells, result_type='count', dist_id='all')

orig_obs, fnl_obs = first_last_counts(multi_labels)

save_dict_as_series(final_diff, 'final_diff')
save_dict_as_series(final_counts, 'final_counts')
save_dict_as_series(not_only_selected_counts, 'not_only_selected_counts')
save_dict_as_series(orig_obs, 'orig_obs')
save_dict_as_series(fnl_obs, 'fnl_obs')

print('a')

