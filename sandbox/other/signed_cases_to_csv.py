import sys
import os
import json
from pandas import DataFrame
from datetime import datetime, timedelta
sys.path.append(os.path.abspath(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\import'))
from pull_events import filenames_in_dir
from export_cbm_dss_to_csv import get_obs_values, get_scan_row, get_obs_headers


# def scan_info_collect(uuid, obs, events, legend, headers=None):
#     sign_offs = []
#     flattened_obses = []
#     for cell_obses in obs.values():
#         if not isinstance(cell_obses, dict):
#             print(f'Unexpected error for case {uuid}, obses: {obses}')
#             continue
#         flattened_obses.extend(cell_obses.get('observations', []))
#
#     # enrich obses with legend by matching legend 'id' with 'observation_id' in flattened_obses
#     for obs in flattened_obses:
#         # Find matching legend entry
#         matching_legend = next((l for l in legend if l.get('id') == obs['observation_id']), {})
#         obs['link'] = matching_legend.get('link', {})
#         obs['obs_type'] = matching_legend.get('obs_type', 'unknown')
#
#     analysis_clumps = []
#     for event in events:
#         if event['event_type'] == 'CLUMP_ANALYSIS_COMPLETED':
#             analysis_clumps = event['payload']['clumps']
#
#     for event in events:
#         if not (event['event_type'] == 'SESSION_PROGRESS_UPDATED' and event['payload']['step'] == 'SIGN_OFF'):
#             continue
#
#         signoff_obses = sorted(event['payload']['payload']['review']['observations'],
#                                key=lambda x: (x['group'], x['obs_type'], x['observation_id']))
#
#
#
#         if headers is not None:
#             headers = ['case_id', 'barcode', 'scan time', 'sign off time', 'scan_id', 'user',
#                        'flags'] + get_obs_headers(signoff_obses, flattened_obses)
#
#         username = event['payload']['payload']['session_info']['reviewer']['username']
#         try:
#             res = []
#             res.append(uuid)
#             res.append()
#             scan_row = get_scan_row(case, signoff_obses, flattened_obses, analysis_clumps, curr_clumps,
#                                     username, signoff_time=event['registered_at'])
#             writer.writerow(scan_row)
#         except Exception as e:
#             print(f'Error processing scan {scan_id}: {e}')


def obs_pair(obs):
    name = obs.get('adjusted_name')
    if name is None:
        name = obs.get('name')

    if obs['obs_type'] in ['tri-scale', 'indicator']:
        if obs['selected_value'] != 'unset':
            val = obs['selected_value']
        else:
            val = None
    elif obs['obs_type'] == 'differential' and obs['observation_id'] != 'total_differential':
        val = obs['result']['percentage']
    else:
        val = obs['result']['count']
    return name, val


def adjst_timestamp(instr, infrmt="%Y-%m-%dT%H:%M:%S.%f", outfrmt="%Y-%m-%d %H:%M:%S", hour_move=0):
    # hour_move > 0 to add time, <0 to move back
    dt = datetime.strptime(instr, infrmt)
    dt = dt + timedelta(hours=hour_move)
    return dt.strftime(outfrmt)


def signed_off_extract(event,
                       # hours reporting arguments
                       outfrmt="%Y-%m-%d %H:%M:%S", hour_move=0):
    pload = event['payload']['payload']
    observations = pload['review']['observations']
    ses_inf = pload['session_info']
    rev = ses_inf['reviewer']
    case_data = ses_inf['case_data']

    rprt = {'Case Name': case_data['reference']}

    # gather observations
    for obs in observations:
        k, v = obs_pair(obs)
        rprt[k] = v

    # gather comments
    for comnt in pload['review']['comments']:
        rprt[comnt['name']] = comnt['message']

    # additional information
    add_inf = {
        'Scan Mode': case_data['scans'][0]['tool'],
        'Scan Time': adjst_timestamp(case_data['scans'][0]['ready_at'], outfrmt=outfrmt, hour_move=hour_move),
        'Sign-Off Time': adjst_timestamp(ses_inf['signed_at'], outfrmt=outfrmt, hour_move=hour_move, infrmt="%Y-%m-%dT%H:%M:%SZ"),
        'Signed-Off By': f"{rev['given_name']} {rev['family_name']}",
        'Scan ID': case_data['scans'][0]['scan_id'],
        'Case ID': case_data['case_id']}
    rprt.update(add_inf)
    return rprt



if __name__ == "__main__":
    session_name = 'Hamburg_validation'
    main_dir = rf'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\{session_name}'
    save_path_pbs = fr'{main_dir}\all_signoffs_pbs.csv'
    save_path_bma = fr'{main_dir}\all_signoffs_bma.csv'

    bma_uuids_dir = fr'{main_dir}\bma'
    pbs_uuids_dir = fr'{main_dir}\pbs'
    events_dir = fr'{main_dir}\events'
    obs_dir = fr'{main_dir}\pbs_obs'

    limit = 500
    hours_diff = -1

    # not in use
    pbs_legend_path = fr'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\import\pbs_legend.json'
    with open(pbs_legend_path) as json_file:
        legend = json.load(json_file)

    rows = []
    for uuid in filenames_in_dir(pbs_uuids_dir):
        print(uuid)
        event_path = os.path.join(events_dir, f'{uuid}.json')
        with open(event_path) as json_file:
            events = json.load(json_file)
            for event in events:
                if event['event_type'] == 'SESSION_PROGRESS_UPDATED' or event['payload'].get('step') == 'SIGN_OFF':
                    rows.append(signed_off_extract(event, hour_move=hours_diff))
    df = DataFrame(rows)
    df.to_csv(save_path_pbs, index=False)

    rows = []
    for uuid in filenames_in_dir(bma_uuids_dir):
        print(uuid)
        event_path = os.path.join(events_dir, f'{uuid}.json')
        with open(event_path) as json_file:
            events = json.load(json_file)
            for event in events:
                if event['event_type'] == 'SESSION_PROGRESS_UPDATED' or event['payload'].get('step') == 'SIGN_OFF':
                    rows.append(signed_off_extract(event, hour_move=hours_diff))
    df = DataFrame(rows)
    df.to_csv(save_path_bma, index=False)


