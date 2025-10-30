import csv
import os
from datetime import timedelta
from time import time
from jose import jwt
import requests
import urllib3

urllib3.disable_warnings()

AUTH_SECRET = os.environ.get('AUTH_SECRET_TOKEN')


def get_auth_headers(swtoken):
    # current_time = time()
    # auth_token = jwt.encode(
    #     {
    #         'email': 'test@example.com',
    #         'claims': {'roles': {'analysis': 'editor'}},
    #         'iat': current_time,
    #         'exp': current_time + timedelta(minutes=3000).seconds,
    #     },
    #     AUTH_SECRET,
    # )

    headers = {'accept': 'application/json',
               'authorization': f'Bearer {swtoken}'}
    return headers


def get_obs_values(obs, obses):
    obs_res = next((o for o in obses if o['observation_id'] == obs['observation_id']), None)
    obs_type = obs['obs_type']
    # print(f'Processing observation {obs["observation_id"]} of type {obs_type}')
    values = []

    if obs_type in {'differential', 'additional'}:
        assert obs['result']['count'] == obs_res['result']['count']

        values += [obs_res['result']['count'], obs_res['result']['percentage'] / 100]

    elif obs_type in {'density'}:
        count = obs.get('result', {}).get('count')
        if count is None and obs_res:
            count = obs_res.get('result', {}).get('count')
        values += [count]

    elif obs_type in {'tri-scale', 'indicator'}:
        values += [obs.get('selected_value')]

        if obs_res:
            result = obs_res.get('result', {})
            percentage = result.get('percentage')
            count = result.get('count')

            if percentage is not None:
                if obs_res.get('link', {}).get('link_type') == 'feature_to_classes':
                    values += [percentage / 100]  # Normalize to 0-1
                else:
                    values += [percentage]  # RBC - Already in 0-1 range
            elif count is not None:
                values += [count]
    else:
        raise ValueError(f'Unsupported obs type: {obs_type}')

    if obs_res and obs_res.get('link', {}).get('link_type') == 'feature_to_classes':
        result = obs_res.get('result', {})
        nom = result.get('nominator') or result.get('numerator') or 0
        denom = result.get('denominator') or 0
        values += [nom, denom]

    return values


def get_scan_row(case, signoff_obses, obses, analysis_clump, curr_clumps, username, signoff_time):
    res = []
    res.append(case['case_id'])
    res.append(case['reference_id'])
    res.append(case['scans'][0]['ready_at'])
    res.append(signoff_time)
    res.append(case['scans'][0]['scan_id'])
    res.append(username)
    res.append(case['flags'])
    for obs in signoff_obses:
        _ = get_obs_values(obs, obses)
        res.extend(_)

    # add platelet satellitism
    platelet_satellitism = next((obs for obs in obses if obs['observation_id'] == 'platelet_satellitism'), None)
    if platelet_satellitism:
        res.append(platelet_satellitism.get('result', {}).get('count', 0))

    # add clumps count
    res.append(len(analysis_clump))
    res.append(len(curr_clumps))

    print(f'Processed scan {case["scans"][0]["scan_id"]}')
    return res


def get_obs_headers(signoff_obses, obses):
    headers = []
    obses_ids = {obs["observation_id"] for obs in obses}
    signoff_obses_ids = {obs["observation_id"] for obs in signoff_obses}
    missing_in_signoff = obses_ids - signoff_obses_ids
    print(f'Found {len(missing_in_signoff)} observations in obses that are not in signoff_obses')

    for obs in signoff_obses:
        obs_res = next((o for o in obses if o['observation_id'] == obs['observation_id']), None)
        obs_id = obs["observation_id"]
        obs_type = obs["obs_type"]

        if obs_type in {'differential', 'additional'}:
            headers += [f'{obs_id}_count', f'{obs_id}_percentage']

        elif obs_type in {'density'}:
            headers += [f'{obs_id}_count']

        elif obs_type in {'tri-scale', 'indicator'}:
            headers += [f'{obs_id}_grade']

            if obs_res:
                if 'result' in obs_res and 'percentage' in obs_res['result']:
                    headers += [f'{obs_id}_percentage']
                elif 'count' in obs_res['result']:
                    headers += [f'{obs_id}_count']

        else:
            raise ValueError(f'Unhandled obs_type {obs_type}')

        if obs_res and obs_res.get('link', {}).get('link_type') == 'feature_to_classes':
            headers += [f'{obs_id}_nominator', f'{obs_id}_denominator']

    headers += ['platelet_satellitism_count', 'analysis_clumps_count', 'curr_clumps_count']
    return headers


# @click.command()
# @click.option('--hosts', required=True, help='Comma-separated list of hosts')
# @click.option('--output_dir', required=True, type=click.Path(), help='Path to the output files directory')
# @click.option('--lim', required=False, default=200, help='Maximum number of cases to export')
def main(hosts, output_dir, swtoken, lim):
    hosts = hosts.split(',')
    print(f'Exporting CBM DSS to CSV for hosts: {hosts}')
    for host in hosts:
        headers = get_auth_headers(swtoken)
        # in CBM trial all cases are auto signed off so this covers all relevant cases
        resp = requests.get(f'https://{host}/review/cases?filters=tools_include:CBM&limit={lim}', headers=headers,
                            verify=False)
        cases = resp.json()['items']
        print(f'Found {len(cases)} cases for host {host}')
        # filter out cases with more than one scan
        filtered_cases = [case for case in cases if len(case['scans']) == 1]
        if len(filtered_cases) != len(cases):
            print(f'Found {len(cases) - len(filtered_cases)} cases with more then one scan, omitting them')

        base_file_name = host.replace(".", "_").replace(":", "_")
        output_file = f'{output_dir}/{base_file_name}_cbm.csv'
        headers_written = False
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            for case in filtered_cases:
                scan_id = case['scans'][0]['scan_id']
                resp = requests.get(f'https://{host}/analysis/scans/{scan_id}/_events', headers=headers, verify=False)
                events = resp.json()
                resp = requests.get(f'https://{host}/observations/scans/{scan_id}?variant=CBM', headers=headers,
                                    verify=False)
                obses = resp.json()
                resp = requests.get(f'https://{host}/observations/legend/blood?flat=true', headers=headers,
                                    verify=False)
                legend = resp.json()
                resp = requests.get(f'https://{host}/clumps/scans/{scan_id}', headers=headers, verify=False)
                curr_clumps = resp.json()['clumps']
                flattened_obses = []
                for cell_obses in obses.values():
                    if not isinstance(cell_obses, dict):
                        print(f'Unexpected error for case {case["case_id"]}, obses: {obses}')
                        continue

                    flattened_obses.extend(cell_obses.get('observations', []))

                # enrich obses with legend by matching legend 'id' with 'observation_id' in flattened_obses
                for obs in flattened_obses:
                    # Find matching legend entry
                    matching_legend = next((l for l in legend if l.get('id') == obs['observation_id']), {})
                    obs['link'] = matching_legend.get('link', {})
                    obs['obs_type'] = matching_legend.get('obs_type', 'unknown')

                analysis_clumps = []
                for event in events:
                    if event['event_type'] == 'CLUMP_ANALYSIS_COMPLETED':
                        analysis_clumps = event['payload']['clumps']

                for event in events:
                    if not (event['event_type'] == 'SESSION_PROGRESS_UPDATED' and event['payload']['step'] == 'SIGN_OFF'):
                        continue

                    signoff_obses = sorted(event['payload']['payload']['review']['observations'],
                                           key=lambda x: (x['group'], x['obs_type'], x['observation_id']))

                    if not headers_written:
                        csv_headers = ['case_id', 'barcode', 'scan time', 'sign off time', 'scan_id', 'user',
                                       'flags'] + get_obs_headers(signoff_obses, flattened_obses)
                        writer.writerow(csv_headers)
                        headers_written = True

                    username = event['payload']['payload']['session_info']['reviewer']['username']
                    try:
                        scan_row = get_scan_row(case, signoff_obses, flattened_obses, analysis_clumps, curr_clumps,
                                                username, signoff_time=event['registered_at'])
                        writer.writerow(scan_row)
                    except Exception as e:
                        print(f'Error processing scan {scan_id}: {e}')

        print(f'Exported CBM data to {output_file}')


if __name__ == '__main__':
    swtoken = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IlNjb3Bpb1JEQHVzZXJzLm5vcmVwbHkuc2NvcGlvbGFicy5jb20iLCJjbGFpbXMiOnsicm9sZXMiOnsic2Nhbm5lciI6ImVkaXRvciIsImFuYWx5c2lzIjoiZWRpdG9yIn0sImNsaW5pYyI6IjA4MDk5NTk2LTNhMmMtNDk0Ni04ZGY3LTJmMGFlYjNkOWY3ZSJ9LCJpYXQiOjE3NjE3MzU4MzQsImV4cCI6MTc2MTczNzYzNH0.Wx2PFsz4eulr8vbSBqiCUqqc77ycmn4GLw6EvhpyjL0"

    dest_dir = os.path.abspath(os.path.dirname(__file__))
    main('maintenance.scopiolabs.com:27599', dest_dir, swtoken, 500)


