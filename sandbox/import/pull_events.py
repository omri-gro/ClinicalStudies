import requests
import os



def request_text(url, headers):
    try:
        response = requests.get(url, headers=headers, verify=False)  # change verify if url not trusted
        if response.status_code == 200:
            return response.text
        else:
            print(f"Request failed for {url}: {response.status_code}")
            print(response.text)
            return None
    except requests.exceptions.RequestException as e:
        print(f"\033[93mError for {url}: {e}\033[0m")
        return None

def event_from_api(domain, uuid, token):
    headers = {"Accept": "application/json", "Authorization": f"Bearer {token}"}
    events_url = fr'https://{domain}/analysis/scans/{uuid}/_events'
    result = request_text(events_url, headers)
    if result is None:
        print(f"\033[91mNo events returned for UUID {uuid}\033[0m")
    return result

def obs_from_api(domain, uuid, token, variant="PBS"):
    headers = {"Accept": "application/json", "Authorization": f"Bearer {token}"}
    events_url = fr'https://{domain}/observations/scans/{uuid}?variant={variant}'
    result = request_text(events_url, headers)
    if result is None:
        print(f"\033[91mNo obs returned for UUID {uuid}\033[0m")
    return result


def filenames_in_dir(dir_path, needed_ext=None):
    # extension should include the '.'
    names_list = []
    for dir_item in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, dir_item)):
            # Split the filename into root and extension
            root, extension = os.path.splitext(dir_item)
            if needed_ext is None or extension == needed_ext:
                names_list.append(root)
    if not names_list:
        print(f"\033[91mNo appropriate files found in {dir_path}\033[0m")
    return names_list


if __name__ == '__main__':
    port = 31303
    swtoken = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IlNjb3Bpb1JEQHVzZXJzLm5vcmVwbHkuc2NvcGlvbGFicy5jb20iLCJjbGFpbXMiOnsicm9sZXMiOnsic2Nhbm5lciI6ImVkaXRvciIsImFuYWx5c2lzIjoiZWRpdG9yIn0sImNsaW5pYyI6IjgxZjJkYTk0LTE0OTMtNDY3Zi1hMzYwLWZkZWYzZjM0ZjBmMCJ9LCJpYXQiOjE3NjQyNDk4NTgsImV4cCI6MTc2NDI1MTY1OH0.cKYv78n8cQI_Akd05F1ZXIB4eAx38seibX8WXYkCsiU"
    batch_name = 'Hamburg_validation'

    machine_domain = f"maintenance.scopiolabs.com:{port}"
    uuids_dir = fr'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\{batch_name}\pbs'
    save_dir = fr'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\{batch_name}\pbs_obs'

    # list of requested scan IDs
    for uuid in filenames_in_dir(uuids_dir, '.json'):
        events_f = obs_from_api(machine_domain, uuid, swtoken)
        if events_f is not None:
            # Save to file with the ID as filename
            out_path = os.path.join(save_dir, f"{uuid}.json")
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(events_f)









