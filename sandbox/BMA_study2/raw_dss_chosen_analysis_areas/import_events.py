import requests
import os
from pathlib import Path
import pandas as pd


# Define the target directory
# Note: Ensure you have write permissions to the root /events directory
output_dir = Path(r"C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\BMA_study2\raw_dss_chosen_analysis_areas\scopiobox1271_events")


def download_events(uuids, session):
    # Create directory if it doesn't exist
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        print(f"Error: Permission denied. Cannot create {output_dir}")
        return

    for uuid in uuids:
        # url = f"https://127.83.194.14/analysis/scans/{uuid}/_events"
        url = f"https://maintenance.scopiolabs.com:36198/analysis/scans/{uuid}/_events"
        output_path = output_dir / f"{uuid}.json"

        print(f"Fetching {uuid}...")

        try:
            # Create a headers dictionary with your copied cookie
            headers = {
                # Replace the string below with your actual cookie string
                "Cookie": f"session_id={session}"
            }

            response = requests.get(
                url,
                headers=headers,  # Pass the browser's credentials
                verify=False,
                timeout=10
            )

            # Raise an error for bad status codes (4xx, 5xx)
            response.raise_for_status()

            # Save the content
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(response.text)
            print(f"Successfully saved to {output_path}")

        except requests.exceptions.RequestException as e:
            print(f"Failed to download {uuid}: {e}")


if __name__ == "__main__":
    # Suppress InsecureRequestWarning if using verify=False
    import urllib3

    df = pd.read_csv('sb1271_case_mapping.csv')
    uuid_list = df['scan'].tolist()

    session = 'userUUID=8225eecf-01db-46d9-b89e-3dfbee07ddc5; sid=s%3Au0FPkJ1tD7F4L6AXIbL2GI4nQLK_fEtA.qB6VIaxIw1H62fIdZHinwPgruirXiu01L8hIsBHobxI; swtoken=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IlNjb3Bpb1JEQHVzZXJzLm5vcmVwbHkuc2NvcGlvbGFicy5jb20iLCJjbGFpbXMiOnsicm9sZXMiOnsic2Nhbm5lciI6ImVkaXRvciIsImFuYWx5c2lzIjoiZWRpdG9yIn0sImNsaW5pYyI6IjE0MzhkOTYxLWI2MjEtNGE1YS05ODZkLWU2MTBhMTI1NzdhOSJ9LCJpYXQiOjE3NzkwMjQzODIsImV4cCI6MTc3OTAyNjE4Mn0.mPGC_p5Ollin29Dc0qASH8ALgc-h4w0UzfNlM_gpbSE'

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    download_events(uuid_list, session)
