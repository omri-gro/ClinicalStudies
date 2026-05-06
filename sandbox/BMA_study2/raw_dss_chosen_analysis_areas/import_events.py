import requests
import os
from pathlib import Path
import pandas as pd


# Define the target directory
# Note: Ensure you have write permissions to the root /events directory
output_dir = Path(r"C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\sandbox\BMA_study2\raw_dss_chosen_analysis_areas\sb1060_events")


def download_events(uuids, session):
    # Create directory if it doesn't exist
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        print(f"Error: Permission denied. Cannot create {output_dir}")
        return

    for uuid in uuids:
        url = f"https://127.182.253.74/analysis/scans/{uuid}/_events"
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

    df = pd.read_csv('sb1060_case_mapping.csv')
    uuid_list = df['scan'].tolist()

    session = 'userUUID=39760423-6ec0-48be-b9ca-7fe2102784e5; sid=s%3AwaBsSvGD2o9rhsjrIk2HcRuYPEsFlVM-.odYXUfiH%2Bw8yFbU8rEtMYRn6rigiz1%2Bj9HfLalpmkx0; swtoken=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IlNjb3Bpb1JEQHVzZXJzLm5vcmVwbHkuc2NvcGlvbGFicy5jb20iLCJjbGFpbXMiOnsicm9sZXMiOnsic2Nhbm5lciI6ImVkaXRvciIsImFuYWx5c2lzIjoiZWRpdG9yIn0sImNsaW5pYyI6IjVlMTEwMWY0LThjYzUtNDQ5NC05ZGE2LWVjNDMwMzFmY2QxMyJ9LCJpYXQiOjE3NzY5NTA2NDgsImV4cCI6MTc3Njk1MjQ0OH0.cixkM3xT9p9tNYtM75DHqOyxi89gno0A6Vs9fqQ4UL4'

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    download_events(uuid_list, session)
