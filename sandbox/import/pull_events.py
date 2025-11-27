import requests
import os



def request_text(url, headers):
    try:
        response = requests.get(url, headers=headers, verify=False)
        if response.status_code == 200:
            return response.text
        else:
            print(f"Request failed for {url}: {response.status_code}")
            print(response.text)
    except requests.exceptions.RequestException as e:
        print(f"\033[93mError for {url}: {e}\033[0m")


if __name__ == '__main__':
    port_number = 20267
    swtoken = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IlNjb3Bpb1JEQHVzZXJzLm5vcmVwbHkuc2NvcGlvbGFicy5jb20iLCJjbGFpbXMiOnsicm9sZXMiOnsic2Nhbm5lciI6ImVkaXRvciIsImFuYWx5c2lzIjoiZWRpdG9yIn0sImNsaW5pYyI6IjgxZjJkYTk0LTE0OTMtNDY3Zi1hMzYwLWZkZWYzZjM0ZjBmMCJ9LCJpYXQiOjE3NjQyMjk2OTQsImV4cCI6MTc2NDIzMTQ5NH0.g9RDzWPs-ERQ3rI5Q5BVc3XOFo8nVm949FEZAEECBqU"

    # list of requested scan IDs




