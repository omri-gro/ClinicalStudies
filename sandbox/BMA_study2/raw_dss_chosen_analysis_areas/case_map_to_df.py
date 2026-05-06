import json
import pandas as pd


def process_case_mapping(file_path):
    # Load the JSON data
    with open(file_path, 'r') as file:
        data = json.load(file)

    rows = []

    # Iterate through each case in the JSON
    for item in data:
        # Filter only for 'romanowsky' review_flow
        if item.get("review_flow") == "romanowsky":
            case_id = item.get("case_id")
            status = item.get("status")
            created_at_raw = item.get("created_at")

            # Handle null values for barcode and reference by replacing with empty strings
            barcode = item.get("barcode") if item.get("barcode") is not None else ""
            reference = item.get("reference") if item.get("reference") is not None else ""

            # Create a new row for each scan in the scans array
            for scan in item.get("scans", []):
                rows.append({
                    "case_id": case_id,
                    "scan": scan,
                    "status": status,
                    "created_at": created_at_raw,
                    "barcode": barcode,
                    "reference": reference
                })

    # Create the DataFrame
    df = pd.DataFrame(rows)

    # Format the 'created_at' column to dd.mm.yyyy HH:MM
    if not df.empty:
        df['created_at'] = pd.to_datetime(df['created_at']).dt.strftime('%d.%m.%Y %H:%M')

    return df


# Example usage:
if __name__ == "__main__":
    df = process_case_mapping('sb1060_case_mapping.json')
    print(df.head(10))
    df.to_csv(f"sb1060_case_mapping.csv", index=False)
