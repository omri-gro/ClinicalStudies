from gcs_read import *
import json
import pandas as pd
from google.cloud import storage


def get_latest_labels_blob(bucket: storage.Bucket, base_prefix: str, uuid: str) -> storage.Blob:
    """
    Searches GCS for the latest labels.json file for a specific UUID.
    """
    if not base_prefix.endswith('/'):
        base_prefix += '/'

    # Target exactly this UUID's directory to minimize API scanning time
    uuid_prefix = f"{base_prefix}{uuid}/"
    blobs = list(bucket.list_blobs(prefix=uuid_prefix))

    # Filter for labels files
    label_blobs = [b for b in blobs if b.name.endswith('labels.json')]

    if not label_blobs:
        return None

    # Since datetime folders sort alphabetically, the last one is the most recent
    latest_blob = sorted(label_blobs, key=lambda b: b.name)[-1]
    return latest_blob


def count_megas_in_roi(
        csv_path: str,
        bucket_name: str,
        gcs_base_prefix: str,
        output_csv_path: str
) -> pd.DataFrame:
    """
    Reads a CSV of ROIs, fetches labels.json from GCS, counts megakaryocytes
    within the ROI, and outputs an updated CSV.
    """
    # 1. Load the CSV
    df = pd.read_csv(csv_path)

    # Initialize GCS client
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    mega_counts = []

    print(f"Processing {len(df)} ROIs from CSV...")

    # 2. Iterate over each row
    for index, row in df.iterrows():
        uuid = str(row['scanUUID']).strip()
        xleft = row['xleft']
        ytop = row['ytop']
        xright = row['xright']
        ybottom = row['ybottom']

        # 3. Fetch the latest labels.json for this UUID
        blob = get_latest_labels_blob(bucket, gcs_base_prefix, uuid)

        if not blob:
            print(f"[{uuid}] WARNING: No labels.json found in GCS. Recording as NaN.")
            mega_counts.append(pd.NA)
            continue

        # Download and parse JSON directly into memory
        labels_text = blob.download_as_text()
        labels_data = json.loads(labels_text)

        # 4. Count megakaryocytes in the ROI boundaries
        count = 0

        # We can look directly at the "megakaryocyte" top-level key in labels.json
        megas = labels_data.get("megakaryocyte", [])

        for cell in megas:
            cx = cell.get("center_x")
            cy = cell.get("center_y")

            if cx is not None and cy is not None:
                # Check if the center of the cell falls within the CSV bounding box
                if xleft <= cx <= xright and ytop <= cy <= ybottom:
                    count += 1

        print(f"[{uuid}] Found {count} megakaryocytes in ROI.")
        mega_counts.append(count)

    # 5. Append the counts as a new column
    df['ROI_Megakaryocyte_Count'] = mega_counts

    # 6. Save the updated DataFrame
    df.to_csv(output_csv_path, index=False)
    print(f"\nDone! Updated data saved to '{output_csv_path}'")

    return df


if __name__ == "__main__":
    INPUT_CSV = "new_rois.csv"
    OUTPUT_CSV = "new_rois_with_counts.csv"
    BUCKET = "scopio_event_sync_hematology_prod"
    BASE_DIR = "events/scopiobox3217"

    updated_df = count_megas_in_roi(INPUT_CSV, BUCKET, BASE_DIR, OUTPUT_CSV)
    print(updated_df.head())
