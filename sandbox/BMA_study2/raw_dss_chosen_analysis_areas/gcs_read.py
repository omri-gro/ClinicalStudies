import json
import yaml
import re
import glob
import os
from typing import Dict, List, Set, Tuple
from collections import Counter, defaultdict
import pandas as pd
from google.cloud import storage

# ==========================================
# 1. Updated Class Name Mapping Function
# ==========================================

# def load_class_mapping(filepath: str) -> Dict[int, str]:
#     """
#     Reads the romawosky_classes.json file and builds a mapping dictionary
#     from class_id to the human-readable display_name.
#     """
#     with open(filepath, 'r') as f:
#         data = json.load(f)
#
#     mapping = {}
#
#     # Iterate through all possible class lists in the JSON
#     for category in ['builtin_classes', 'custom_classes', 'placeholder_classes']:
#         for cls in data.get(category, []):
#             # Use adjusted_display_name if the user renamed it, otherwise fallback to default
#             name = cls.get('adjusted_display_name') or cls.get('display_name')
#             mapping[cls['id']] = name
#
#     return mapping

def load_class_mapping(filepath: str) -> Tuple[Dict[int, str], Set[str]]:
    """
    Reads the label_classes.yml file.
    Returns:
        - A mapping dictionary from class_id to display_name.
        - A set of display_names that have `in_differential: true`.
    """
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    mapping = {}
    differential_classes = set()

    for cls in data:
        class_id = cls.get('id')
        display_name = cls.get('display_name')

        if class_id is not None and display_name:
            mapping[class_id] = display_name
            if cls.get('in_differential') is True:
                differential_classes.add(display_name)

    return mapping, differential_classes



# ==========================================
# 2. Core Analysis Functions
# ==========================================
def extract_case_metadata(events: List[dict]) -> dict:
    """Recursively searches for key metadata nested deeply in the events list."""
    def find_first_key(obj, key):
        if isinstance(obj, dict):
            if key in obj:
                return obj[key]
            for k, v in obj.items():
                res = find_first_key(v, key)
                if res is not None:
                    return res
        elif isinstance(obj, list):
            for item in obj:
                res = find_first_key(item, key)
                if res is not None:
                    return res
        return None

    # Time of sign-off
    signed_at_raw = find_first_key(events, 'signed_at')
    signed_at = None
    if signed_at_raw:
        # Transforms ISO format '2026-01-30T23:19:23Z' to '2026-01-30 23:19:23'
        signed_at = signed_at_raw.replace('T', ' ').replace('Z', '').split('.')[0]

    # Barcode
    scan_info = find_first_key(events, 'scan_info')
    barcode = scan_info.get('slide_id') if isinstance(scan_info, dict) else None

    # Case ID
    case_data = find_first_key(events, 'case_data')
    case_id = case_data.get('case_id') if isinstance(case_data, dict) else None

    # Reviewer name
    reviewer_name = None
    for e in events:
        if e.get("event_type") == "SESSION_PROGRESS_UPDATED" and e.get("payload", {}).get("step") == "SIGN_OFF":
            reviewer = find_first_key(e, 'reviewer')
            if isinstance(reviewer, dict):
                given = reviewer.get('given_name', '')
                family = reviewer.get('family_name', '')
                reviewer_name = f"{given} {family}".strip()
            break

    return {
        "Sign-off Time": signed_at,
        "Barcode": barcode,
        "Case ID": case_id,
        "Reviewer": reviewer_name
    }


def get_particle_count(events: List[dict]) -> int:
    """
    Scans events from earliest to latest to find the feature_detection_results
    and returns the length of the particles list.
    """
    # Iterate backwards (from the earliest event at the bottom of the list, upwards)
    for event in reversed(events):
        try:
            # Attempt the deep dictionary lookup
            particles = event['payload']['analysis']['cv_params']['feature_detection_results']['analysis_data']['region_data']['particles']

            # If successful and it's a list, return its length immediately
            if isinstance(particles, list):
                return len(particles)

        except (KeyError, TypeError):
            # KeyError: A key in the chain is missing
            # TypeError: An intermediate value was None instead of a dict
            continue

    # Fallback if the particles list is entirely missing from the file
    return 0

def get_approved_rois_at_signoff(events: List[dict]) -> Tuple[Set[str], Dict[str, Tuple[int, int, int, int]], int]:
    sign_off_idx = next(
        (i for i, e in enumerate(events)
         if e.get("event_type") == "SESSION_PROGRESS_UPDATED"
         and e.get("payload", {}).get("step") == "SIGN_OFF"),
        None
    )

    if sign_off_idx is None:
        raise ValueError("No SIGN_OFF event found in the events log.")

    approved_roi_ids: Set[str] = set()
    for event in events[sign_off_idx:]:
        if event.get("event_type") == "SELECTED_REGIONS_UPDATED":
            selections = event.get("payload", {}).get("selection", [])
            approved_roi_ids = {
                roi["roi_id"] for roi in selections
                if roi.get("selection_status") == "SELECTED"
            }
            break

    roi_bounds: Dict[str, Tuple[int, int, int, int]] = {}
    investigator_roi_count = 0

    for event in events:
        if event.get("event_type") == "ANALYSIS_REGION_UPDATED":
            payload = event.get("payload", {})
            roi_id = payload.get("roi_id")
            bounds = payload.get("bounds")
            roi_type = payload.get("params", {}).get("roi_type")

            if roi_id in approved_roi_ids and roi_id not in roi_bounds and bounds and roi_type != "BLANK":
                x_min = bounds["xleft"]
                y_min = bounds["ytop"]
                x_max = x_min + bounds["width"]
                y_max = y_min + bounds["height"]
                roi_bounds[roi_id] = (x_min, y_min, x_max, y_max)

                # Check if investigator generated it
                if payload.get("auto_suggested") is not True:
                    investigator_roi_count += 1

    return approved_roi_ids, roi_bounds, investigator_roi_count

def get_rois_without_signoff(events: List[dict]) -> Tuple[Set[str], Dict[str, Tuple[int, int, int, int]], int]:
    """
    Extracts all ROIs in their latest configuration without relying on a SIGN_OFF event or SELECTED_REGIONS_UPDATED.
    """
    approved_roi_ids: Set[str] = set()
    roi_bounds: Dict[str, Tuple[int, int, int, int]] = {}
    investigator_roi_count = 0

    # Events are reverse-chronological. The first time we see an ROI ID, it is its final state.
    for event in events:
        if event.get("event_type") == "ANALYSIS_REGION_UPDATED":
            payload = event.get("payload", {})
            roi_id = payload.get("roi_id")
            bounds = payload.get("bounds")
            roi_type = payload.get("params", {}).get("roi_type")

            if roi_id and roi_id not in roi_bounds and bounds and roi_type != "BLANK":
                approved_roi_ids.add(roi_id)
                x_min = bounds["xleft"]
                y_min = bounds["ytop"]
                x_max = x_min + bounds["width"]
                y_max = y_min + bounds["height"]
                roi_bounds[roi_id] = (x_min, y_min, x_max, y_max)

                # Check if investigator generated it
                if payload.get("auto_suggested") is not True:
                    investigator_roi_count += 1

    return approved_roi_ids, roi_bounds, investigator_roi_count


def get_user_reclassifications_at_signoff(events: List[dict]) -> Dict[str, int]:
    sign_off_idx = next(
        (i for i, e in enumerate(events)
         if e.get("event_type") == "SESSION_PROGRESS_UPDATED"
         and e.get("payload", {}).get("step") == "SIGN_OFF"),
        0
    )

    reclassifications: Dict[str, int] = {}
    for event in events[sign_off_idx:]:
        if event.get("event_type") == "REGION_LABELS_USER_UPDATED":
            labels = event.get("payload", {}).get("labels", [])
            for label in labels:
                blob_id = label.get("blob_id")
                class_id = label.get("extra", {}).get("class_id")

                if blob_id and class_id is not None and blob_id not in reclassifications:
                    reclassifications[blob_id] = class_id

    return reclassifications


def get_user_reclassifications_without_signoff(events: List[dict]) -> Dict[str, int]:
    """
    Extracts the latest user reclassifications by scanning from the beginning (latest) of the file.
    """
    reclassifications: Dict[str, int] = {}

    for event in events:
        if event.get("event_type") == "REGION_LABELS_USER_UPDATED":
            labels = event.get("payload", {}).get("labels", [])
            for label in labels:
                blob_id = label.get("blob_id")
                class_id = label.get("extra", {}).get("class_id")

                if blob_id and class_id is not None and blob_id not in reclassifications:
                    reclassifications[blob_id] = class_id

    return reclassifications



def calculate_final_classifications(
        labels_data: dict,
        roi_bounds: Dict[str, Tuple[int, int, int, int]],
        reclassifications: Dict[str, int],
        with_reclass=False
) -> Dict[int, int]:
    class_counts = Counter()
    approved_boxes = list(roi_bounds.values())

    for cell_type, cells in labels_data.items():
        for cell in cells:
            cx = cell.get("center_x")
            cy = cell.get("center_y")
            blob_id = cell.get("blob_id")
            original_class_id = cell.get("extra", {}).get("class_id")

            if cx is None or cy is None or original_class_id is None:
                continue

            in_approved_area = any(
                x_min <= cx <= x_max and y_min <= cy <= y_max
                for x_min, y_min, x_max, y_max in approved_boxes
            )

            # megakaryocytes are counted regardless of roi
            if in_approved_area or original_class_id == 5003:
                if with_reclass:
                    final_class_id = reclassifications.get(blob_id, original_class_id)
                else:
                    final_class_id = original_class_id
                class_counts[final_class_id] += 1

    return dict(class_counts)


def process_single_scan(events_data: List[dict], labels_data: dict, with_reclass=False) -> Dict[int, int]:
    """Helper to run the analysis logic on loaded JSON objects."""
    try:
        approved_roi_ids, roi_bounds, inv_roi_count = get_approved_rois_at_signoff(events_data)
        reclassifications = get_user_reclassifications_at_signoff(events_data)
        counts = counts = calculate_final_classifications(labels_data, roi_bounds, reclassifications, with_reclass=with_reclass)
        metadata = extract_case_metadata(events_data)
        particles_count = get_particle_count(events_data)

        return {
            **metadata,
            "Approved ROIs": len(approved_roi_ids),
            "Investigator ROIs": inv_roi_count,
            "Particles": particles_count,
            **counts
        }
    except ValueError as e:
        print(f"Skipping a scan: {e}")
        return {}

def process_single_scan_no_signoff(events_data: List[dict], labels_data: dict, with_reclass=False) -> Dict[int, int]:
    """Helper to run the analysis logic on a single scan lacking a SIGN_OFF event."""
    try:
        approved_roi_ids, roi_bounds, inv_roi_count = get_rois_without_signoff(events_data)
        reclassifications = get_user_reclassifications_without_signoff(events_data)
        counts = calculate_final_classifications(labels_data, roi_bounds, reclassifications, with_reclass=with_reclass)
        metadata = extract_case_metadata(events_data)
        particles_count = get_particle_count(events_data)

        return {
            **metadata,
            "Approved ROIs": len(approved_roi_ids),
            "Investigator ROIs": inv_roi_count,
            "Particles": particles_count,
            **counts
        }
    except Exception as e:
        print(f"Skipping a scan due to error: {e}")
        return {}

# ==========================================
# 3. Fast GCS Traversal and DataFrame Generation
# ==========================================
def get_latest_signoff_blobs(bucket_name: str, prefix: str, folder_name: str = "SIGN_OFF") -> Dict[str, Dict[str, storage.Blob]]:
    """
    Scans the bucket and groups blobs. For each UUID that has a SIGNOFF folder,
    it identifies the latest datetime and returns references to its JSON blobs.
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # Ensure our base path ends with a slash for proper prefix matching
    if not prefix.endswith('/'):
        prefix += '/'

    # We will use a regex to parse the blob name.
    # Expected structure: some_path/<UUID>/<folder_name>/<datetime>/<filename>
    # Group 1: UUID, Group 2: Datetime, Group 3: Filename
    pattern = re.compile(rf"^{re.escape(prefix)}/?([^/]+)/{folder_name}/([^/]+)/(events\.json|labels\.json)$")

    # 1. Use delimiter='/' to get only the immediate "folders" under base_prefix
    iterator = bucket.list_blobs(prefix=prefix, delimiter='/')

    # We must iterate through the page to populate the `iterator.prefixes` property
    list(iterator)
    uuid_prefixes = iterator.prefixes

    print(f"Found {len(uuid_prefixes)} UUID directories starting with {prefix}. Skipping irrelevant folders and jumping to SIGN_OFFs...")

    latest_blobs_per_uuid = {}

    for uuid_prefix in uuid_prefixes:
        # Extract just the UUID string from the path
        uuid = uuid_prefix[len(prefix):].strip('/')

        # The exact, targeted path to the SIGN_OFF folder
        signoff_prefix = f"{uuid_prefix}{folder_name}/"

        # 2. List ONLY the files inside this specific UUID's <folder_name> folder
        # We don't need a delimiter here because we only care about what's inside SIGN_OFF
        signoff_blobs = list(bucket.list_blobs(prefix=signoff_prefix))

        # If the list is empty, this UUID didn't have a <folder_name> folder. We skip it instantly!
        if not signoff_blobs:
            continue

        datetimes_dict = defaultdict(dict)

        # Regex to capture the datetime and the filename
        # Expected: some_path/<UUID>/<folder_name>/<datetime>/<filename>
        pattern = re.compile(rf"^{re.escape(signoff_prefix)}([^/]+)/(events\.json|labels\.json)$")

        for blob in signoff_blobs:
            match = pattern.match(blob.name)
            if match:
                dt_str, filename = match.groups()
                datetimes_dict[dt_str][filename] = blob

        if datetimes_dict:
            # Standard string sorting easily finds the most recent datetime
            latest_dt = max(datetimes_dict.keys())
            files = datetimes_dict[latest_dt]

            # Ensure both required files are present
            if "events.json" in files and "labels.json" in files:
                latest_blobs_per_uuid[uuid] = files

    print(f"Discovered {len(latest_blobs_per_uuid)} valid {folder_name} cases.")
    return latest_blobs_per_uuid


# ==========================================
# 4. DataFrame Cleaning Function
# ==========================================
def clean_diff(
        df: pd.DataFrame,
        differential_classes: Set[str],
        exclude_unclassified: bool = False,
        include_metadata: bool = True
) -> pd.DataFrame:
    df = df.copy()
    working_diff_classes = set(differential_classes)

    # Track Metadata Columns
    metadata_cols = ["Barcode", "Particles", "Case ID", "Reviewer", "Sign-off Time", "Approved ROIs", "Investigator ROIs"]

    # 1. Combine Eosinophils
    if "Pro Eosinophil" in df.columns:
        if "Eosinophil" not in df.columns:
            df["Eosinophil"] = 0
        df["Eosinophil"] += df["Pro Eosinophil"]
        df.drop(columns=["Pro Eosinophil"], inplace=True)
    working_diff_classes.discard("Pro Eosinophil")
    working_diff_classes.add("Eosinophil")

    # 2. Exclude Unclassified logic
    if exclude_unclassified:
        working_diff_classes.discard("Unclassified")

    # 3. Add Total Nucleated
    existing_diff_cols = [col for col in working_diff_classes if col in df.columns]
    df["Total Nucleated"] = df[existing_diff_cols].sum(axis=1)

    # 4. Convert numbers to Percentages
    totals = df["Total Nucleated"].replace(0, pd.NA)
    for col in existing_diff_cols:
        df[col] = (df[col] / totals * 100).fillna(0)

    # 5. Remove specific columns
    cols_to_remove = [c for c in ["Large Granular Lymphocyte", "Atypical Lymphocyte", "Smudge Cell", "Platelet"] if
                      c in df.columns]
    df.drop(columns=cols_to_remove, inplace=True)

    # 6. Column Ordering Logic
    final_diff_cols = sorted([col for col in working_diff_classes if col in df.columns])

    non_diff_cols = sorted([
        col for col in df.columns
        if col not in final_diff_cols
           and col != "Total Nucleated"
           and col not in metadata_cols
    ])

    # Base configuration: Total Nucleated, Diffs (A-Z), Non-Diffs (A-Z)
    final_order = ["Total Nucleated"] + final_diff_cols + non_diff_cols

    if include_metadata:
        existing_meta = [m for m in metadata_cols if m in df.columns and m != "Barcode"]
        # Add Barcode at the front and everything else at the tail
        if "Barcode" in df.columns:
            final_order = ["Barcode"] + final_order
        final_order = final_order + existing_meta
    else:
        existing_meta = [m for m in metadata_cols if m in df.columns]
        df.drop(columns=existing_meta, inplace=True)

    return df[final_order]

# ==========================================
# 5. Main Orchestration Functions
# ==========================================
def compile_study_results(bucket_name: str, prefix: str, classes_yml_path: str,
                          exclude_unclassified=False, include_metadata=True, with_reclass=False,
                          ) -> pd.DataFrame:
    """
    Main orchestration function to read GCS files directly into memory,
    analyze them, and compile the final Pandas DataFrame.
    """
    scan_blobs = get_latest_signoff_blobs(bucket_name, prefix)

    all_scan_results = {}

    for uuid, files in scan_blobs.items():
        print(f"Processing UUID: {uuid}...")

        # Read file contents directly from GCS to memory (no disk IO)
        events_text = files["events.json"].download_as_text()
        labels_text = files["labels.json"].download_as_text()

        # Parse JSON
        events_data = json.loads(events_text)
        labels_data = json.loads(labels_text)

        # Run analysis and store the resulting dictionary mapped to the UUID
        counts = process_single_scan(events_data, labels_data, with_reclass=with_reclass)
        if counts:  # Ignore if empty (e.g. no valid sign-off found in the file)
            all_scan_results[uuid] = counts

    print("\nCompiling final Pandas DataFrame...")
    # from_dict creates a dataframe where dictionary keys (UUIDs) become the index/rows.
    # Columns are automatically created from the inner dictionary keys (class_ids).
    df = pd.DataFrame.from_dict(all_scan_results, orient='index')

    # Fill missing values with 0 (since not all scans will have every class_id)
    # and convert floats to integers
    numeric_cols = df.select_dtypes(include=['number']).columns
    df[numeric_cols] = df[numeric_cols].fillna(0).astype(int)

    # Rename Columns using mapping dictionary
    print("Loading class dictionary from YAML...")
    class_mapping, differential_classes = load_class_mapping(classes_yml_path)
    df.columns = [class_mapping.get(c, c if isinstance(c, str) else f"Unknown_Class_{c}") for c in df.columns]
    df.index.name = "UUID"

    final_df = clean_diff(df, differential_classes,
                          exclude_unclassified=exclude_unclassified, include_metadata=include_metadata)

    return final_df


def compile_local_events_study(
        local_events_dir: str,
        bucket_name: str,
        gcs_labels_prefix: str,
        classes_yml_path: str,
        exclude_unclassified=False,
        include_metadata=True,
        with_reclass=False,
        labels_dir_name='SCAN_READY'
) -> pd.DataFrame:
    """
    Iterates over local `<UUID>.json` files, fetches their matching `labels.json`
    from GCS, and compiles the final Pandas DataFrame.
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # Ensure standard prefix formatting
    if not gcs_labels_prefix.endswith('/'):
        gcs_labels_prefix += '/'

    all_scan_results = {}

    # Find all local JSON files
    local_files = glob.glob(os.path.join(local_events_dir, "*.json"))
    print(f"Found {len(local_files)} local event files in {local_events_dir}.")

    for local_filepath in local_files:
        # Extract UUID from the filename (e.g., "path/to/123-uuid.json" -> "123-uuid")
        uuid = os.path.splitext(os.path.basename(local_filepath))[0]
        print(f"Processing UUID: {uuid}...")

        # 1. Load the local events file
        with open(local_filepath, 'r') as f:
            events_data = json.load(f)

        # 2. Hunt for the latest labels.json in GCS for this UUID
        # Because we don't know the exact datetime folder without SIGN_OFF,
        # we list all blobs under the UUID prefix and take the latest labels.json
        uuid_prefix = f"{gcs_labels_prefix}{uuid}/{labels_dir_name}/"
        blobs = list(bucket.list_blobs(prefix=uuid_prefix))

        label_blobs = [b for b in blobs if b.name.endswith('labels.json')]
        if not label_blobs:
            print(f"  -> WARNING: No labels.json found in GCS for UUID {uuid}'s {labels_dir_name}. Skipping.")
            continue

        # Standard alphabetical sorting guarantees the latest datetime folder is last
        latest_label_blob = sorted(label_blobs, key=lambda b: b.name)[-1]

        labels_text = latest_label_blob.download_as_text()
        labels_data = json.loads(labels_text)

        # 3. Process data using the no-signoff logic
        counts = process_single_scan_no_signoff(events_data, labels_data, with_reclass=with_reclass)
        if counts:
            all_scan_results[uuid] = counts

    print("\nCompiling final Pandas DataFrame...")

    df = pd.DataFrame.from_dict(all_scan_results, orient='index')

    # Safe Integer casting only for numeric columns (avoids crashing on metadata strings)
    numeric_cols = df.select_dtypes(include=['number']).columns
    df[numeric_cols] = df[numeric_cols].fillna(0).astype(int)

    # Rename Columns using mapping dictionary (ignores strings)
    print("Loading class dictionary from YAML...")
    class_mapping, differential_classes = load_class_mapping(classes_yml_path)
    df.columns = [
        class_mapping.get(c, c if isinstance(c, str) else f"Unknown_Class_{c}")
        for c in df.columns
    ]
    df.index.name = "UUID"

    # Run the final polish
    final_df = clean_diff(
        df,
        differential_classes,
        exclude_unclassified=exclude_unclassified,
        include_metadata=include_metadata
    )

    return final_df

def compile_local_events_study(
        local_events_dir: str,
        bucket_name: str,
        gcs_labels_prefix: str,
        classes_yml_path: str,
        exclude_unclassified=False,
        include_metadata=True,
        with_reclass=False,
        labels_dir_name='SCAN_READY',
        uuids_list=[]
) -> pd.DataFrame:
    """
    Iterates over local `<UUID>.json` files, fetches their matching `labels.json`
    from GCS, and compiles the final Pandas DataFrame.
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # Ensure standard prefix formatting
    if not gcs_labels_prefix.endswith('/'):
        gcs_labels_prefix += '/'

    all_scan_results = {}

    # Find all local JSON files
    local_files = glob.glob(os.path.join(local_events_dir, "*.json"))
    print(f"Found {len(local_files)} local event files in {local_events_dir}.")

    for local_filepath in local_files:
        # Extract UUID from the filename (e.g., "path/to/123-uuid.json" -> "123-uuid")
        uuid = os.path.splitext(os.path.basename(local_filepath))[0]
        if uuids_list and uuid not in uuids_list:
            continue

        print(f"Processing UUID: {uuid}...")

        # 1. Load the local events file
        with open(local_filepath, 'r') as f:
            events_data = json.load(f)

        # 2. Hunt for the latest labels.json in GCS for this UUID
        # Because we don't know the exact datetime folder without SIGN_OFF,
        # we list all blobs under the UUID prefix and take the latest labels.json
        uuid_prefix = f"{gcs_labels_prefix}{uuid}/{labels_dir_name}/"
        blobs = list(bucket.list_blobs(prefix=uuid_prefix))

        label_blobs = [b for b in blobs if b.name.endswith('labels.json')]
        if not label_blobs:
            print(f"  -> WARNING: No labels.json found in GCS for UUID {uuid}'s {labels_dir_name}. Skipping.")
            continue

        # Standard alphabetical sorting guarantees the latest datetime folder is last
        latest_label_blob = sorted(label_blobs, key=lambda b: b.name)[-1]

        labels_text = latest_label_blob.download_as_text()
        labels_data = json.loads(labels_text)

        # 3. Process data using the no-signoff logic
        counts = process_single_scan_no_signoff(events_data, labels_data, with_reclass=with_reclass)
        if counts:
            all_scan_results[uuid] = counts

    print("\nCompiling final Pandas DataFrame...")

    df = pd.DataFrame.from_dict(all_scan_results, orient='index')

    # Safe Integer casting only for numeric columns (avoids crashing on metadata strings)
    numeric_cols = df.select_dtypes(include=['number']).columns
    df[numeric_cols] = df[numeric_cols].fillna(0).astype(int)

    # Rename Columns using mapping dictionary (ignores strings)
    print("Loading class dictionary from YAML...")
    class_mapping, differential_classes = load_class_mapping(classes_yml_path)
    df.columns = [
        class_mapping.get(c, c if isinstance(c, str) else f"Unknown_Class_{c}")
        for c in df.columns
    ]
    df.index.name = "UUID"

    # Run the final polish
    final_df = clean_diff(
        df,
        differential_classes,
        exclude_unclassified=exclude_unclassified,
        include_metadata=include_metadata
    )

    return final_df



if __name__ == "__main__":
    machine = "scopiobox3217"
    BUCKET = "scopio_event_sync_hematology_prod"
    BASE_DIR = "events"
    CLASSES_FILE = "label_classes.yml"  # Path to your local classes mapping file

    not_signed = False  # if True then search does not require sign-off and just takes latest chosen analysis areas

    bucket_name = BUCKET
    PREFIX = fr'{BASE_DIR}/{machine}'



    # Run the pipeline
    if machine == "scopiobox1060" or not_signed:
        LOCAL_EVENTS_DIR = f"./{machine}_events"
        final_df = compile_local_events_study(
            local_events_dir=LOCAL_EVENTS_DIR,
            bucket_name=BUCKET,
            gcs_labels_prefix=PREFIX,
            classes_yml_path=CLASSES_FILE,
            exclude_unclassified=False,
            include_metadata=True,
            with_reclass=True
        )
    else:
        final_df = compile_study_results(BUCKET, PREFIX, CLASSES_FILE,
                                         exclude_unclassified=False, include_metadata=True, with_reclass=False)

    # Display the final results
    print("\nFinal Results DataFrame:")
    print(final_df.head())

    final_df.to_csv(f"{machine}_raw_cbm.csv")





