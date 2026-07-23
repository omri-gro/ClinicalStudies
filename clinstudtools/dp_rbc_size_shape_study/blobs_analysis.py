import os
import glob
import json
import ast
import io
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from google.cloud import storage

# ==========================================
# 1. Configuration & Constants
# ==========================================
MARGIN_OF_ERROR = 30.0  # Pixels for AI vs Reviewer matching
TINY_MARGIN = 10.0  # Pixels for duplicate reviewer click deduplication
EXTENDED_ROI_PIXELS = 15.0  # Pixels to extend the ROI boundary for FN detection

MODEL_VERSION = "pbs-v3.24"
SCAN_RESULTS_PATH = r"S:\talm\cbm_clinical_trial\Final_Run\RGB\pbs-3.25\json_structured"  # Directory containing the scan_id.json files
BBOX_FROM_JSON = True
USE_CACHED_FULL_SCAN = True  # Toggle to False to force re-parsing of JSONs

# If a class is not in this dict, script will use the Reviewer label for both.
MORPHOLOGY_ALIASES = {
    "Bite cell": "Bite",
    "Helmet cell": "Helmet",
    "Blister cell": "Blister",
    "Schistocyte": "Schistocytes",
}

# multi-label limitations on microcytes and macrocytes based on labeling policy
APPLY_SHAPE_EXCLUSIONS_FOR_SIZE = True
MACROCYTE_EXCLUSIONS = {
    "Sickle", "Bite", "Blister", "Spherocyte", "Acanthocyte", "Helmet",
    "Teardrop", "Schistocytes", "Poikilocyte", "Elliptocyte"
}
MICROCYTE_EXCLUSIONS = MACROCYTE_EXCLUSIONS.union({"Stomatocyte"})

# Translates the JSON keys to standard names so exclusions trigger correctly
JSON_TO_STD_MORPH = {
    "sickle": "Sickle",
    "bite": "Bite",
    "blister": "Blister",
    "spherocytes": "Spherocyte",
    "spur": "Acanthocyte",
    "acanthocytes": "Acanthocyte",
    "tear_drop": "Teardrop",
    "schistocytes": "Schistocyte",
    "poikilocytes": "Poikilocyte",
    "elliptocytes": "Elliptocyte",
    "stomatocytes": "Stomatocyte",
    "helmet": "Helmet",
    "burr": "Echinocyte"
}


# ==========================================
# 2. Data Loading
# ==========================================
def load_mapping(filepath, bucket=None):
    """Loads the tasks mapping CSV/Excel from local disk or GCS."""
    if bucket:
        blob = bucket.blob(filepath)
        content = blob.download_as_bytes()
        if filepath.endswith('.csv'):
            return pd.read_csv(io.BytesIO(content))
        return pd.read_excel(io.BytesIO(content))
    else:
        if filepath.endswith('.csv'):
            return pd.read_csv(filepath)
        return pd.read_excel(filepath)


def load_rois_from_json(json_path, bucket=None):
    """Parses ROIs from the scans JSON."""
    if bucket:
        blob = bucket.blob(json_path)
        data = json.loads(blob.download_as_text())
    else:
        with open(json_path, 'r') as f:
            data = json.load(f)

    roi_list = []
    for scan in data.get('scans', []):
        scan_id = scan.get('scan_id')
        res_mm_pix = scan.get('resolution_mm_pix')
        for roi in scan.get('rois', []):
            if roi.get('roi_deleted', False):
                continue
            # Filter out null tasks
            tasks = [t for t in roi.get('tasks', []) if t is not None]
            roi_list.append({
                'scan_id': scan_id, 'roi_id': roi.get('roi_id'), 'roi_name': roi.get('name'),
                'x': roi.get('x'), 'y': roi.get('y'), 'w': roi.get('w'), 'h': roi.get('h'),
                'resolution_mm_pix': res_mm_pix,
                'tasks': tasks
            })
    return pd.DataFrame(roi_list)


def load_blobs(directory_path, bucket=None):
    """Loads all blob CSVs without upfront evaluation to save time/memory."""
    df_list = []

    if bucket:
        # directory_path acts as a prefix in GCS
        blobs = bucket.list_blobs(prefix=directory_path)
        for blob in blobs:
            if "_blobs_" in blob.name and blob.name.endswith(".csv"):
                content = blob.download_as_bytes()
                df_list.append(pd.read_csv(io.BytesIO(content), low_memory=False))
    else:
        all_files = glob.glob(os.path.join(directory_path, "*_blobs_*.csv"))
        for file in all_files:
            df_list.append(pd.read_csv(file, low_memory=False))

    if not df_list:
        raise FileNotFoundError("No blobs CSV files found.")

    blobs_df = pd.concat(df_list, ignore_index=True)
    # Fill NAs to safely allow string operations later
    blobs_df['ai_morphologies'] = blobs_df['ai_morphologies'].fillna('[]')
    return blobs_df


def load_scan_json_blobs(json_path, base_class):
    """Extracts base_class detections from the local scan_id.json into a compatible DataFrame."""
    if not os.path.exists(json_path):
        print(f"\033[31mNo {json_path} file!\033[0m")
        return pd.DataFrame(), None

    with open(json_path, 'r') as f:
        data = json.load(f)

    # Convert um to mm to match calculate_size_based_morphologies logic
    res_mm_pix = data.get("um_per_pixel", 0.1659) / 1000.0

    rows = []
    for d in data.get("detections", []):
        cls_info = d.get("classification") or {}
        if cls_info.get("name", "").upper() == base_class.upper():
            bounds = d.get("bounds", {})
            x = bounds.get("x", 0)
            y = bounds.get("y", 0)
            w = bounds.get("width", 0)
            h = bounds.get("height", 0)

            morphs = d.get("morphologies")
            if morphs:
                # Map internal lowercase json keys to standard names
                base_morph_list = [JSON_TO_STD_MORPH.get(k, k) for k in morphs.keys()]
            else:
                base_morph_list = []

            rows.append({
                "bounding_box": f"[{x}, {y}, {w}, {h}]",
                "base_morph_list": base_morph_list
            })
    return pd.DataFrame(rows), res_mm_pix


# ==========================================
# 3. Geometry & Logic Helpers
# ==========================================
def get_blobs_in_roi(blobs_df, roi_series, expand_by=0):
    """Filters blobs that fall within the ROI bounding box."""
    if blobs_df.empty or roi_series is None:
        return blobs_df

    x_min, y_min = roi_series['x'] - expand_by, roi_series['y'] - expand_by
    x_max, y_max = roi_series['x'] + roi_series['w'] + expand_by, roi_series['y'] + roi_series['h'] + expand_by

    mask = (
            (blobs_df['center_x'] >= x_min) & (blobs_df['center_x'] <= x_max) &
            (blobs_df['center_y'] >= y_min) & (blobs_df['center_y'] <= y_max)
    )
    return blobs_df[mask].copy()


def check_roi_intersection(roi1, roi2):
    """Checks if two ROIs overlap geometrically."""
    if roi1 is None or roi2 is None:
        return False

    left1, right1 = roi1['x'], roi1['x'] + roi1['w']
    top1, bottom1 = roi1['y'], roi1['y'] + roi1['h']

    left2, right2 = roi2['x'], roi2['x'] + roi2['w']
    top2, bottom2 = roi2['y'], roi2['y'] + roi2['h']

    # If one rectangle is on left side of other, or above other, they don't intersect
    if right1 <= left2 or left1 >= right2 or bottom1 <= top2 or top1 >= bottom2:
        return False
    return True


def safe_eval_morphologies(morph_series):
    """Safely evaluates the morphology strings only for the requested subset."""

    def _eval(val):
        try:
            parsed = ast.literal_eval(val)
            return parsed if isinstance(parsed, list) else []
        except:
            return []
    return morph_series.apply(_eval)


def calculate_size_based_morphologies(bbox_str, res_mm_pix, mode, task_row, ai_morph_list):
    """Computes AI morphology based on bounding box physical size for size-based modes."""
    if pd.isna(bbox_str) or res_mm_pix is None:
        return []
    try:
        bbox = ast.literal_eval(bbox_str)
        w, h = bbox[2], bbox[3]

        if mode == "PLT_SIZE":
            size_px = max(w, h)
        elif mode == "RBC_SIZE":
            size_px = (w + h) / 2.0
        else:
            return []

        # Convert pixels to microns (resolution is mm/pix. 1mm = 1000 microns)
        size_um = size_px * res_mm_pix * 1000.0

        morphs = []
        if mode == "PLT_SIZE":
            giant_def = task_row.get("PLT giant definition")
            large_def = task_row.get("PLT large definition")
            if pd.notna(giant_def) and size_um > giant_def:
                morphs.append("PLT giant")
            elif pd.notna(large_def) and size_um > large_def:
                morphs.append("PLT large")

        elif mode == "RBC_SIZE":
            macro_def = task_row.get("RBC macrocyte definition")
            micro_def = task_row.get("RBC microcyte definition")

            # Check multi-label exclusions
            has_macro_exclusion = APPLY_SHAPE_EXCLUSIONS_FOR_SIZE and any(shape in ai_morph_list for shape in MACROCYTE_EXCLUSIONS)
            has_micro_exclusion = APPLY_SHAPE_EXCLUSIONS_FOR_SIZE and any(shape in ai_morph_list for shape in MICROCYTE_EXCLUSIONS)

            if pd.notna(macro_def) and size_um > macro_def:
                if not has_macro_exclusion:
                    morphs.append("RBC macrocyte")
            elif pd.notna(micro_def) and size_um < micro_def:
                if not has_micro_exclusion:
                    morphs.append("RBC microcyte")

        return morphs
    except (ValueError, SyntaxError, TypeError, IndexError):
        return []


def deduplicate_reviewer_blobs(rev_blobs, tiny_margin):
    """Removes duplicate blobs clicked too close together by the same reviewer."""
    if len(rev_blobs) <= 1:
        return rev_blobs, 0
    coords = rev_blobs[['center_x', 'center_y']].values
    dist_matrix = cdist(coords, coords)
    keep_indices = []
    dropped_count = 0
    for i in range(len(coords)):
        if any(dist_matrix[i, j] < tiny_margin for j in keep_indices):
            dropped_count += 1
        else:
            keep_indices.append(i)
    return rev_blobs.iloc[keep_indices], dropped_count


def calculate_confusion_metrics(rev_blobs, ai_blobs_ext, target_morphs, margin):
    """Matches reviewer marks with AI blobs to calculate confusion metrics."""
    # Ensure target_morphs is always a list for consistent iteration
    if isinstance(target_morphs, str):
        target_morphs = [target_morphs]

    metrics = {'TPs': 0, 'FPs': 0, 'FNs': 0, 'empty_blobs': 0}

    if rev_blobs.empty:
        metrics['FNs'] = sum([any(tm in m for tm in target_morphs) for m in ai_blobs_ext['ai_morph_list']])
        return metrics

    rev_coords = rev_blobs[['center_x', 'center_y']].values
    ai_coords = ai_blobs_ext[['center_x', 'center_y']].values if not ai_blobs_ext.empty else np.array([])
    matched_ai_indices = set()

    for rev_coord in rev_coords:
        if len(ai_coords) == 0:
            metrics['empty_blobs'] += 1
            continue

        distances = cdist([rev_coord], ai_coords)[0]
        valid_matches = np.where(distances <= margin)[0]
        valid_matches = [idx for idx in valid_matches if idx not in matched_ai_indices]

        if not valid_matches:
            metrics['empty_blobs'] += 1
        else:
            best_match_idx = valid_matches[np.argmin(distances[valid_matches])]
            matched_ai_indices.add(best_match_idx)

            # Check if ANY of the valid target morphs are in this matched AI blob
            if any(tm in ai_blobs_ext.iloc[best_match_idx]['ai_morph_list'] for tm in target_morphs):
                metrics['TPs'] += 1
            else:
                metrics['FPs'] += 1

    # FNs: AI detected ANY of the classes in the extended ROI, but it wasn't matched
    for j in range(len(ai_blobs_ext)):
        if j not in matched_ai_indices and target_morphs in ai_blobs_ext.iloc[j]['ai_morph_list']:
            metrics['FNs'] += 1

    return metrics


# ==========================================
# 4. Core Pipeline
# ==========================================
def process_study(mapping_df, rois_df, blobs_df, morphologies, required_roi_names=None, cached_scan_cells=None):
    results = []
    all_scan_cells_list = []
    processed_scans = set()   # Prevents duplicating cells if a scan has multiple tasks
    if required_roi_names:
        rois_df = rois_df[rois_df['roi_name'].isin(required_roi_names) | (rois_df['roi_name'] == BASE_CLASS)]

    # Ensure ai_version has no NaNs for safe string comparison
    blobs_df['ai_version'] = blobs_df['ai_version'].fillna('unknown')

    for _, task_row in mapping_df.iterrows():
        scan_id = task_row['scan_id']
        task_id = task_row['task_id']
        reviewer_name = task_row.get('Reviewer Name', task_row.get('Reviewer', 'Unknown'))

        scan_blobs = blobs_df[blobs_df['scan_id'] == scan_id]
        scan_rois = rois_df[rois_df['scan_id'] == scan_id]

        # --- Analysis 1: Full Scan AI ---
        if STUDY_MODE in ["PLT_SIZE", "RBC_SIZE"] and BBOX_FROM_JSON:
            # Route A: Compute dynamically from local JSON bounding boxes
            # Check if we have this scan in our cached CSV
            if cached_scan_cells is not None and scan_id in cached_scan_cells['scan_id'].values:
                base_ai_blobs = cached_scan_cells[cached_scan_cells['scan_id'] == scan_id].copy()
                # Grab the resolution directly from the ROI definitions since we bypassed the JSON
                json_res_mm_pix = scan_rois.iloc[0]['resolution_mm_pix'] if not scan_rois.empty else 0.1659 / 1000.0
            else:
                json_path = os.path.join(SCAN_RESULTS_PATH, scan_id, "results.json").replace("\\", "/")
                base_ai_blobs, json_res_mm_pix = load_scan_json_blobs(json_path, JSON_CLASS)

                if not base_ai_blobs.empty and scan_id not in processed_scans:
                    cache_df = base_ai_blobs.copy()
                    cache_df['scan_id'] = scan_id
                    # Calculate the numeric sizes for our regression database
                    cache_df['Cell_Size_Max_um'] = cache_df['bounding_box'].apply(
                        lambda b: max(ast.literal_eval(b)[2],
                                      ast.literal_eval(b)[3]) * json_res_mm_pix * 1000.0 if pd.notna(b) else None
                    )
                    cache_df['Cell_Size_Mean_um'] = cache_df['bounding_box'].apply(
                        lambda b: ((ast.literal_eval(b)[2] + ast.literal_eval(b)[3]) / 2.0) * json_res_mm_pix * 1000.0 if pd.notna(b) else None
                    )
                    all_scan_cells_list.append(cache_df[
                                                   ['scan_id', 'bounding_box', 'Cell_Size_Max_um', 'Cell_Size_Mean_um', 'base_morph_list']])
                    processed_scans.add(scan_id)

            cells_in_scan = len(base_ai_blobs)

            if not base_ai_blobs.empty:
                base_ai_blobs['ai_morph_list'] = base_ai_blobs.apply(
                    lambda row: calculate_size_based_morphologies(
                        row['bounding_box'], json_res_mm_pix, STUDY_MODE, task_row, row['base_morph_list']
                    ), axis=1
                )
            else:
                base_ai_blobs['ai_morph_list'] = pd.Series(dtype=object)

        else:
            # Route B: Use standard Blobs_df string morphologies
            # AI ROI processing
            base_ai_roi_df = scan_rois[(scan_rois['roi_name'] == BASE_CLASS)]
            base_ai_roi = base_ai_roi_df.iloc[0] if not base_ai_roi_df.empty else None

            # Grab model's blobs from model-chosen ROI
            base_ai_blobs = scan_blobs[(scan_blobs['label_name'] == BASE_CLASS) & (scan_blobs['ai_version'] != MODEL_VERSION)]
            if base_ai_roi is not None:
                base_ai_blobs = get_blobs_in_roi(base_ai_blobs, base_ai_roi)
            cells_in_scan = len(base_ai_blobs)


        # --- Analysis 2: Reviewers' ROIs ---
        rev_rois = scan_rois[scan_rois['tasks'].apply(lambda t: task_id in t)]
        rev_roi = None
        if not rev_rois.empty:

            # cases where multiple ROIs were drawn within the same task
            if len(rev_rois) > 1:
                # Find all ROIs that actually contain blobs
                rois_with_blobs = []
                for _, r in rev_rois.iterrows():
                    if not scan_blobs[scan_blobs['roi_id'] == r['roi_id']].empty:
                        rois_with_blobs.append(r)

                # Check how many populated ROIs we found
                if len(rois_with_blobs) > 1:
                    roi_names = [r['roi_name'] for r in rois_with_blobs]

                    # Prioritize the ROI whose name ends with '1'
                    rois_ending_in_1 = [r for r in rois_with_blobs if str(r['roi_name']).strip().endswith('1')]

                    if rois_ending_in_1:
                        rev_roi = rois_ending_in_1[0]
                        chosen_name = rev_roi['roi_name']
                        print(f"WARNING: Multiple populated ROIs found for scan {scan_id} (Task: {task_id}). "
                              f"ROI names: {roi_names}. Proceeding with '{chosen_name}' as it ends with '1'.")
                    else:
                        rev_roi = rois_with_blobs[0]
                        chosen_name = rev_roi['roi_name']
                        print(f"WARNING: Multiple populated ROIs found for scan {scan_id} (Task: {task_id}). "
                              f"ROI names: {roi_names}. Proceeding with the first one: '{chosen_name}'.")

                    rev_roi = rois_with_blobs[0]
                elif len(rois_with_blobs) == 1:
                    rev_roi = rois_with_blobs[0]
                else:
                    # Fallback: all multiple ROIs are empty
                    rev_roi = rev_rois.iloc[0]

            else:
                rev_roi = rev_rois.iloc[0]

        # --- Analysis 3: BLOB-TO-BLOB AI ('pbs-v3.24' Model)
        # Identify valid AI ROI IDs for this scan (name ends with 'model', case-insensitive)
        valid_model_rois = scan_rois[scan_rois['roi_name'].fillna('').str.lower().str.endswith('model') | scan_rois['roi_name'].fillna('').str.lower().str.endswith('study')]
        valid_model_roi_ids = valid_model_rois['roi_id'].tolist()

        # Check if the AI actually covered the Reviewer's area
        ai_missing_in_roi = False
        if rev_roi is not None:
            has_coverage = any(check_roi_intersection(rev_roi, row) for _, row in valid_model_rois.iterrows())
            if not has_coverage:
                ai_missing_in_roi = True
                print(f"Notice: AI did not run on the area chosen by reviewer in scan {scan_id} (Task: {task_id}).")

        # Extract blobs specific to the new model
        model_ai_blobs = scan_blobs[
            (scan_blobs['label_name'] == BASE_CLASS) &
            (scan_blobs['ai_version'] == MODEL_VERSION) &
            (scan_blobs['roi_id'].isin(valid_model_roi_ids))
        ]

        # Extract Extended Reviewer AI Blobs ONLY ONCE per scan/task
        ai_in_ext_roi = pd.DataFrame()
        if rev_roi is not None:
            # Spatially filter the model's blobs using the Reviewer's ROI
            ai_in_ext_roi = get_blobs_in_roi(model_ai_blobs, rev_roi, expand_by=EXTENDED_ROI_PIXELS)
            if ai_in_ext_roi.empty:
                ai_in_ext_roi['base_morph_list'] = pd.Series(dtype=object)
                ai_in_ext_roi['ai_morph_list'] = pd.Series(dtype=object)
            else:
                ai_in_ext_roi['base_morph_list'] = safe_eval_morphologies(ai_in_ext_roi['ai_morphologies'])
                if STUDY_MODE in ["PLT_SIZE", "RBC_SIZE"]:
                    # Use Size-Based Calculations
                    res_mm_pix = rev_roi['resolution_mm_pix']
                    ai_in_ext_roi['ai_morph_list'] = ai_in_ext_roi.apply(
                        lambda row: calculate_size_based_morphologies(
                            row['bounding_box'], res_mm_pix, STUDY_MODE, task_row, row['base_morph_list']
                        ), axis=1
                    )
                else:
                    # Standard shape evaluation
                    ai_in_ext_roi['ai_morph_list'] = ai_in_ext_roi['base_morph_list']

        for morphology in morphologies:
            # Define target labels (handling standard vs synthetic combined classes)
            if morphology == "Helmet&Schisto":
                target_rev_labels = ["Schistocyte", "Helmet cell"]
                target_ai_labels = [
                    MORPHOLOGY_ALIASES.get("Schistocyte", "Schistocyte"),
                    MORPHOLOGY_ALIASES.get("Helmet cell", "Helmet cell")
                ]
            elif morphology == "PLT large&giant":
                target_rev_labels = ["PLT large", "PLT giant"]
                target_ai_labels = ["PLT large", "PLT giant"]
            else:
                target_rev_labels = [morphology]
                target_ai_labels = [MORPHOLOGY_ALIASES.get(morphology, morphology)]

            res = {
                'Morphology': morphology, 'SampleID': task_row.get('SampleID'),
                'Site': task_row.get('Site'), 'Internal Count': task_row.get('Internal Count'),
                'Reviewer': task_row.get('Reviewer'), 'Reviewer Name': reviewer_name,
                'Task Name': task_row.get('Task Name'),
                'New task': task_row.get('New task name'), 'scan_id': scan_id, 'task_id': task_id,
                'roi_id': rev_roi['roi_id'] if rev_roi is not None else None,
                'roi_name': rev_roi['roi_name'] if rev_roi is not None else None,
                'roi_size': (rev_roi['w'] * rev_roi['h']) if rev_roi is not None else None,
                'cells_in_scan': cells_in_scan
            }

            # Full scan stats
            if STUDY_MODE in ["PLT_SIZE", "RBC_SIZE"] and BBOX_FROM_JSON:
                # Route A (JSON): Search through the actual Python lists we generated
                class_in_scan = sum(any(label in m for label in target_ai_labels) for m in
                                    base_ai_blobs['ai_morph_list']) if cells_in_scan else 0
            else:
                # Route B (CSV): Fast regex search on the raw string column
                ai_regex = '|'.join([f'"{label}"' for label in target_ai_labels])
                class_in_scan = base_ai_blobs['ai_morphologies'].str.contains(ai_regex, regex=True, na=False).sum() if cells_in_scan else 0

            res['class_in_scan'] = class_in_scan
            res['percent_in_scan_ai'] = round((100 * class_in_scan / cells_in_scan), 2) if cells_in_scan else None

            if rev_roi is None:
                res.update({'cells_in_roi': None, 'class_in_roi_ai': None, 'class_in_roi_rev': None,
                            'percent_in_roi_ai': None, 'percent_in_roi_rev': None, 'TPs': None,
                            'FPs': None, 'FNs': None, 'empty_blobs': None, 'double_blobs': None})
            else:
                # Reviewer explicitly marked blobs (matches ANY of the target reviewer labels)
                rev_all_blobs = scan_blobs[
                    (scan_blobs['roi_id'] == rev_roi['roi_id']) &
                    (scan_blobs['label_name'].isin(target_rev_labels))
                    ]
                rev_blobs_dedup, double_blobs = deduplicate_reviewer_blobs(rev_all_blobs, TINY_MARGIN)
                class_in_roi_rev = len(rev_blobs_dedup)

                if ai_missing_in_roi:
                    # Output None for all AI-dependent ROI columns
                    res.update({
                        'cells_in_roi': None, 'class_in_roi_ai': None,
                        'class_in_roi_rev': class_in_roi_rev, 'double_blobs': double_blobs,
                        'percent_in_roi_ai': None, 'percent_in_roi_rev': None,
                        'TPs': None, 'FPs': None, 'FNs': None, 'empty_blobs': None
                    })
                else:
                    # Filter down the already-evaluated extended AI blobs to just the strict ROI boundaries
                    ai_in_roi = get_blobs_in_roi(ai_in_ext_roi, rev_roi, expand_by=0)
                    cells_in_roi = len(ai_in_roi)

                    # AI ROI logic (matches ANY of the target AI labels)
                    class_in_roi_ai = sum([any(label in m for label in target_ai_labels) for m in
                                           ai_in_roi['ai_morph_list']]) if cells_in_roi else 0

                    res.update({
                        'cells_in_roi': cells_in_roi, 'class_in_roi_ai': class_in_roi_ai,
                        'class_in_roi_rev': class_in_roi_rev, 'double_blobs': double_blobs,
                        'percent_in_roi_ai': round((100 * class_in_roi_ai / cells_in_roi), 2) if cells_in_roi else None,
                        'percent_in_roi_rev': round((100 * class_in_roi_rev / cells_in_roi), 2) if cells_in_roi else None
                    })

                    # Confusion metrics use the pre-evaluated `ai_in_ext_roi` subset
                    metrics = calculate_confusion_metrics(rev_blobs_dedup, ai_in_ext_roi, target_ai_labels, MARGIN_OF_ERROR)
                    res.update(metrics)

            results.append(res)

    # Save or update the master database of all scan cells
    if all_scan_cells_list:
        new_cells_df = pd.concat(all_scan_cells_list, ignore_index=True)

        # If we had existing cached cells, merge the new ones with them
        if cached_scan_cells is not None:
            new_cells_df = pd.concat([cached_scan_cells, new_cells_df], ignore_index=True)

        new_cells_df.to_csv(rf"side_results/full_scan_cells_{STUDY_MODE}.csv", index=False)
        print(f"Updated full_scan_cells_{STUDY_MODE}.csv with new scan data.")

    cols_order = ['Morphology', 'SampleID', 'Site', 'Internal Count', 'Reviewer', 'Reviewer Name', 'Task Name', 'New task',
                  'scan_id', 'task_id', 'roi_id', 'roi_name', 'roi_size', 'cells_in_scan', 'cells_in_roi',
                  'class_in_scan', 'class_in_roi_ai', 'class_in_roi_rev', 'percent_in_scan_ai',
                  'percent_in_roi_ai', 'percent_in_roi_rev', 'TPs', 'FPs', 'FNs', 'empty_blobs', 'double_blobs']
    return pd.DataFrame(results)[cols_order]


if __name__ == "__main__":
    STUDY_MODE = "PLT_SIZE"  # RBC_Shape, RBC_SIZE or PLT_SIZE
    USE_GCS = True  # toggle this to False is blobs csv and scans json are in local directory
    GCS_BUCKET_NAME = "scopio_labeling_tool_datasets_eur"

    # If GCS, DATA_DIR is the prefix path in the bucket. If local, it's the folder path.
    DATA_DIR = "PLT_Size_Study/2026-07-16_16:35:07.541866+00:00/"
    # DATA_DIR = "RBC_Shape_Study/2026-07-16_07:02:49.736033+00:00/"
    # DATA_DIR = "RBC_Size_Study/2026-07-16_16:28:34.192300+00:00/"
    MAPPING_DIR = "./mapping"  # location of tasks mapping

    MAPPING_FILE = os.path.join(MAPPING_DIR, f"{STUDY_MODE}_tasks_mapping.csv")
    JSON_FILE = os.path.join(DATA_DIR, "scans_0.json").replace("\\", "/")

    # Initialize storage
    bucket = None
    if USE_GCS:
        print(f"Connecting to GCS Bucket: {GCS_BUCKET_NAME}...")
        client = storage.Client()
        bucket = client.bucket(GCS_BUCKET_NAME)

    # define cells to look for
    BASE_CLASS = "PLT" if STUDY_MODE == "PLT_SIZE" else "RBC"
    JSON_CLASS = "platelet" if STUDY_MODE == "PLT_SIZE" else "rbc"
    if STUDY_MODE == "PLT_SIZE":
        REQUESTED_MORPHOLOGIES = ["PLT large", "PLT giant", "PLT large&giant"]
    elif STUDY_MODE == "RBC_SIZE":
        REQUESTED_MORPHOLOGIES = ["RBC macrocyte", "RBC microcyte"]
    else:
        REQUESTED_MORPHOLOGIES = ["Bite cell", "Helmet cell", "Spherocyte", "Schistocyte", "Blister cell", "Helmet&Schisto"]

    # load data
    mapping_df = load_mapping(MAPPING_FILE)  # mapping does not sit in the bucket, so don't include the bucket argument
    rois_df = load_rois_from_json(JSON_FILE, bucket=bucket)
    blobs_df = load_blobs(DATA_DIR, bucket=bucket)

    blobs_df = blobs_df[blobs_df['roi_id'].isin(rois_df['roi_id'])]  # Purge any blobs belonging to ROIs that were marked as deleted in the JSON

    # Load cached full-scan cells if available
    cached_df = None
    cache_file = rf"side_results/full_scan_cells_{STUDY_MODE}.csv"
    if USE_CACHED_FULL_SCAN and os.path.exists(cache_file):
        print(f"Loading cached scan cells from {cache_file}...")
        cached_df = pd.read_csv(cache_file)
        # Convert string representations of lists back to Python lists upfront to save loop time
        cached_df['base_morph_list'] = cached_df['base_morph_list'].apply(
            lambda x: ast.literal_eval(x) if pd.notna(x) else [])

    # Process Study
    output_df = process_study(mapping_df, rois_df, blobs_df, REQUESTED_MORPHOLOGIES)
    output_df.to_csv(rf"results/{STUDY_MODE}_results.csv", index=False)
