import os
import ast
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from google.cloud import storage

# Import the data loading and geometry tools we already built
from blobs_analysis import (
    load_mapping,
    load_rois_from_json,
    load_blobs,
    get_blobs_in_roi,
    deduplicate_reviewer_blobs,
    MARGIN_OF_ERROR,
    TINY_MARGIN,
    MODEL_VERSION
)


def calculate_cell_sizes(bbox_str, res_mm_pix):
    """Parses the bounding box and returns both Max and Mean sizes in microns."""
    if pd.isna(bbox_str) or res_mm_pix is None:
        return None, None
    try:
        bbox = ast.literal_eval(bbox_str)
        w, h = bbox[2], bbox[3]

        # Convert pixels to microns (resolution is mm/pix. 1mm = 1000 microns)
        size_max_um = max(w, h) * res_mm_pix * 1000.0
        size_mean_um = ((w + h) / 2.0) * res_mm_pix * 1000.0

        return size_max_um, size_mean_um
    except (ValueError, SyntaxError, TypeError, IndexError):
        return None, None


def extract_cell_level_data(mapping_df, rois_df, blobs_df, study_mode):
    """
    Builds a cell-level DataFrame where every row is an AI-detected candidate cell,
    mapped to the specific action the reviewer took on it.
    """
    cell_data = []
    base_class = "PLT" if study_mode == "PLT_SIZE" else "RBC"

    # Identify valid model ROIs
    valid_model_rois = rois_df[rois_df['roi_name'].fillna('').str.lower().str.endswith('model') |
                               rois_df['roi_name'].fillna('').str.lower().str.endswith('study')]
    valid_model_roi_ids = valid_model_rois['roi_id'].tolist()

    # Isolate AI model blobs
    model_ai_blobs = blobs_df[
        (blobs_df['label_name'] == base_class) &
        (blobs_df['ai_version'] == MODEL_VERSION) &
        (blobs_df['roi_id'].isin(valid_model_roi_ids))
        ]

    for _, task_row in mapping_df.iterrows():
        scan_id = task_row['scan_id']
        task_id = task_row['task_id']
        reviewer_name = task_row.get('Reviewer Name', task_row.get('Reviewer', 'Unknown'))

        scan_blobs = blobs_df[blobs_df['scan_id'] == scan_id]
        scan_rois = rois_df[rois_df['scan_id'] == scan_id]

        # --- 1. Identify the Reviewer's ROI ---
        rev_rois = scan_rois[scan_rois['tasks'].apply(lambda t: task_id in t)]
        if rev_rois.empty:
            continue

        rev_roi = None
        if len(rev_rois) > 1:
            rois_with_blobs = [r for _, r in rev_rois.iterrows() if
                               not scan_blobs[scan_blobs['roi_id'] == r['roi_id']].empty]
            if rois_with_blobs:
                rois_ending_in_1 = [r for r in rois_with_blobs if str(r['roi_name']).strip().endswith('1')]
                rev_roi = rois_ending_in_1[0] if rois_ending_in_1 else rois_with_blobs[0]
            else:
                rev_roi = rev_rois.iloc[0]
        else:
            rev_roi = rev_rois.iloc[0]

        res_mm_pix = rev_roi['resolution_mm_pix']

        # --- 2. Get AI Cells in Reviewer's Area ---
        # We don't extend the ROI here because we want strictly the cells the reviewer was evaluating
        ai_in_roi = get_blobs_in_roi(model_ai_blobs, rev_roi, expand_by=0)
        if ai_in_roi.empty:
            continue

        # Extract coordinates for matching
        ai_coords = ai_in_roi[['center_x', 'center_y']].values

        # Initialize all AI cells as "Ignored"
        ai_actions = ["Ignored"] * len(ai_in_roi)

        # --- 3. Get Reviewer's Clicks ---
        rev_all_blobs = scan_blobs[scan_blobs['roi_id'] == rev_roi['roi_id']]
        if not rev_all_blobs.empty:
            rev_blobs_dedup, _ = deduplicate_reviewer_blobs(rev_all_blobs, TINY_MARGIN)
            rev_coords = rev_blobs_dedup[['center_x', 'center_y']].values

            # --- 4. Match Clicks to AI Cells ---
            # We use the same greedy matching logic from the confusion matrix
            matched_ai_indices = set()
            for i, rev_coord in enumerate(rev_coords):
                distances = cdist([rev_coord], ai_coords)[0]
                valid_matches = np.where(distances <= MARGIN_OF_ERROR)[0]
                valid_matches = [idx for idx in valid_matches if idx not in matched_ai_indices]

                if valid_matches:
                    best_match_idx = valid_matches[np.argmin(distances[valid_matches])]
                    matched_ai_indices.add(best_match_idx)

                    # Log the reviewer's specific class label for this cell
                    rev_label = rev_blobs_dedup.iloc[i]['label_name']
                    ai_actions[best_match_idx] = f"Marked as {rev_label}"

        # --- 5. Compile Cell-Level Rows ---
        for j, (_, ai_blob) in enumerate(ai_in_roi.iterrows()):
            size_max, size_mean = calculate_cell_sizes(ai_blob.get('bounding_box'), res_mm_pix)

            cell_data.append({
                'scan_id': scan_id,
                'task_id': task_id,
                'roi_id': rev_roi['roi_id'],
                'Reviewer_Name': reviewer_name,
                'blob_id': ai_blob.get('blob_id'),
                'Cell_Size_Max_um': size_max,
                'Cell_Size_Mean_um': size_mean,
                'Reviewer_Action': ai_actions[j]
            })

    return pd.DataFrame(cell_data)


def sweep_empirical_thresholds(df, size_col, positive_actions, min_val=0.0, max_val=20.0, step=0.1,
                               direction='greater'):
    """
    Sweeps a threshold to find the optimal size (in microns) that maximizes the F1-score.
    'direction' can be 'greater' (e.g., Macrocytes) or 'less' (e.g., Microcytes).
    """
    if df.empty:
        return None, 0.0, 0.0, 0.0

    y_true = df['Reviewer_Action'].isin(positive_actions).values
    sizes = df[size_col].values

    if not y_true.any():
        return None, 0.0, 0.0, 0.0

    thresholds = np.arange(min_val, max_val + step, step)

    best_f1 = 0.0
    best_th = None
    best_prec = 0.0
    best_rec = 0.0

    for th in thresholds:
        # NEW: Apply directional logic
        if direction == 'greater':
            y_pred = sizes > th
        elif direction == 'less':
            y_pred = sizes < th
        else:
            raise ValueError("Direction must be 'greater' or 'less'.")

        tp = (y_true & y_pred).sum()
        fp = (~y_true & y_pred).sum()
        fn = (y_true & ~y_pred).sum()

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        if f1 > best_f1:
            best_f1 = f1
            best_th = th
            best_prec = precision
            best_rec = recall

    return best_th, best_f1, best_prec, best_rec


def analyze_reviewer_thresholds(cell_df, target_actions, min_val=0.0, max_val=20.0, step=0.1, direction='greater'):
    """
    Groups data by reviewer and runs the directional sweeping logic.
    """
    results = []

    for reviewer, rev_df in cell_df.groupby('Reviewer_Name'):
        # Evaluate using Cell_Size_Max_um
        max_th, max_f1, max_prec, max_rec = sweep_empirical_thresholds(
            rev_df, 'Cell_Size_Max_um', target_actions, min_val, max_val, step, direction
        )

        # Evaluate using Cell_Size_Mean_um
        mean_th, mean_f1, mean_prec, mean_rec = sweep_empirical_thresholds(
            rev_df, 'Cell_Size_Mean_um', target_actions, min_val, max_val, step, direction
        )

        results.append({
            'Reviewer_Name': reviewer,
            'Total_Cells_Evaluated': len(rev_df),
            'Target_Marks_By_Reviewer': rev_df['Reviewer_Action'].isin(target_actions).sum(),
            'Sweep_Direction': direction,

            # Results using Max(w, h)
            'Optimal_Max_Threshold_um': round(max_th, 2) if max_th is not None else None,
            'Max_Calc_F1': round(max_f1, 4),
            'Max_Calc_Precision': round(max_prec, 4),
            'Max_Calc_Recall': round(max_rec, 4),

            # Results using Mean(w, h)
            'Optimal_Mean_Threshold_um': round(mean_th, 2) if mean_th is not None else None,
            'Mean_Calc_F1': round(mean_f1, 4),
            'Mean_Calc_Precision': round(mean_prec, 4),
            'Mean_Calc_Recall': round(mean_rec, 4)
        })

    return pd.DataFrame(results)


if __name__ == "__main__":
    # Example execution
    STUDY_MODE = "PLT_SIZE"  # RBC_Shape, RBC_SIZE or PLT_SIZE
    USE_GCS = True  # toggle this to False is blobs csv and scans json are in local directory
    GCS_BUCKET_NAME = "scopio_labeling_tool_datasets_eur"

    # If GCS, DATA_DIR is the prefix path in the bucket. If local, it's the folder path.
    DATA_DIR = "PLT_Size_Study/2026-07-16_16:35:07.541866+00:00/"
    MAPPING_DIR = "./mapping"  # location of tasks mapping

    MAPPING_FILE = os.path.join(MAPPING_DIR, f"{STUDY_MODE}_tasks_mapping.csv")
    JSON_FILE = os.path.join(DATA_DIR, "scans_0.json").replace("\\", "/")

    # 1. Load Data using main script tools
    # Initialize storage
    bucket = None
    if USE_GCS:
        print(f"Connecting to GCS Bucket: {GCS_BUCKET_NAME}...")
        client = storage.Client()
        bucket = client.bucket(GCS_BUCKET_NAME)

    mapping_df = load_mapping(MAPPING_FILE)  # mapping does not sit in the bucket, so don't include the bucket argument
    rois_df = load_rois_from_json(JSON_FILE, bucket=bucket)
    blobs_df = load_blobs(DATA_DIR, bucket=bucket)
    blobs_df = blobs_df[blobs_df['roi_id'].isin(rois_df['roi_id'])]  # Purge any blobs belonging to ROIs that were marked as deleted in the JSON

    # 2. Extract Cell Data
    cell_df = extract_cell_level_data(mapping_df, rois_df, blobs_df, STUDY_MODE)

    # Save cell-level data
    output_name = rf"side_results/behavioral_cell_data_{STUDY_MODE}.csv"
    cell_df.to_csv(output_name, index=False)
    print(f"Extraction complete! {len(cell_df)} candidate cells saved to {output_name}.")


    # 3. Run Empirical Threshold Analysis
    # Define what constitutes a "Positive" action for the threshold we are testing.
    # E.g., To find the threshold for "Large", we count if they marked it Large OR Giant.
    sweep_configs = []

    if STUDY_MODE == "PLT_SIZE":
        sweep_configs = [
            {
                'name': 'PLT Large (includes Giant)',
                'targets': ["Marked as PLT large", "Marked as PLT giant"],
                'direction': 'greater',
                'min': 2.0, 'max': 15.0
            }
        ]
    elif STUDY_MODE == "RBC_SIZE":
        sweep_configs = [
            {
                'name': 'RBC Macrocyte',
                'targets': ["Marked as RBC macrocyte"],
                'direction': 'greater',
                'min': 5.0, 'max': 15.0
            },
            {
                'name': 'RBC Microcyte',
                'targets': ["Marked as RBC microcyte"],
                'direction': 'less',
                'min': 2.0, 'max': 10.0
            }
        ]

    all_results = []
    for config in sweep_configs:
        res_df = analyze_reviewer_thresholds(
            cell_df,
            target_actions=config['targets'],
            min_val=config['min'],
            max_val=config['max'],
            step=0.1,
            direction=config['direction']
        )
        # Insert the class name at the front so we know which row is which
        res_df.insert(0, 'Target_Class', config['name'])
        all_results.append(res_df)

    final_threshold_df = pd.concat(all_results, ignore_index=True)

    # Save the threshold analysis
    threshold_output_name = rf"side_results/reviewer_empirical_thresholds_{STUDY_MODE}.csv"
    final_threshold_df.to_csv(threshold_output_name, index=False)
    print(f"Threshold sweeping complete! Results saved to {threshold_output_name}.")
