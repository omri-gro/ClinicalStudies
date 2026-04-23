from __future__ import annotations

import json
import os
from typing import Dict, List, Set, Tuple
from collections import Counter


def load_json(filepath: str) -> dict | list:
    """Helper function to load a JSON file."""
    with open(filepath, 'r') as file:
        return json.load(file)


def get_approved_rois_at_signoff(events: List[dict]) -> Tuple[Set[str], Dict[str, Tuple[int, int, int, int]]]:
    """
    Parses the events log to find the final state of approved ROIs at sign-off.

    Returns:
        A tuple containing:
        - A set of approved `roi_id` strings.
        - A dictionary mapping `roi_id` to its bounding box (x_min, y_min, x_max, y_max).
    """
    # 1. Find the index of the final (latest) SIGN_OFF event.
    # Since the file is reverse-chronological, the first match is the final sign-off.
    sign_off_idx = next(
        (i for i, e in enumerate(events)
         if e.get("event_type") == "SESSION_PROGRESS_UPDATED"
         and e.get("payload", {}).get("step") == "SIGN_OFF"),
        None
    )

    if sign_off_idx is None:
        raise ValueError("No SIGN_OFF event found in the events log.")

    # 2. Find the final state of selected ROIs prior to (or at) the sign-off.
    approved_roi_ids: Set[str] = set()
    for event in events[sign_off_idx:]:
        if event.get("event_type") == "SELECTED_REGIONS_UPDATED":
            selections = event.get("payload", {}).get("selection", [])
            # Extract all ROIs that the reviewer explicitly left as "SELECTED"
            approved_roi_ids = {
                roi["roi_id"] for roi in selections
                if roi.get("selection_status") == "SELECTED"
            }
            break  # Found the state at the time of sign-off

    # 3. Extract the physical coordinates for these selected ROIs.
    # We scan all events. The first time we encounter an ANALYSIS_REGION_UPDATED for
    # a given roi_id, it represents its latest coordinate state.
    roi_bounds: Dict[str, Tuple[int, int, int, int]] = {}
    for event in events:
        if event.get("event_type") == "ANALYSIS_REGION_UPDATED":
            payload = event.get("payload", {})
            roi_id = payload.get("roi_id")
            bounds = payload.get("bounds")

            # Record bounds if it's an approved ROI and we haven't mapped it yet
            if roi_id in approved_roi_ids and roi_id not in roi_bounds and bounds:
                x_min = bounds["xleft"]
                y_min = bounds["ytop"]
                x_max = x_min + bounds["width"]
                y_max = y_min + bounds["height"]
                roi_bounds[roi_id] = (x_min, y_min, x_max, y_max)

    return approved_roi_ids, roi_bounds


def calculate_approved_classifications(labels_data: dict, roi_bounds: Dict[str, Tuple[int, int, int, int]]) -> Dict[
    int, int]:
    """
    Calculates the number of cells per model class_id inside the approved ROIs.
    """
    class_counts = Counter()

    # Extract just the bounding box values for faster iteration over cells
    approved_boxes = list(roi_bounds.values())

    # Iterate over all types in labels.json (e.g., 'nucleated', 'megakaryocyte')
    for cell_type, cells in labels_data.items():
        for cell in cells:
            cx = cell.get("center_x")
            cy = cell.get("center_y")
            class_id = cell.get("extra", {}).get("class_id")

            if cx is None or cy is None or class_id is None:
                continue

            # Check if the cell's center coordinates fall within ANY of the approved ROIs
            in_approved_area = any(
                x_min <= cx <= x_max and y_min <= cy <= y_max
                for x_min, y_min, x_max, y_max in approved_boxes
            )

            if in_approved_area:
                class_counts[class_id] += 1

    return dict(class_counts)


def analyze_case(labels_file: str, events_file: str):
    """Main pipeline execution for a single case/scan."""
    # Load JSON files
    events_data = load_json(events_file)
    labels_data = load_json(labels_file)

    # Extract Region definitions and states
    approved_roi_ids, roi_bounds = get_approved_rois_at_signoff(events_data)

    # Count model classifications using only those within the approved ROIs
    classification_counts = calculate_approved_classifications(labels_data, roi_bounds)

    # Print or log the results
    print(f"Total Approved ROIs: {len(approved_roi_ids)}")
    print(f"Classification Counts inside Approved ROIs:\n{json.dumps(classification_counts, indent=2)}")

    return classification_counts


if __name__ == "__main__":
    jsons_dir = 'script_prep_raw_files-OHSUR051TW'
    labels_path = os.path.join(jsons_dir, "labels.json")
    events_path = os.path.join(jsons_dir, "events.json")

    # Example execution based on the files provided
    analyze_case(labels_path, events_path)

