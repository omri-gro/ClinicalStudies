"""
version of results_aggregation specific for parameters included in 2026 CBM submission
all parameters are collected, but the ones to appear in the submission appear first
"""
import os
import json
import re
import traceback
from datetime import datetime
from collections import Counter
from tqdm import tqdm
import pandas as pd

# ==========================================
# CONFIGURATION
# ==========================================
# Expected parameters for final CBM submission
SUBMISSION_COLUMNS = [
    'Barcode', 'Scan UUID', 'Site', 'Creation_Date', 'Total WBC', 'Total RBC',
    'Blast', 'Monocyte', 'Basophil', 'Eosinophil', 'nRBC',
    'Segmented Neutrophil', 'Band Neutrophil', 'Metamyelocyte',
    'Myelocyte', 'Promyelocyte', 'Lymphocyte', 'LGL',
    'Reactive Lymphocyte', 'Abnormal Lymphocyte', 'Smudge',
    'Plasma cell', 'Hypochromia', 'Sickle cells', 'Stomatocytes', 'Target cells',
    'Tear drop cells', 'Polychromasia', 'Parasites', 'Spherocytes',
    'Schistocytes', 'Macrocytes', 'Microcytes', 'Large Platelets',
    'Platelet Clumps', 'Platelet Satellitism', 'Platelets Estimate',
    'Unclassified WBC'
]

# Classification groupings for relative percentage calculations
WBC_TYPES = {
    'Segmented Neutrophil', 'Band Neutrophil', 'Metamyelocyte', 'Myelocyte',
    'Promyelocyte', 'Blast', 'Monocyte', 'Basophil', 'Eosinophil',
    'Lymphocyte', 'Large Granular Lymphocyte', 'Atypical Lymphocyte',
    'Plasma Cell', 'Aberrant Lymphocyte', 'Tier Two Aberrant Lymphocyte',
    'Hairy Cell', 'Sezary Cell', 'Unclassified WBC'
}

# Togglable sets for calculation methods
PER_100_WBC_TYPES = {'Smudge Cell', 'Normoblast', 'Dohle Bodies', 'Pelger Cell', 'Auer Rods'}
# PER_10_FOV_TYPES = {'Smudge Cell'}
PER_10_FOV_TYPES = {}

# Renaming dictionary mapping raw JSON terms to final requested terms
JSON_TO_FINAL_MAP = {
    'Normoblast': 'nRBC',
    'Atypical Lymphocyte': 'Reactive Lymphocyte',
    'Large Granular Lymphocyte': 'LGL',
    'Smudge Cell': 'Smudge',
    'Plasma Cell': 'Plasma cell',
    'hypochromatic': 'Hypochromia',
    'sickle': 'Sickle cells',
    'target': 'Target cells',
    'tear_drop': 'Tear drop cells',
    'polychromatic': 'Polychromasia',
    'micro_organisms': 'Parasites',
    'stomatocytes': 'Stomatocytes',
    'spherocytes': 'Spherocytes',
    'schistocytes': 'Schistocytes'
}

# Size-based Morphology Definitions
SIZE_CRITERIA = {
    'Macrocytes': {
        'target_class': 'rbc',
        'threshold': 8.5,
        'operator': 'greater',
        'exclude_morphologies': {'sickle', 'bite', 'spherocytes', 'ovalocytes',
                                 'blister', 'spur', 'poikilocytes', 'schistocytes', 'elliptocytes', 'helmet'}
    },
    'Microcytes': {
        'target_class': 'rbc',
        'threshold': 6.5,
        'operator': 'less',
        'exclude_morphologies': {'sickle', 'bite', 'spherocytes', 'ovalocytes', 'stomatocytes',
                                 'blister', 'spur', 'poikilocytes', 'schistocytes', 'elliptocytes', 'helmet'}
    },
    'Large Platelets': {
        'target_class': 'platelet',
        'threshold': 4.0,
        'operator': 'greater',
        'exclude_morphologies': {}
    }
}

SIZE_CRITERIA = {
    'Macrocytes': {
        'target_class': 'rbc',
        'threshold': 9.5,
        'operator': 'greater',
        'exclude_morphologies': {'sickle', 'bite', 'spherocytes', 'ovalocytes',
                                 'blister', 'spur', 'poikilocytes', 'schistocytes', 'elliptocytes', 'helmet'}
    },
    'Microcytes': {
        'target_class': 'rbc',
        'threshold': 6.0,
        'operator': 'less',
        'exclude_morphologies': {'sickle', 'bite', 'spherocytes', 'ovalocytes', 'stomatocytes',
                                 'blister', 'spur', 'poikilocytes', 'schistocytes', 'elliptocytes', 'helmet'}
    },
    'Large Platelets': {
        'target_class': 'platelet',
        'threshold': 6.0,
        'operator': 'greater',
        'exclude_morphologies': {}
    }
}


# ==========================================
# PARSER LOGIC
# ==========================================

def get_datetime_from_log(log_path, time_format="%Y-%m-%d %H:%M:%S"):
    """Extracts the completion datetime from the scan's log file."""
    if not os.path.exists(log_path):
        return ""

    try:
        with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
            first_line = f.readline()
        match = re.search(r"\[(.*?)\]", first_line)
        if match:
            raw_dt = match.group(1)
            return datetime.strptime(raw_dt, "%d-%m-%Y %H:%M:%S.%f").strftime(time_format)
    except (ValueError, AttributeError, IOError):
        pass
    return ""


def extract_raw_counts(detections):
    """Single-pass iteration to tally classes, morphologies, and size metrics."""
    counts = Counter()
    total_wbc, total_rbc, total_plt = 0, 0, 0

    for det in detections:
        classification = det.get('classification') or {}
        class_name = classification.get('name', '')
        counts[class_name] += 1

        # Track base totals
        if class_name in WBC_TYPES:
            total_wbc += 1
        elif class_name == 'rbc':
            total_rbc += 1
        elif class_name == 'platelet':
            total_plt += 1

        # Track morphologies
        morphs = set((det.get('morphologies') or {}).keys())
        for m in morphs:
            counts[m] += 1

        # Evaluate size-based criteria
        attributes = det.get('attributes') or {}
        diameter = attributes.get('diameter_um')
        if diameter is not None:
            for criteria_name, config in SIZE_CRITERIA.items():
                if class_name == config['target_class']:
                    if not morphs.intersection(config['exclude_morphologies']):
                        op = config['operator']
                        if (op == 'greater' and diameter > config['threshold']) or \
                                (op == 'less' and diameter < config['threshold']):
                            counts[criteria_name] += 1

    return counts, total_wbc, total_rbc, total_plt


def calculate_wbc_metrics(counts, total_wbc):
    """Calculates all percentages relative to the Total WBC count."""
    metrics = {}
    wbc_denominator = max(total_wbc, 1)  # Prevent ZeroDivisionError

    # Calculate standard WBC parts and extra WBC morphologies
    for cell in WBC_TYPES.union(PER_100_WBC_TYPES):
        std_name = JSON_TO_FINAL_MAP.get(cell, cell)
        metrics[std_name] = round(counts[cell] / wbc_denominator * 100, 2)

    # Calculate sum parameter
    abnormal_count = (counts.get('Aberrant Lymphocyte', 0) +
                      counts.get('Tier Two Aberrant Lymphocyte', 0) +
                      counts.get('Sezary Cell', 0))

    # Handle Hairy Cell logic
    hairy_cell_count = counts.get('Hairy Cell', 0)

    if CONDITIONAL_HAIRY_CELL:
        hairy_cell_pct = (hairy_cell_count / wbc_denominator) * 100
        if hairy_cell_pct > 2.0:
            abnormal_count += hairy_cell_count
        else:
            # Reassign to regular Lymphocytes and recalculate the metric
            lymph_count = counts.get('Lymphocyte', 0) + hairy_cell_count
            metrics['Lymphocyte'] = round(lymph_count / wbc_denominator * 100, 2)
    else:
        # Default behavior: always abnormal
        abnormal_count += hairy_cell_count

    metrics['Abnormal Lymphocyte'] = round(abnormal_count / wbc_denominator * 100, 2)

    return metrics


def calculate_fov_metrics(counts, jd):
    """Calculates parameters based on the physical High Power Field (HPF) area scanned."""
    metrics = {}

    # Safely extract um_per_pixel (fallback to 0.2016 from standard scans)
    um_per_pixel = jd.get('um_per_pixel', 0.2016)
    pix2hpf = (um_per_pixel ** 2) / (195 ** 2)

    # Safely sum up the areas of all regions by cell type
    areas_hpf = {}
    for region in jd.get('regions', []):
        rtype = region.get('cell_type')
        bounds = region.get('bounds') or {}
        area_px = bounds.get('width', 0) * bounds.get('height', 0)
        areas_hpf[rtype] = areas_hpf.get(rtype, 0) + (area_px * pix2hpf)

    # We default to the 'wbc' region area for these elements. Prevent ZeroDivisionError.
    wbc_area_hpf = max(areas_hpf.get('wbc', 1), 1)

    metrics['WBC Analysis FOVs'] = round(wbc_area_hpf, 2)

    for cell in PER_10_FOV_TYPES:
        std_name = JSON_TO_FINAL_MAP.get(cell, cell)
        count_val = counts.get(cell, 0)
        metrics[std_name] = round((count_val / wbc_area_hpf) * 10, 2)

    return metrics


def parse_json(jd, scan_uuid, site):
    """Orchestrates parsing for a single scan JSON."""
    cv_info = jd.get('cv_info') or {}
    profile = jd.get('profile') or {}

    result = {
        'Scan UUID': scan_uuid,
        'Site': site,
        'Barcode': jd.get('barcode', ''),
        'CVI Version': cv_info.get('cvi_version', ''),
        'Scan Mode': profile.get('scan_profile', '')
    }

    # 1. Tally single-pass metrics
    counts, total_wbc, total_rbc, total_plt = extract_raw_counts(jd.get('detections', []))

    result['Total WBC'] = total_wbc
    result['Total RBC'] = total_rbc

    # 2. Extract specific RBC traits populated in "distributions" by the system
    for dist in jd.get('distributions', []):
        if dist['cell_type'] == 'rbc':
            for obs in dist['observations']:
                raw_name = obs['name']
                final_name = JSON_TO_FINAL_MAP.get(raw_name, raw_name)
                result[final_name] = round(obs['percentage'], 2)

    # 3. Add Calculated Relatives & FOV Metrics
    result.update(calculate_wbc_metrics(counts, total_wbc))
    result.update(calculate_fov_metrics(counts, jd))

    # 4. Add Absolute Counts / Plt Metrics
    result['Platelets Estimate'] = counts['platelet']
    result['Platelet Clumps'] = counts['plt_clumps'] + counts['plt_clump'] + counts['platelet_clumps']
    result['Platelet Satellitism'] = counts['Satellitisms']
    result['Dirt'] = counts['Dirt']  # Keeping Dirt for the full_results

    # 5. Add custom size-based parameters
    result['Macrocytes'] = round(counts['Macrocytes'] / max(total_rbc, 1) * 100, 2)
    result['Microcytes'] = round(counts['Microcytes'] / max(total_rbc, 1) * 100, 2)
    result['Large Platelets'] = round(counts['Large Platelets'] / max(total_plt, 1) * 100, 2)

    # 6. Dump all remaining raw frequencies into the dictionary for 'full_results'
    for k, v in counts.items():
        if k not in result and k not in JSON_TO_FINAL_MAP:
            result[f"raw_count_{k}"] = v

    return result


def process_site(site_name, site_path, all_results):
    """Processes all JSON files for a given site directory with a progress bar."""
    # Check for 'pbs' subdirectory; if not present, assume the site_path is the search directory
    search_dir = os.path.join(site_path, 'pbs') if os.path.exists(os.path.join(site_path, 'pbs')) else site_path
    logs_dir = os.path.join(site_path, 'logs')

    if not os.path.exists(search_dir):
        print(f"Warning: Directory not found -> {search_dir}")
        return

    tasks = []

    # Map out files based on the chosen directory structure
    if USE_SUBDIR_STRUCTURE:
        for item_name in os.listdir(search_dir):
            item_path = os.path.join(search_dir, item_name)

            if os.path.isdir(item_path):
                scan_uuid = item_name
                json_path = os.path.join(item_path, 'results.json')

                if os.path.exists(json_path):
                    # Check for log file in the scan's subdirectory first, fallback to logs folder
                    log_path = os.path.join(item_path, f'{scan_uuid}.log')
                    if not os.path.exists(log_path):
                        log_path = os.path.join(logs_dir, f'{scan_uuid}.log')

                    tasks.append((scan_uuid, json_path, log_path))
    else:
        for filename in os.listdir(search_dir):
            if filename.endswith('.json'):
                scan_uuid = os.path.splitext(filename)[0]
                json_path = os.path.join(search_dir, filename)
                log_path = os.path.join(logs_dir, f'{scan_uuid}.log')

                tasks.append((scan_uuid, json_path, log_path))

    if not tasks:
        print(f"No valid JSON files found in -> {search_dir}")
        return

    # Process mapped tasks with the progress bar
    for scan_uuid, file_path, log_path in tqdm(tasks, desc=f"Processing {site_name}"):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                jd = json.load(f)

            scan_data = parse_json(jd, scan_uuid, site_name)
            scan_data['Creation_Date'] = get_datetime_from_log(log_path)

            all_results.append(scan_data)

        except Exception as e:
            tqdm.write(f"\nError parsing scan ID {scan_uuid} at site {site_name}: {e}")
            tqdm.write(traceback.format_exc())
            continue

# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == '__main__':
    run_name = "rnr_run"

    # False: The directory contains flat files named {scan_uuid}.json
    # True: The directory contains subdirectories named {scan_uuid}, each holding a 'results.json'
    USE_SUBDIR_STRUCTURE = False

    # True: Hairy Cell > 2% -> Abnormal Lymphocyte | Hairy Cell <= 2% -> regular Lymphocyte
    # False: Hairy Cell is always counted as Abnormal Lymphocyte
    CONDITIONAL_HAIRY_CELL = True

    parent_dir = r"C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing"
    output_dir = r"C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\CBM_Study\results\scans_analysis"

    # SITES can contain either folder names (which will be appended to PARENT_DIR)
    # OR full absolute paths to specific directories.
    sites = ["sb1024", "sb1108", "sb1114", "sb1127", "sb1132", "sb1134", "sb3058", "sb3130", "sb3184", "sb3334"]
    # sites = [
    #     r"S:\talm\cbm_clinical_trial\Final_Run\RGB\pbs-3.25\json_structured",
    #     r"S:\talm\cbm_clinical_trial\Final_Run\Amber\pbs-3.25\json_structured"
    # ]

    all_results = []
    os.makedirs(output_dir, exist_ok=True)

    for site_entry in sites:
        # Determine if the entry is a full absolute path or a relative folder name
        if os.path.isabs(site_entry):
            site_name = os.path.basename(site_entry)
            site_path = site_entry
        else:
            site_name = site_entry
            site_path = os.path.join(parent_dir, site_entry)

        process_site(site_name, site_path, all_results)

    if all_results:
        # Construct and clean master DataFrame
        df_full = pd.DataFrame(all_results)
        sort_cols = [c for c in ['Barcode', 'Site', 'Creation_Date'] if c in df_full.columns]
        df_full.sort_values(by=sort_cols, inplace=True)

        # Extract strict subset for CBM submission
        df_submission = pd.DataFrame()
        for col in SUBMISSION_COLUMNS:
            df_submission[col] = df_full[col] if col in df_full.columns else ""

        # Save to disk
        full_output_path = os.path.join(output_dir, f"{run_name}_full_results_changed_thresh.csv")
        cbm_output_path = os.path.join(output_dir, f"{run_name}_cbm_results_thresh.csv")

        df_full.to_csv(full_output_path, index=False)
        df_submission.to_csv(cbm_output_path, index=False)

        print(f"Data collection complete! Found {len(all_results)} scans.")
        print(f"Full results saved to: {full_output_path}")
        print(f"CBM specific results saved to: {cbm_output_path}")
    else:
        print("No JSON files were successfully parsed.")

