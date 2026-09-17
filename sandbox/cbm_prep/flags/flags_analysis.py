import pandas as pd
import numpy as np
from sklearn.utils import resample

# ==========================================
# 1. TOGGLES & CONFIGURATION
# ==========================================
# FILE_PATH = 'Flag Sensitivity Analysis - comp - final.csv'
FILE_PATH = 'Flag Sensitivity Analysis - comp - promyelo1%.csv'

# Toggle Variables
CALC_95_CI = False  # If True, calculates 95% CI using bootstrapping
SHOW_CONFUSION_MATRIX = True  # If True, displays TP, TN, FP, FN counts
INCLUDE_NORMALS = True  # If True, includes 'Normal' samples as True Negatives

# Column Definitions
SITE_COL = 'Site'
NORMAL_COL = 'Normal'
TECH_COL = 'Tech Flag'
LEUKO_COL = 'Extreme Leukopenia'

# Clinical Flags (CBM Test Method)
CBM_CLINICAL_FLAGS = [
    'Blast.1', 'Plasma', 'Promyelocyte.1', 'Abnormal Lym',
    'Unclass', 'Parasites', 'Schistocytes',
    'Schisto > 0.5\nCBM',
    'Any Clinical Flag'
]

# Mapping CBM Flags to their respective Reference Arms
# (Plasma, Unclass, and Any Clinical Flag are omitted here as they have no matching ref arm in BL-BS)
FLAG_MAPPINGS = {
    'Extreme Leukopenia': ['Extreme Leukopenia (WBC Auto=<2)'],
    'Plasma': ['Plasma cells > 1\nmanual'],
    'Promyelocyte.1': ['Promyelocyte > 1\nmanual', 'Promyelocyte > 1\nmanual or OMR'],
    'Blast.1': ['Blast >1\nmanual', 'Blast >1 manual or OMR'],
    'Abnormal Lym': ['Abnorm lym >1\nmanual'],
    'Parasites': ['Parasites > 0\nClV'],
    'Schistocytes': ['Schisto > 1\nDP'],
    'Schisto > 0.5 CBM': ['Schisto > 0.5 DP'],
    'Unclass': ['Any other clinical flag', 'Any other WBC clinical flag', 'Any WBC clinical flag manual'],
    'Any Clinical Flag': ['Any clinical flag manual in any referene method']
}


# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================
def parse_bool(val):
    if pd.isna(val):
        return np.nan
    val_str = str(val).strip().lower()
    if val_str in ['1', '1.0', 'true', 'yes']: return True
    if val_str in ['0', '0.0', 'false', 'no']: return False
    return np.nan


def calc_metrics(tp, tn, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    ppv = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    return sens, spec, ppv


def bootstrap_ci(df, cbm_col, ref_col, n_iterations=1000):
    sens_list, spec_list, ppv_list = [], [], []
    for _ in range(n_iterations):
        sample = resample(df)
        tp = ((sample[cbm_col] == True) & (sample[ref_col] == True)).sum()
        tn = ((sample[cbm_col] == False) & (sample[ref_col] == False)).sum()
        fp = ((sample[cbm_col] == True) & (sample[ref_col] == False)).sum()
        fn = ((sample[cbm_col] == False) & (sample[ref_col] == True)).sum()

        se, sp, pv = calc_metrics(tp, tn, fp, fn)
        sens_list.append(se)
        spec_list.append(sp)
        ppv_list.append(pv)

    def get_ci(m_list):
        valid = [m for m in m_list if not np.isnan(m)]
        if len(valid) < n_iterations * 0.5:
            return np.nan, np.nan
        return np.percentile(valid, 2.5), np.percentile(valid, 97.5)

    return get_ci(sens_list), get_ci(spec_list), get_ci(ppv_list)


def print_rates(df, condition_col, cohort_name):
    # Overall
    total = len(df)
    count = df[condition_col].sum()
    pct = (count / total) * 100 if total > 0 else 0
    print(f"\n--- {condition_col} Rates ({cohort_name}) ---")
    print(f"Overall: {count}/{total} ({pct:.2f}%)")

    # Per Site
    sites = sorted(df[SITE_COL].dropna().unique())
    for site in sites:
        sub = df[df[SITE_COL] == site]
        s_total = len(sub)
        s_count = sub[condition_col].sum()
        s_pct = (s_count / s_total) * 100 if s_total > 0 else 0
        print(f"  {site}: {s_count}/{s_total} ({s_pct:.2f}%)")


def format_metric(val, ci_tuple=None):
    if pd.isna(val): return "-"
    base = f"{val * 100:.2f}%"
    if ci_tuple and CALC_95_CI:
        lower, upper = ci_tuple
        base += f" [{lower * 100:.2f}-{upper * 100:.2f}]"
    return base


# ==========================================
# 3. DATA LOADING & PREPROCESSING
# ==========================================
print("Loading data...")
df = pd.read_csv(FILE_PATH, low_memory=False)

# Convert all relevant columns to boolean mappings
bool_cols = [TECH_COL, LEUKO_COL, NORMAL_COL] + CBM_CLINICAL_FLAGS
for refs in FLAG_MAPPINGS.values():
    bool_cols.extend(refs)

for col in set(bool_cols):
    if col in df.columns:
        df[col] = df[col].apply(parse_bool)

# Define Hierarchical Cohorts
df_all = df.copy()
df_no_tech = df_all[df_all[TECH_COL] == False].copy()
df_clin_eligible = df_no_tech[df_no_tech[LEUKO_COL] == False].copy()

# ==========================================
# 4. PREVALENCE RATES
# ==========================================
print_rates(df_all, TECH_COL, "All Scans")
print_rates(df_no_tech, LEUKO_COL, "Excluding Tech Flags")

print("\n--- Clinical Flag Rates (Excl. Tech & Leukopenia) ---")
print(f"Total Eligible Scans: {len(df_clin_eligible)}")
for flag in CBM_CLINICAL_FLAGS:
    if flag in df_clin_eligible.columns:
        count = df_clin_eligible[flag].sum()
        pct = (count / len(df_clin_eligible)) * 100 if len(df_clin_eligible) > 0 else 0
        print(f"  {flag}: {count} ({pct:.2f}%)")

# ==========================================
# 5. SENSITIVITY, SPECIFICITY, PPV
# ==========================================
print("\n" + "=" * 50)
print("PERFORMANCE METRICS (CBM vs. REFERENCE)")
print("=" * 50)

for cbm_col, ref_cols in FLAG_MAPPINGS.items():
    if cbm_col not in df.columns: continue

    for ref_col in ref_cols:
        if ref_col not in df.columns: continue

        print(f"\n>>> Analyzing: '{cbm_col}' vs '{ref_col}'")

        # 1. Base cohort selection based on hierarchy
        if cbm_col == 'Extreme Leukopenia':
            eval_df = df_no_tech.copy()
        else:
            eval_df = df_clin_eligible.copy()


        # 2. Inject Normals or Drop missing references
        def resolve_ref(row):
            if INCLUDE_NORMALS and row[NORMAL_COL] == True:
                return False  # Normal cohort is strictly True Negative
            return row[ref_col]


        eval_df['resolved_ref'] = eval_df.apply(resolve_ref, axis=1)
        eval_df = eval_df.dropna(subset=[cbm_col, 'resolved_ref'])


        # 3. Calculation Engine
        def evaluate_cohort(sub_df, label):
            if len(sub_df) == 0: return

            tp = ((sub_df[cbm_col] == True) & (sub_df['resolved_ref'] == True)).sum()
            tn = ((sub_df[cbm_col] == False) & (sub_df['resolved_ref'] == False)).sum()
            fp = ((sub_df[cbm_col] == True) & (sub_df['resolved_ref'] == False)).sum()
            fn = ((sub_df[cbm_col] == False) & (sub_df['resolved_ref'] == True)).sum()

            se, sp, pv = calc_metrics(tp, tn, fp, fn)

            se_ci, sp_ci, pv_ci = None, None, None
            if CALC_95_CI:
                se_ci, sp_ci, pv_ci = bootstrap_ci(sub_df, cbm_col, 'resolved_ref')

            out = f"{label:<12} | N={len(sub_df):<4} | Sens: {format_metric(se, se_ci):<18} | Spec: {format_metric(sp, sp_ci):<18} | PPV: {format_metric(pv, pv_ci):<18}"
            if SHOW_CONFUSION_MATRIX:
                out += f" | (TP:{tp} TN:{tn} FP:{fp} FN:{fn})"
            print(out)


        # 4. Execute Overall & Per Site
        evaluate_cohort(eval_df, "OVERALL")
        sites = sorted(eval_df[SITE_COL].dropna().unique())
        for site in sites:
            evaluate_cohort(eval_df[eval_df[SITE_COL] == site], f"  {site}")

# ==========================================
# 6. UNCLASS FLAG DEEP DIVE
# ==========================================
print("\n" + "=" * 50)
print("UNCLASS FLAG DEEP DIVE")
print("=" * 50)

if 'Unclass' in df_clin_eligible.columns:
    # Isolate only the scans that received the 'Unclass' flag
    df_unclass = df_clin_eligible[df_clin_eligible['Unclass'] == True].copy()
    total_unclass = len(df_unclass)

    print(f"Total scans with 'Unclass' flag (excluding Tech & Leukopenia): {total_unclass}\n")

    if total_unclass > 0:
        # 1. Received another CBM clinical flag
        other_clin_cols = [c for c in CBM_CLINICAL_FLAGS if
                           c in df_unclass.columns and c not in ['Unclass', 'Any Clinical Flag']]
        df_unclass['Another_Clin_Flag'] = df_unclass[other_clin_cols].any(axis=1)

        count_another = df_unclass['Another_Clin_Flag'].sum()
        pct_another = (count_another / total_unclass) * 100
        print(f"> Received ANOTHER clinical flag (CBM):")
        print(f"  {count_another} / {total_unclass} ({pct_another:.2f}%)")

        # 2. Received another WBC-related clinical flag (CBM)
        wbc_clin_cols = ['Blast.1', 'Promyelocyte.1', 'Abnormal Lym']
        wbc_clin_cols = [c for c in wbc_clin_cols if c in df_unclass.columns]
        df_unclass['Another_WBC_Flag'] = df_unclass[wbc_clin_cols].any(axis=1)

        count_wbc = df_unclass['Another_WBC_Flag'].sum()
        pct_wbc = (count_wbc / total_unclass) * 100
        print(f"\n> Received a WBC-RELATED clinical flag (CBM):")
        print(f"  {count_wbc} / {total_unclass} ({pct_wbc:.2f}%)")

        # 3. Would have received a flag by Reference Arms
        ref_flag_cols = []
        for ref_list in FLAG_MAPPINGS.values():
            ref_flag_cols.extend(ref_list)
        ref_flag_cols = [c for c in set(ref_flag_cols) if c in df_unclass.columns]

        df_unclass['Any_Ref_Flag'] = df_unclass[ref_flag_cols].any(axis=1)
        count_ref = df_unclass['Any_Ref_Flag'].sum()
        pct_ref = (count_ref / total_unclass) * 100
        print(f"\n> Triggered ANY clinical/leukopenia flag in the Reference Arms:")
        print(f"  {count_ref} / {total_unclass} ({pct_ref:.2f}%)")

        # Breakdown of Reference Arm Triggers
        print("\n  Breakdown of Reference Flags for these 'Unclass' samples:")
        for ref_col in ref_flag_cols:
            count = df_unclass[ref_col].sum()
            if count > 0:
                clean_ref_name = ref_col.replace('\n', ' ')
                pct = (count / total_unclass) * 100
                print(f"    - {clean_ref_name}: {count} ({pct:.2f}%)")

        # ------------------------------------------
        # Custom Composite Reference Flags
        # ------------------------------------------

        # A) Any manual/OMR clinical flag
        def check_manual_omr_flag(row):
            # Blast
            if parse_bool(row.get('Blasts>1|manual')) == True or parse_bool(row.get('Blasts>1|OMR')) == True: return True
            # Promyelocyte
            if parse_bool(row.get('Promyelo>1|manual')) == True or parse_bool(row.get('Promyelo>1|OMR')) == True: return True
            # Abnorm Lym
            if parse_bool(row.get('Abnorm lym >1|manual')) == True or parse_bool(row.get('Abnorm lym >1|OMR')) == True: return True
            # Plasma
            if parse_bool(row.get('Plasma>1|OMR')) == True or parse_bool(row.get('Plasma>1|OMR')) == True: return True

            return False


        df_unclass['Any_manual_OMR_clinical_flag'] = df_unclass.apply(check_manual_omr_flag, axis=1)
        count_man_omr = df_unclass['Any_manual_OMR_clinical_flag'].sum()
        pct_man_omr = (count_man_omr / total_unclass) * 100 if total_unclass else 0

        print(f"\n> Triggered ANY manual/OMR clinical flag:")
        print(f"  {count_man_omr} / {total_unclass} ({pct_man_omr:.2f}%)")


        # B) Any reference arm clinical flag
        def check_any_ref_flag(row):
            # If it already triggered the manual/OMR composite, it triggers this one
            if row['Any_manual_OMR_clinical_flag'] == True: return True
            # Schisto DP
            if parse_bool(row.get('Schisto > 0.5\nDP')) == True: return True
            # Parasites ClV
            if parse_bool(row.get('Parasites > 0.85\nClV')) == True: return True
            # Parasites OMR
            if parse_bool(row.get('Parasites|OMR')) == True: return True
            return False

    df_unclass['Any_ref_arm_clinical_flag'] = df_unclass.apply(check_any_ref_flag, axis=1)
    count_any_ref = df_unclass['Any_ref_arm_clinical_flag'].sum()
    pct_any_ref = (count_any_ref / total_unclass) * 100 if total_unclass else 0

    print(f"\n> Triggered ANY reference arm clinical flag:")
    print(f"  {count_any_ref} / {total_unclass} ({pct_any_ref:.2f}%)")

    # C) Any reference arm clinical flag or leukopenia
    def check_any_ref_or_leuk_flag(row):
        # If it already triggered the manual/OMR composite, it triggers this one
        if row['Any_ref_arm_clinical_flag'] == True: return True
        # Leukopenia Auto CBC
        if parse_bool(row.get('Extreme Leukopenia (WBC Auto=<2)')) == True: return True
        return False

    df_unclass['Any_ref_arm_flag'] = df_unclass.apply(check_any_ref_or_leuk_flag, axis=1)
    count_any_ref_flag = df_unclass['Any_ref_arm_flag'].sum()
    pct_any_ref_flag = (count_any_ref_flag / total_unclass) * 100 if total_unclass else 0

    print(f"\n> Triggered ANY reference arm clinical flag or leukopenia:")
    print(f"  {count_any_ref_flag} / {total_unclass} ({pct_any_ref_flag:.2f}%)")


    # D) Any reference arm clinical flag or leukopenia or other cbm clinical flag
    def check_any_flag(row):
        # If it already triggered the manual/OMR composite, it triggers this one
        if row['Any_ref_arm_flag'] == True: return True
        # Leukopenia Auto CBC
        if parse_bool(row.get('Any other clinical flag')) == True: return True
        return False


    df_unclass['Any_flag'] = df_unclass.apply(check_any_flag, axis=1)
    count_any_flag = df_unclass['Any_flag'].sum()
    pct_any_flag = (count_any_flag / total_unclass) * 100 if total_unclass else 0

    print(f"\n> Triggered ANY reference arm clinical flag or leukopenia or other CBM clinical flag:")
    print(f"  {count_any_flag} / {total_unclass} ({pct_any_flag:.2f}%)")

