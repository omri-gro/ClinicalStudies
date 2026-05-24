import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
import statsmodels.formula.api as smf
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import os
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

warnings.simplefilter('ignore', ConvergenceWarning)

# --- CONFIGURATION ---
FILE_REPEATABILITY = r"raw/bma_repeatability.csv"
FILE_REPRODUCIBILITY = r"raw/bma_reproducibility.csv"
OUTPUT_DIR = "results"
PLOT_DIR = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# (Threshold, SD_Limit, CV_Limit, Evaluation_Mode)
AC_CONFIG = {
    # Tier 1: Ultra-Rare
    'Basophil': (2.0, 0.5, 30.0, 'HYBRID'),
    'Mast Cell': (2.0, 0.5, 30.0, 'HYBRID'),
    # Tier 2: Minor / Diagnostic
    'Blast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Promyelocyte': (5.0, 2.0, 30.0, 'HYBRID'),
    'Plasma Cell': (5.0, 2.0, 30.0, 'HYBRID'),
    'Erythroblast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Basophilic Normoblast': (5.0, 2.0, 30.0, 'HYBRID'),
    'Monocyte': (5.0, 2.0, 30.0, 'HYBRID'),
    'Eosinophil': (5.0, 2.0, 30.0, 'HYBRID'),
    # Tier 3: Intermediate
    'Myelocyte': (10.0, 3.0, 25.0, 'HYBRID'),
    'Metamyelocyte': (10.0, 3.0, 25.0, 'HYBRID'),
    'Band Neutrophil': (10.0, 3.0, 25.0, 'HYBRID'),
    'Normoblast': (10.0, 3.0, 25.0, 'HYBRID'),
    'Polychromatophilic Normoblast': (10.0, 3.0, 25.0, 'HYBRID'),
    # Tier 4: Major Populations
    'Lymphocyte': (20.0, 5.0, 25.0, 'HYBRID'),
    'Segmented Neutrophil': (20.0, 5.0, 25.0, 'HYBRID'),
    # Safety fallback (If a name doesn't match perfectly)
    'default': (5.0, 2.0, 30.0, 'HYBRID')
}


def evaluate_status(mean_val, sd, cv, param):
    if param not in AC_CONFIG:
        return 'Pass'
    threshold, sd_limit, cv_limit, mode = AC_CONFIG[param]
    if mode == 'HYBRID':
        if mean_val <= threshold:
            return 'Pass' if sd <= sd_limit else 'Fail'
        else:
            return 'Pass' if cv <= cv_limit else 'Fail'
    return 'Pass'

def calculate_ci(sd, mean, df):
    if df <= 0 or sd == 0: return "-", "-"
    alpha = 0.05
    lower_factor = np.sqrt(df / stats.chi2.ppf(1 - alpha / 2, df))
    upper_factor = np.sqrt(df / stats.chi2.ppf(alpha / 2, df))
    sd_lower, sd_upper = sd * lower_factor, sd * upper_factor
    cv_lower, cv_upper = (sd_lower / mean) * 100, (sd_upper / mean) * 100
    return f"{sd_lower:.2f}-{sd_upper:.2f}", f"{cv_lower:.2f}%-{cv_upper:.2f}%"


def process_repeatability():
    if not os.path.exists(FILE_REPEATABILITY):
        return []
    df_rep = pd.read_csv(FILE_REPEATABILITY)
    param_cols = [c for c in df_rep.columns if c not in ['Sample', 'Day', 'Run', 'Scan']]

    results = []
    for param in param_cols:
        for sample, group in df_rep.groupby('Sample'):
            N = len(group[param].dropna())
            mean_val = group[param].mean()
            if mean_val == 0 or len(group) < 5 or group[param].std() == 0:
                results.append({'Parameter': param, 'Sample': sample, 'N': N, 'Mean': 0.0, 'constant': True})
                continue

            v_scan = group.groupby(['Day', 'Run'])[param].var().mean()
            v_run = max(0, group.groupby(['Day', 'Run'])[param].mean().var() - v_scan / 2)
            v_day = max(0, group.groupby('Day')[param].mean().var() - (v_run / 2 + v_scan / 4))
            v_wl = v_scan + v_run + v_day

            sd_scan, sd_run, sd_day, sd_wl = np.sqrt(v_scan), np.sqrt(v_run), np.sqrt(v_day), np.sqrt(v_wl)
            cv_scan, cv_run, cv_day, cv_wl = (sd_scan / mean_val) * 100, (sd_run / mean_val) * 100, (
                        sd_day / mean_val) * 100, (sd_wl / mean_val) * 100
            df_rep_val = len(group) - len(group[['Day', 'Run']].drop_duplicates())
            df_wl_val = len(group) - 1
            ci_sd_rep, ci_cv_rep = calculate_ci(sd_scan, mean_val, df_rep_val)
            ci_sd_wl, ci_cv_wl = calculate_ci(sd_wl, mean_val, df_wl_val)
            status = evaluate_status(mean_val, sd_wl, cv_wl, param)

            results.append({
                'Parameter': param, 'Sample': sample, 'N': N, 'Mean': mean_val, 'constant': False,
                'Rep_SD': sd_scan, 'Rep_CV': cv_scan, 'BR_SD': sd_run, 'BR_CV': cv_run,
                'BD_SD': sd_day, 'BD_CV': cv_day, 'Total_SD': sd_wl, 'Total_CV': cv_wl,
                'Status': status, 'DF_Rep': df_rep_val, 'Rep_SD_CI': ci_sd_rep, 'Rep_CV_CI': ci_cv_rep,
                'DF_Total': df_wl_val, 'Total_SD_CI': ci_sd_wl, 'Total_CV_CI': ci_cv_wl
            })
    pd.DataFrame(results).to_csv(os.path.join(OUTPUT_DIR, "Repeatability_Results.csv"), index=False)
    return results


def process_reproducibility():
    if not os.path.exists(FILE_REPRODUCIBILITY):
        return []
    df_repro = pd.read_csv(FILE_REPRODUCIBILITY)
    param_cols = [c for c in df_repro.columns if c not in ['Sample', 'Machine', 'Day', 'Scan']]

    results = []
    for param in param_cols:
        for sample, group in df_repro.groupby('Sample'):
            N = len(group[param].dropna())
            mean_val = group[param].mean()
            if mean_val == 0 or len(group) < 5 or group[param].std() == 0:
                results.append({'Parameter': param, 'Sample': sample, 'N': N, 'Mean': 0.0, 'constant': True})
                continue

            v_scan = group.groupby(['Machine', 'Day'])[param].var().mean()
            v_day = max(0, group.groupby(['Machine', 'Day'])[param].mean().var() - v_scan / 4)
            v_machine = max(0, group.groupby('Machine')[param].mean().var() - (v_day / 5 + v_scan / 20))
            v_repro = v_scan + v_day + v_machine

            sd_scan, sd_day, sd_machine, sd_repro = np.sqrt(v_scan), np.sqrt(v_day), np.sqrt(v_machine), np.sqrt(
                v_repro)
            cv_scan, cv_day, cv_machine, cv_repro = (sd_scan / mean_val) * 100, (sd_day / mean_val) * 100, (
                        sd_machine / mean_val) * 100, (sd_repro / mean_val) * 100
            df_rep_val = len(group) - len(group[['Machine', 'Day']].drop_duplicates())
            df_repro_val = len(group) - 1
            ci_sd_rep, ci_cv_rep = calculate_ci(sd_scan, mean_val, df_rep_val)
            ci_sd_repro, ci_cv_repro = calculate_ci(sd_repro, mean_val, df_repro_val)
            status = evaluate_status(mean_val, sd_repro, cv_repro, param)

            results.append({
                'Parameter': param, 'Sample': sample, 'N': N, 'Mean': mean_val, 'constant': False,
                'Rep_SD': sd_scan, 'Rep_CV': cv_scan, 'BD_SD': sd_day, 'BD_CV': cv_day,
                'BS_SD': sd_machine, 'BS_CV': cv_machine, 'Total_SD': sd_repro, 'Total_CV': cv_repro,
                'Status': status, 'DF_Rep': df_rep_val, 'Rep_SD_CI': ci_sd_rep, 'Rep_CV_CI': ci_cv_rep,
                'DF_Total': df_repro_val, 'Total_SD_CI': ci_sd_repro, 'Total_CV_CI': ci_cv_repro
            })
    pd.DataFrame(results).to_csv(os.path.join(OUTPUT_DIR, "Reproducibility_Results.csv"), index=False)
    return results

def plot_precision_profiles(data, study_type="Reproducibility"):
    df = pd.DataFrame(data)
    df = df[df['constant'] == False]
    plot_paths = {}

    for param in df['Parameter'].unique():
        param_data = df[df['Parameter'] == param]
        if param not in AC_CONFIG: continue

        threshold, sd_limit, cv_limit, _ = AC_CONFIG[param]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        below_th = param_data[param_data['Mean'] <= threshold]
        above_th = param_data[param_data['Mean'] > threshold]
        max_x = max(param_data['Mean'].max() * 1.15, threshold * 1.5)
        if np.isnan(max_x) or max_x == 0: max_x = threshold * 2

        # Plot 1: SD vs Mean
        if not below_th.empty: ax1.scatter(below_th['Mean'], below_th['Total_SD'], color='#1f77b4', edgecolor='k', s=55,
                                           zorder=4, label='Active (Evaluated by SD)')
        if not above_th.empty: ax1.scatter(above_th['Mean'], above_th['Total_SD'], color='#a6a6a6', edgecolor='k', s=55,
                                           alpha=0.25, zorder=3, label='Muted (Evaluated by %CV)')
        ax1.plot([0, threshold], [sd_limit, sd_limit], color='#d62728', linestyle='-', linewidth=2.5, zorder=5,
                 label=f'SD Limit (≤ {sd_limit})')
        ax1.fill_between([0, threshold], 0, sd_limit, color='#2ca02c', alpha=0.07, zorder=1)
        ax1.axvline(x=threshold, color='#7f7f7f', linestyle=':', linewidth=1.5, zorder=2)
        ax1.set_title(f'{param} - {study_type} SD', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Mean Concentration (%)');
        ax1.set_ylabel('Total SD')
        ax1.set_xlim(0, max_x)
        ax1.set_ylim(0, max(sd_limit * 1.6, (below_th['Total_SD'].max() if not below_th.empty else 0) * 1.3, 0.5))
        ax1.legend(loc='upper right', fontsize=8.5);
        ax1.grid(True, linestyle='--', alpha=0.4)

        # Plot 2: CV vs Mean
        if not below_th.empty: ax2.scatter(below_th['Mean'], below_th['Total_CV'], color='#a6a6a6', edgecolor='k', s=55,
                                           alpha=0.25, zorder=3, label='Muted (Evaluated by SD)')
        if not above_th.empty: ax2.scatter(above_th['Mean'], above_th['Total_CV'], color='#2ca02c', edgecolor='k', s=55,
                                           zorder=4, label='Active (Evaluated by %CV)')
        ax2.plot([threshold, max_x], [cv_limit, cv_limit], color='#d62728', linestyle='-', linewidth=2.5, zorder=5,
                 label=f'CV Limit (≤ {cv_limit}%)')
        ax2.fill_between([threshold, max_x], 0, cv_limit, color='#2ca02c', alpha=0.07, zorder=1)
        ax2.axvline(x=threshold, color='#7f7f7f', linestyle=':', linewidth=1.5, zorder=2)
        ax2.set_title(f'{param} - {study_type} %CV', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Mean Concentration (%)');
        ax2.set_ylabel('Total %CV')
        ax2.set_xlim(0, max_x)
        ax2.set_ylim(0, max(cv_limit * 2.0, (above_th['Total_CV'].max() if not above_th.empty else 0) * 1.3, 60.0))
        ax2.legend(loc='upper right', fontsize=8.5);
        ax2.grid(True, linestyle='--', alpha=0.4)

        plt.tight_layout()
        img_path = os.path.join(PLOT_DIR, f"{param.replace(' ', '_')}_{study_type}.png")
        plt.savefig(img_path, dpi=150)
        plt.close()
        plot_paths[param] = img_path
    return plot_paths

def set_repeat_table_header(row):
    """ Helper to force MS Word to repeat table headers across pages """
    tr = row._tr
    trPr = tr.get_or_add_trPr()
    tblHeader = OxmlElement('w:tblHeader')
    tblHeader.set(qn('w:val'), "true")
    trPr.append(tblHeader)
    return row

def set_landscape(doc):
    section = doc.add_section()
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = section.page_height, section.page_width
    return section


def build_variance_table(doc, data, study_type):
    if not data: return
    doc.add_heading(f'{study_type} - Variance Components', level=3)
    table = doc.add_table(rows=2, cols=13)
    table.style = 'Table Grid'

    # Repeat first two rows across pages
    set_repeat_table_header(table.rows[0])
    set_repeat_table_header(table.rows[1])

    h1 = table.rows[0].cells
    for i, t in enumerate(['Parameter', 'Sample', 'N', 'Mean']): h1[i].text = t
    h1[4].text = 'Repeatability'
    h1[4].merge(h1[5])
    if study_type == 'Repeatability':
        h1[6].text = 'Between-Run'
        h1[6].merge(h1[7])
        h1[8].text = 'Between-Day'
        h1[8].merge(h1[9])
        h1[10].text = 'Within-Lab (Total)'
        h1[10].merge(h1[11])
    else:
        h1[6].text = 'Between-Day'
        h1[6].merge(h1[7])
        h1[8].text = 'Between-Site'
        h1[8].merge(h1[9])
        h1[10].text = 'Reproducibility (Total)'
        h1[10].merge(h1[11])
    h1[12].text = 'Status'

    h2 = table.rows[1].cells
    for i in range(4, 12, 2): h2[i].text = 'SD'; h2[i + 1].text = '%CV'

    for row in data:
        cells = table.add_row().cells
        cells[0].text, cells[1].text, cells[2].text, cells[3].text = row['Parameter'], str(row['Sample']), str(
            row['N']), f"{row['Mean']:.2f}"
        if row.get('constant'):
            cells[4].text = "Dependent variable is constant."
            cells[4].merge(cells[11])
            cells[12].text = "N/A"
        else:
            cells[4].text = f"{row['Rep_SD']:.2f}"
            cells[5].text = f"{row['Rep_CV']:.1f}%"
            cells[6].text = f"{row.get('BR_SD', row.get('BD_SD')):.2f}"
            cells[7].text = f"{row.get('BR_CV', row.get('BD_CV')):.1f}%"
            cells[8].text = f"{row.get('BD_SD', row.get('BS_SD')):.2f}"
            cells[9].text = f"{row.get('BD_CV', row.get('BS_CV')):.1f}%"
            cells[10].text = f"{row['Total_SD']:.2f}"
            cells[11].text = f"{row['Total_CV']:.1f}%"
            cells[12].text = row['Status']


def build_ci_table(doc, data, study_type):
    if not data: return
    doc.add_heading(f'{study_type} - Degrees of Freedom & 95% Confidence Intervals', level=3)
    table = doc.add_table(rows=2, cols=10)
    table.style = 'Table Grid'

    # Repeat first two rows across pages
    set_repeat_table_header(table.rows[0])
    set_repeat_table_header(table.rows[1])

    h1 = table.rows[0].cells
    for i, t in enumerate(['Parameter', 'Sample', 'N', 'Mean']): h1[i].text = t
    h1[4].text = 'DF'
    h1[5].text = 'Repeatability 95% CI'
    h1[5].merge(h1[6])
    h1[7].text = 'DF'
    h1[8].text = f'{"Within-Lab" if study_type == "Repeatability" else "Reproducibility"} 95% CI'
    h1[8].merge(h1[9])

    h2 = table.rows[1].cells
    h2[5].text = 'SD CI'
    h2[6].text = '%CV CI'
    h2[8].text = 'SD CI'
    h2[9].text = '%CV CI'

    for row in data:
        cells = table.add_row().cells
        cells[0].text, cells[1].text, cells[2].text, cells[3].text = row['Parameter'], str(row['Sample']), str(
            row['N']), f"{row['Mean']:.2f}"
        if row.get('constant'):
            cells[4].text = "-"
            cells[5].text = "-"
            cells[6].text = "-"
            cells[7].text = "-"
            cells[8].text = "-"
            cells[9].text = "-"
        else:
            cells[4].text = str(row['DF_Rep'])
            cells[5].text = row['Rep_SD_CI']
            cells[6].text = row['Rep_CV_CI']
            cells[7].text = str(row['DF_Total'])
            cells[8].text = row['Total_SD_CI']
            cells[9].text = row['Total_CV_CI']


def generate_docx(rep_data, repro_data):
    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Arial'
    style.font.size = Pt(8)

    doc.add_heading('Results', level=1)

    # --- CONSOLIDATED TABLES (LANDSCAPE) ---
    set_landscape(doc)

    if rep_data:
        doc.add_heading('1. Repeatability Data Tables', level=2)
        build_variance_table(doc, rep_data, 'Repeatability')
        doc.add_paragraph("")  # Spacing
        build_ci_table(doc, rep_data, 'Repeatability')
        doc.add_page_break()

    if repro_data:
        doc.add_heading('2. Reproducibility Data Tables', level=2)
        build_variance_table(doc, repro_data, 'Reproducibility')
        doc.add_paragraph("")  # Spacing
        build_ci_table(doc, repro_data, 'Reproducibility')
        doc.add_page_break()

    # --- PRECISION PROFILES (PORTRAIT) ---
    section = doc.add_section()
    section.orientation = WD_ORIENT.PORTRAIT
    section.page_width, section.page_height = section.page_height, section.page_width

    doc.add_heading('3. Precision Profiles', level=2)
    doc.add_paragraph(
        "The plots below display Total Precision (SD and %CV) mapped against the Mean for each sample. Acceptance Criteria limits are denoted by the dashed lines.")

    if rep_data:
        doc.add_heading('3.1 Repeatability Profiles', level=3)
        rep_plots = plot_precision_profiles(rep_data, "Repeatability")
        for param, img_path in rep_plots.items():
            doc.add_paragraph(f"{param}")
            doc.add_picture(img_path, width=Inches(6.0))

    if repro_data:
        doc.add_heading('3.2 Reproducibility Profiles', level=3)
        repro_plots = plot_precision_profiles(repro_data, "Reproducibility")
        for param, img_path in repro_plots.items():
            doc.add_paragraph(f"{param}")
            doc.add_picture(img_path, width=Inches(6.0))

    doc.save(os.path.join(OUTPUT_DIR, "RnR_Consolidated_Report.docx"))


if __name__ == "__main__":
    print("Processing Repeatability...")
    rep_data = process_repeatability()

    print("Processing Reproducibility...")
    repro_data = process_reproducibility()

    print("Generating Word Document with Plots & Consolidated Tables...")
    if rep_data or repro_data:
        generate_docx(rep_data, repro_data)
        print(f"Done! Check the '{OUTPUT_DIR}' folder for results.")
