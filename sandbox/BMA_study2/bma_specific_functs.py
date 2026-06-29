# bma_specific_functs.py

import pandas as pd
import os
import sys
trgt_dict = os.path.abspath(r'../cbm_prep')
sys.path.append(trgt_dict)
from sandbox import MetadataBundle, raw_bma_to_df, _ensure_list
from stats_sandbox import get_equivocal_zone_masks, binary_classification_metrics_bootstrap

def removed_for_arbitration(df_raw, df_arb, arbitrator):
    arbitrators = _ensure_list(arbitrator)

    # check which samples went to arbitration
    df = pd.merge(df_raw, df_arb, on=['Site', 'SampleID', 'Method'], how='left', indicator=True)

    # keep reviews which were not sent to arbitration and reviewed by regular reviewer, or ones sent to arbitration and reviewed by arbitrator
    df = df[((df['_merge'] == 'left_only') & ~(df['Investigator'].isin(arbitrators))) | ((df['_merge'] == 'both') & (df['Investigator'].isin(arbitrators)))]

    return df


def generate_fda_equivocal_report(methd_comp, decision_dict, level_a='REF', level_b='TEST', dim_col='Method'):
    """
    Study-specific orchestrator for the BMA FDA MDL request.
    Generates strict and equivocal-adjusted metrics.
    """
    rows = []
    for variable, cutoffs in decision_dict.items():
        if not cutoffs:
            continue

        try:
            # Safely extract matched data using the API
            x, y, ids = methd_comp.get_pairwise_data(level_a, level_b, variable, dim_col)

            # Setup stratification (by Site if available)
            sites = ids.get_level_values('Site').to_numpy() if 'Site' in ids.names else None
            df_strat = pd.DataFrame({'Site': sites}) if sites is not None else pd.DataFrame()
            strat_cols = ['y_true', 'Site'] if sites is not None else ['y_true']

        except Exception as e:
            print(f"Skipping {variable} in FDA MDL Report: {e}")
            continue

        for cutoff in cutoffs:
            if cutoff == 0:
                continue

            # 1. Get binary arrays and the exclusion mask from generic stats tool
            y_true_bin, y_pred_bin, exclude_mask = get_equivocal_zone_masks(x, y, cutoff)
            df_strat['y_true'] = y_true_bin

            # 2. Strict Metrics (All Data)
            strict_metrics = binary_classification_metrics_bootstrap(
                y_true_bin, y_pred_bin, stratify_cols=strat_cols, data_df=df_strat
            )

            # 3. Adjusted Metrics (Equivocal Dropped)
            y_true_adj = y_true_bin[~exclude_mask]
            y_pred_adj = y_pred_bin[~exclude_mask]
            df_adj = df_strat[~exclude_mask].copy()

            # Only calculate if we still have data left after dropping
            if len(y_true_adj) > 0:
                adj_metrics = binary_classification_metrics_bootstrap(
                    y_true_adj, y_pred_adj, stratify_cols=strat_cols, data_df=df_adj
                )
            else:
                adj_metrics = None

            # 4. Format Row
            row = {
                "Variable": variable,
                "Cutoff (%)": cutoff,
                "N (Strict)": len(x),
                "N (Adjusted)": len(y_true_adj),
                "Equivocal Cases Excluded": exclude_mask.sum(),

                "Strict TP": strict_metrics["tp"],
                "Strict FP": strict_metrics["fp"],
                "Strict FN": strict_metrics["fn"],
                "Strict TN": strict_metrics["tn"],
                "Strict Sensitivity (%)": f"{strict_metrics['sensitivity']['value'] * 100:.1f} ({strict_metrics['sensitivity']['ci'][0] * 100:.1f}-{strict_metrics['sensitivity']['ci'][1] * 100:.1f})" if pd.notna(
                    strict_metrics['sensitivity']['value']) else "NA",
                "Strict Specificity (%)": f"{strict_metrics['specificity']['value'] * 100:.1f} ({strict_metrics['specificity']['ci'][0] * 100:.1f}-{strict_metrics['specificity']['ci'][1] * 100:.1f})" if pd.notna(
                    strict_metrics['specificity']['value']) else "NA",
                "Strict OPA (%)": f"{strict_metrics['agreement']['value'] * 100:.1f} ({strict_metrics['agreement']['ci'][0] * 100:.1f}-{strict_metrics['agreement']['ci'][1] * 100:.1f})",
            }

            if adj_metrics:
                row.update({
                    "Adjusted TP": adj_metrics["tp"],
                    "Adjusted FP": adj_metrics["fp"],
                    "Adjusted FN": adj_metrics["fn"],
                    "Adjusted TN": adj_metrics["tn"],
                    "Adjusted Sensitivity (%)": f"{adj_metrics['sensitivity']['value'] * 100:.1f} ({adj_metrics['sensitivity']['ci'][0] * 100:.1f}-{adj_metrics['sensitivity']['ci'][1] * 100:.1f})" if pd.notna(
                        adj_metrics['sensitivity']['value']) else "NA",
                    "Adjusted Specificity (%)": f"{adj_metrics['specificity']['value'] * 100:.1f} ({adj_metrics['specificity']['ci'][0] * 100:.1f}-{adj_metrics['specificity']['ci'][1] * 100:.1f})" if pd.notna(
                        adj_metrics['specificity']['value']) else "NA",
                    "Adjusted OPA (%)": f"{adj_metrics['agreement']['value'] * 100:.1f} ({adj_metrics['agreement']['ci'][0] * 100:.1f}-{adj_metrics['agreement']['ci'][1] * 100:.1f})",
                })

            rows.append(row)

    return pd.DataFrame(rows)


