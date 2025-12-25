import os
import json
import re
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


pbs_wbc_types = [
    "Segmented Neutrophil",
    "Band Neutrophil",
    "Metamyelocyte",
    "Myelocyte",
    "Promyelocyte",
    "Blast",
    "Plasma Cell",
    "Lymphocyte",
    "Large Granular Lymphocyte",
    "Atypical Lymphocyte",
    "Aberrant Lymphocyte",
    "Hairy Cell",
    "Sezary Cell",
    "Monocyte",
    "Basophil",
    "Eosinophil",
    "Unclassified WBC"]

pbs_other_types = [
    "Normoblast",
    "Smudge Cell",
    "platelet",
    "plt_clump",
    "agglutination",
    "rouleaux"]

plt_morphologies = ["giant platelet", "large platelet",  "agranular platelet"]
wbc_morphologies = ['Satellitisms', 'Dohle Bodies', 'Pelger Cell', 'Auer Rods']
rbc_agg_types = ["agglutination", "rouleaux"]

rbc_morphology_types = ['Microcytes', 'Stomatocytes', 'Bite cells', 'Tear drop cells', 'Hypochromia', 'Blister cells',
                        'Ovalocytes', 'Burr Cells', 'Target cells', 'Macrocytes', 'Parasites', 'Elliptocytes',
                        'Spherocytes', 'Helmet cells', 'Pappenheimer', 'Howell-Jolly', 'Schistocytes',
                        'Basophilic stippling', 'Poikilocytosis', 'Poikilocytes', 'Polychromasia',
                        'Sickle cells', 'Acanthocytes', 'RBC Agglutination', 'Rouleaux']

rbc_dists = ['spherocytes', 'ovalocytes', 'microcytes', 'tear_drop', 'stomatocytes', 'bite', 'burr',
             'basophilic_stippling', 'target', 'blister', 'macrocytes', 'hypochromatic', 'howell_jolly', 'spur',
             'schistocytes', 'anisocytosis', 'elliptocytes', 'helmet', 'micro_organisms', 'pappenheimer',
             'poikilocytosis', 'poikilocytes', 'polychromatic', 'sickle']

name_convert = {"platelet_clumps": "plt_clump", "Hairy cell": "Hairy Cell"}

bma_wbc_types = [
    "Blast",
    "Promyelocyte",
    "Myelocyte",
    "Metamyelocyte",
    "Band Neutrophil",
    "Segmented Neutrophil",
    "Plasma Cell",
    "Lymphocyte",
    "Erythroblast",
    "Basophilic Normoblast",
    "Polychromatophilic Normoblast",
    "Normoblast",
    "Monocyte",
    "Eosinophil",
    "Pro eosinophil",
    "Basophil",
    "Mast Cell"]

bma_other_types = [
    "megakaryocyte",
    "stripped",
    "BCD",
    "Negative_in_bma"]

pbs_excluded = ['Unclassified', 'Dirt']
bma_excluded = ['Unclassified']

col_names = {
    'barcode': 'Barcode',
    'totalWBC': 'Total WBC',
    'Segmented Neutrophil': 'Segmented Neutrophil',
    'Band Neutrophil': 'Band Neutrophil',
    'Metamyelocyte': 'Metamyelocyte',
    'Myelocyte': 'Myelocyte',
    'Promyelocyte': 'Promyelocyte',
    'Blast': 'Blast',
    'Plasma Cell': 'Plasma Cell',
    'Lymphocyte': 'Lymphocyte',
    'Large Granular Lymphocyte': 'Large Granular Lymphocyte',
    'Atypical Lymphocyte': 'Atypical Lymphocyte',
    'Aberrant Lymphocyte': 'Aberrant Lymphocyte',
    "Hairy Cell": "Hairy Cell",
    "Sezary Cell": 'Sezary Cell',
    'Monocyte': 'Monocyte',
    'Basophil': 'Basophil',
    'Eosinophil': 'Eosinophil',
    'Normoblast': 'NRBC',
    'Smudge Cell': 'Smudge Cell',
    'Unclassified': 'Unclassified',
    'platelet': 'PLT',
    'plt_clump': 'PLT Clumps',
    'giant platelet': 'Giant Platelet',
    'large platelet': 'Large Platelet',
    'agranular platelet': 'Agranular Platelet',
    'totalRBC': 'Total RBC',
    'spherocytes': 'Spherocytes',
    'ovalocytes': 'Ovalocytes',
    'microcytes': 'Microcytes',
    'tear_drop': 'Tear drop cells',
    'stomatocytes': 'Stomatocytes',
    'bite': 'Bite cells',
    'burr': 'Burr Cells',
    'basophilic_stippling': 'Basophilic stippling',
    'target': 'Target cells',
    'blister': 'Blister cells',
    'macrocytes': 'Macrocytes',
    'hypochromatic': 'Hypochromia',
    'howell_jolly': 'Howell-Jolly',
    'spur': 'Acanthocytes',
    'schistocytes': 'Schistocytes',
    'anisocytosis': 'Anisocytosis',
    'elliptocytes': 'Elliptocytes',
    'helmet': 'Helmet cells',
    'micro_organisms': 'Parasites',
    'pappenheimer': 'Pappenheimer',
    'poikilocytosis': 'Poikilocytosis',
    'poikilocytes': 'Poikilocytes',
    'polychromatic': 'Polychromasia',
    'sickle': 'Sickle cells',
    'agglutination': 'RBC Agglutination',
    'rouleaux': 'Rouleaux',
    'rbc_area_hpf': 'RBC Analysis Area',
    'clump_area_hpf': 'Clump Search Area',
    'mode': 'Scan Mode',
    'ScanID': 'ScanID'}


cbm_columns_order = ['Barcode', 'Total WBC'] + pbs_wbc_types +\
                    ["NRBC", "Smudge Cell"] +\
                    wbc_morphologies +\
                    ['PLT', 'PLT Clumps', 'Giant Platelet', 'Large Platelet', 'Agranular Platelet', 'Total RBC'] +\
                    rbc_morphology_types +\
                    ['RBC Analysis Area', 'Clump Search Area', 'Scan Mode', 'ScanID', 'Creation_Date']

cbm_columns_order = ['ScanID', 'Barcode', 'Total WBC'] + pbs_wbc_types


def threshold_change_parse_wbc(jd, orig_cell, dst_cell, thresh, num_dig=2):
    result = {}
    result['barcode'] = 'None'
    if 'barcode' in jd:
        result['barcode'] = jd['barcode']
    count = {vv: 0 for vv in pbs_wbc_types}
    totalWBC = 0
    for lb in jd["detections"]:
        cell_name = lb['classification']['name']
        if cell_name not in pbs_wbc_types:
            continue
        elif cell_name == orig_cell:
            score = lb['classification']['score']
            # change classification according to newly defined threshold
            if score < thresh:
                cell_name = dst_cell

        # Name conversion - discrepancy from old versions
        if cell_name in name_convert:
            cell_name = name_convert[cell_name]

        count[cell_name] += 1
        totalWBC += 1
    result['totalWBC'] = totalWBC
    for vv in pbs_wbc_types:
        result[vv] = round(count[vv] / (totalWBC + 1e-15) * 100, num_dig)
    return result


def cbm_json_parse(jd, num_dig, result=None):
    if not result:
        result = {}

    if 'barcode' in jd:
        result['barcode'] = jd['barcode']

    count = {vv: 0 for vv in pbs_wbc_types + pbs_other_types + plt_morphologies + wbc_morphologies + rbc_agg_types}
    totalWBC, totalRBC = 0, 0
    for lb in jd["detections"]:
        cell_name = lb['classification']['name']


    # scan info
    result['um_per_pixel'] = jd['um_per_pixel']
    result['mode'] = jd['profile']['scan_profile']




def threshold_change_parse(jd, orig_cell, dst_cell, thresh, num_dig=2,
                           result: dict = {}):
    result['barcode'] = 'None'
    if 'barcode' in jd:
        result['barcode'] = jd['barcode']
    result['um_per_pixel'] = jd['um_per_pixel']

    result['mode'] = jd['profile']['scan_profile']

    count = {vv: 0 for vv in pbs_wbc_types + pbs_other_types + plt_morphologies + wbc_morphologies + rbc_agg_types}
    totalWBC, totalRBC = 0, 0
    for lb in jd["detections"]:
        cell_name = lb['classification']['name']
        score = lb['classification']['score']

        # change classification according to newly defined threshold
        if score is None:
            score = 0
        if cell_name == orig_cell and score < thresh:
            cell_name = dst_cell

            # Name conversion - discrepancy from old versions
            if cell_name in name_convert:
                cell_name = name_convert[cell_name]
            if cell_name in pbs_wbc_types:
                count[cell_name] += 1
                totalWBC += 1
                if lb['morphologies'] is not None:
                    morf = lb['morphologies'].keys()
                    for key in morf:
                        count[key] += 1
            elif cell_name in pbs_other_types:
                count[cell_name] += 1
                if cell_name in ['platelet']:
                    if lb['morphologies'] is not None:
                        morf = lb['morphologies'].keys()
                        for key in morf:
                            count[key] += 1
            elif cell_name in ['rbc_aggregation']:
                if lb['morphologies'] is not None:
                    morf = lb['morphologies'].keys()
                    for key in morf:
                        count[key] += 1
            elif cell_name in ['rbc']:
                totalRBC += 1
            elif cell_name not in pbs_excluded:
                print(f'\033[91mError - {cell_name} is unknown key!!!\033[0m')

    result['totalWBC'] = totalWBC
    for vv in pbs_wbc_types + pbs_other_types[:2]:  # Handle also Smudged and NRBC
        result[vv] = round(count[vv] / (totalWBC + 1e-15) * 100, num_dig)
    for vv in wbc_morphologies[:3]:  # NUT morphologies in percent from NUT
        result[vv] = round(count[vv] / (count['Segmented Neutrophil'] + 1e-15) * 100, num_dig)
    for vv in wbc_morphologies[3:]:  # Blast morphologies in percent from Blast
        result[vv] = round(count[vv] / (count['Blast'] + 1e-15) * 100, num_dig)
    for vv in pbs_other_types[2:]:
        result[vv] = count[vv]
    for vv in plt_morphologies:
        result[vv] = round(count[vv] / (count['platelet'] + 1e-15) * 100, num_dig)
    result['totalRBC'] = totalRBC

    rbc_enabled, rbc_morphologies = False, []
    for rbc_dist in jd['distributions']:
        if rbc_dist['cell_type'] == 'rbc':
            rbc_enabled = True
            for vv in rbc_dist['observations']:
                result[vv['name']] = round(vv['percentage'], num_dig)
                rbc_morphologies.append(vv['name'])
        break
    for region in jd['regions']:
        if region['cell_type'] == 'platelet_clumps':
            pix2hpf = result['um_per_pixel'] ** 2 / (195 ** 2)
            result['clump_area_hpf'] = round(region['bounds']['width'] * region['bounds']['height'] * pix2hpf, 1)
            if 'plt_clump' in count.keys():
                result['plt_clump'] = round(100 * result['plt_clump'] / result['clump_area_hpf'], num_dig)
        if region['cell_type'] == 'rbc':
            pix2hpf = result['um_per_pixel'] ** 2 / (195 ** 2)
            result['rbc_area_hpf'] = round(region['bounds']['width'] * region['bounds']['height'] * pix2hpf, 1)
        if region['cell_type'] == 'rbc_agg':
            pix2hpf = result['um_per_pixel'] ** 2 / (195 ** 2)
            result['rbc_agg_area_hpf'] = round(region['bounds']['width'] * region['bounds']['height'] * pix2hpf, 1)

    if rbc_enabled is False:

        print('no RBC in results')

    else:
        for vv in rbc_agg_types:
            if result.get(vv) == None:
                pass
            else:
                if 'rbc_agg_area_hpf' in result.keys():
                    result[vv] = round(100 * count[vv] / result['rbc_agg_area_hpf'], num_dig)
                elif 'rbc_area_hpf' in result.keys():
                    result[vv] = round(100 * count[vv] / result['rbc_area_hpf'], num_dig)
                    print(f'{vv} detections found but no rbc_agg region')
        # fallback for rbc morphologies not appearing in distributions
        for vv in rbc_dists:
            if vv not in result.keys():
                result[vv] = 0
    return result


def results_df_from_jsons_dir_threshold_change(json_directory, **kwargs):
    PBS_results = []
    for filename in os.listdir(json_directory):
        if filename.endswith('.json'):
            file_path = os.path.join(json_directory, filename)
            jd = json.load(open(file_path))
            result = threshold_change_parse_wbc(jd, **kwargs)  # change back to threshold_change_parse if need more than wbc

            scan_uuid = os.path.splitext(filename)[0]
            result['ScanID'] = scan_uuid

            PBS_results.append(result)
            print(scan_uuid)

    PBS_df = pd.DataFrame(PBS_results)  # Create a Pandas DataFrame
    PBS_df.rename(columns=col_names, inplace=True)

    for col in cbm_columns_order:
        if col not in PBS_df.columns:
            print(f'{col} not in any of jsons')
            PBS_df[col] = ''
    PBS_df = PBS_df.loc[:, cbm_columns_order]

    return PBS_df








