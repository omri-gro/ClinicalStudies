import pandas as pd
import os
import json


cols_order = ['Case Name', ]

def rprt_snapshot_read(path, morphs, suggestions=True):
    data = json.load(open(path))
    out_dict = data['subreports']['scansSubreports'][0]['subreports']['rbc']['abnormalities']
    for obs in out_dict:
        if out_dict[obs] in ['nd', 'unset']:
            out_dict[obs] = 0
        else:
            out_dict[obs] = int(out_dict[obs])

    out_dict['Case Name'] = data['casePage']['caseData']['case_no']
    out_dict['Case'] = data['casePage']['caseData']['id']
    out_dict['Signed Off At'] = data['signedAt'][:-1]
    out_dict['Reviewer'] = data['signedBy']
    out_dict['UUID'] = data['scan']['scan_no']

    if suggestions:
        for obs in data['rbcObservations']:
            disp_name = obs['display_name']
            if disp_name in morphs:
                dss_grade = obs['suggested']
                out_dict[f'{disp_name} DSS grade'] = 0 if dss_grade == None else int(dss_grade)
                dss_prc = obs['suggestionInfo']
                out_dict[f'{disp_name} DSS percent'] = 0 if dss_prc == '' else float(dss_prc[:-1])

    return out_dict


if __name__ == "__main__":
    suggestions = True
    morphs_dict = {'schistocytes': 'Schistocytes', 'bite_cells': 'Bite cells', 'blister_cells': 'Blister cells', 'basophilic_stippling': 'Basophilic stippling'}
    morphs = list(morphs_dict.values())
    cols_order = ['Case Name'] + morphs + ['Reviewer', 'Case', 'Signed Off At', 'UUID']
    if suggestions:
        morphs_dss_grades = [f'{i} DSS grade' for i in morphs]
        morphs_dss_prc = [f'{i} DSS percent' for i in morphs]
        cols_order = cols_order + morphs_dss_grades + morphs_dss_prc

    # Define the project root dynamically
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../", "../"))
    rprts_path = os.path.join(project_root, 'data', 'imported', 'RBC_extension', 'reports')

    rows_list = []

    case_nums = os.listdir(rprts_path)
    for case in case_nums:
        case_path = os.path.join(rprts_path, case)
        case_conts = os.listdir(case_path)
        if len(case_conts) != 2:
            print(f'Case {case} does not include 2 files')
        else:
            for case_f in case_conts:
                if os.path.splitext(case_f)[-1] == '.json':
                    path = os.path.join(case_path, case_f)
                    rprt_row = rprt_snapshot_read(path, morphs, suggestions)
                    rows_list.append(rprt_row)

    df = pd.DataFrame(rows_list)
    df.rename(columns=morphs_dict, inplace=True)
    df = df[cols_order]
    df.sort_values(by=['Case Name', 'Reviewer'], inplace=True)
    save_path = os.path.join(os.getcwd(), 'reports_aggregate.xlsx')
    df.to_excel(save_path, index=False)




