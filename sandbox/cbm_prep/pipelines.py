from sandbox import *


def clv_pipe(path, site, metadata, method="ClV",
             sheet_name='Sheet1', dir=None, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'],
             only_mean=True, min_inv=0):
    """ Pipeline for preparing ClV results from GoogleForms results/aggregation of PDF results"""
    clv_raw_df = raw_to_df(path, site, method, sheet_name, dir)
    df = stnd_names(clv_raw_df, metadata.alias_map)
    df = diff_from_total(df, metadata, diff_cells="RBC morphology", total_count="TotalRBC")
    df = diff_from_total(df, metadata, diff_cells="PLT morphology", total_count="TotalPLT")
    df = pivot_long(df, id_vars=id_vars)
    df = add_mean_investigator(df, method, min_inv)
    if only_mean:
        df = df.query("Investigator=='Mean Investigator'")
    df = add_grade_column(df, metadata)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)
    return df



""" Old pipeline functions - kept for compatibility"""
def short_pipe(df, metadata, id_vars=["SampleID", "Site", "Method", "FileName"]):
    # standardize column names
    df = stnd_names(df, metadata.alias_map)

    # print warning if WBCs in differential don't add up to ~100
    # need to reach decision if this should be df = check_diff_sum instead when organizing functions
    df = check_diff_sum(df, metadata, tolerance=5)

    df = pivot_long(df, id_vars=id_vars)
    df = add_grade_column(df, metadata)
    df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    print(df)

    return df


def medium_pipe(file_name, site, method, metadata, sheet_name='Sheet1', dir=None,
                id_vars=["SampleID", "Site", "Method", "FileName"], stnrd_id=True):
    if stnrd_id:
        df = raw_to_df(file_name, site, method, sheet_name, dir)
    else:
        df = raw_bma_to_df(file_name, site, method, sheet_name, dir)
    df = short_pipe(df, metadata, id_vars=id_vars)

    return df

def bma_prep_pipeline(file_name, site, method, metadata, sheet_name='Sheet1', dir=None,
                      id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator']):
    df = raw_bma_to_df(file_name, site, method, sheet_name, dir)
    # standardize column names
    df = stnd_names(df, metadata.alias_map)

    if method == 'TEST':  # in future represent this as site rules
        df = calc_diff(df, metadata, diff_cells="NDC", additional_cells="NDC-like")
        df = calc_diff(df, metadata, diff_cells="NDC lineage")
    else:
        # print warning if WBCs in differential don't add up to ~100
        check_diff_sum(df, metadata, tolerance=5, diff_cells="NDC")
    df = pivot_long(df, id_vars=id_vars)
    df = add_grade_column(df, metadata)
    # df = add_pos_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade"], how="all")
    df = df.dropna(subset="Value", how="all")

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = create_derived_variables_long(df, metadata)
    # reconsider performing only after concatenation of all long dataframes
    return df

