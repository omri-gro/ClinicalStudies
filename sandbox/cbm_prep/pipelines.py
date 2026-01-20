from sandbox import *


# temporary, until LongComparisonData split from MethodComparator
def filter_by_df(df,
                 filtering_source: Union[str, Path, pd.DataFrame],
                 filtering_cols: Optional[Sequence[str]] = None,
                 include_rows=False,
                 both=False):
    """
    Using a dataframe of samples to exclude/include (or a file containing one), create new MethodComparator only without/with those rows.
    filtering_cols - Sequence of column names to use for filtering (e.g., ['Site', 'SampleID', 'Investigator', 'Method']),
                     else infer from the filtering df.
    include_rows - Change to True so resulting df will have only samples included in the filtering_source, instead of those that don't  (overwrites both).
    both - Change to True to get both the MethodComparator where specified samples excluded and the one with only those samples included.
         - need to test this functionality in future
    """
    df2 = filtering_source if isinstance(filtering_source, pd.DataFrame) else raw_to_df(filtering_source)
    if df2.empty:
        print(f'\033[34mNothing to filter out!\033[0m')
        return df
    if not filtering_cols:
        filtering_cols = set(df2.columns)
        # if want to use "FileName" as one of filtering columns, define filtering_cols directly and use DataFrame for filtering_source
        filtering_cols.discard("FileName")

    df1 = df.copy()
    common = list(set(df1.columns) & set(filtering_cols))  # columns existing in both

    df1_com = df1[common]
    df2_com = df2[common]

    # for cases where SampleID is int in one df and str in other, try converting df2's SampleID formats to df1's
    # will need re-adjusting once SampleID which are not just numbers with zeros before
    if pd.api.types.is_numeric_dtype(df2_com['SampleID'].dtype) and not pd.api.types.is_numeric_dtype(
            df1_com['SampleID'].dtype):
        df2_com = df2_com[df2_com["SampleID"].notna()]
        num0s = len(df1_com["SampleID"].astype(str).iloc[0])
        ids_in_format = df2_com["SampleID"].astype(int).astype(str).str.zfill(num0s)
        df2_com = df2_com.assign(SampleID=ids_in_format)

    # Build a set-like MultiIndex of unique df2 keys, test membership for df1
    keys_df2 = pd.MultiIndex.from_frame(df2_com.drop_duplicates())
    mask = pd.MultiIndex.from_frame(df1_com).isin(keys_df2)  # samples appearing in both are in True

    if both:
        incl_df = df1.loc[mask].copy()
        exclude_df = df1.loc[~mask].copy()
        return incl_df, exclude_df

    mask = mask if include_rows else ~mask
    out_df = df1.loc[mask].copy()

    # return a new df with filtered rows
    return out_df



def clv_pipe(path, site, metadata, method="ClV",
             sheet_name='Sheet1', dir=None, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'],
             only_mean=True, min_inv=0, mean_inv=True, grade_pos=True, drv_vars=True):
    """ Pipeline for preparing ClV results from GoogleForms results/aggregation of PDF results"""
    clv_raw_df = raw_to_df(path, site, method, sheet_name, dir)
    df = stnd_names(clv_raw_df, metadata.alias_map)
    df = diff_from_total(df, metadata, diff_cells="RBC morphology", total_count="TotalRBC")
    df = diff_from_total(df, metadata, diff_cells="PLT morphology", total_count="TotalPLT")
    df = pivot_long(df, id_vars=id_vars)
    if mean_inv:
        df = add_mean_investigator(df, method, min_inv)
        if only_mean:
            df = df.query("Investigator=='Mean Investigator'")
    if grade_pos:
        df = add_grade_column(df, metadata)
        df = add_pos_column(df, metadata)
        df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    else:
        df = df.dropna(subset=["Value"], how='all')
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    if drv_vars:
        df = create_derived_variables_long(df, metadata)
    return df


def mean_manual_pipe(path, site, metadata, method="manual",
                     sheet_name='Sheet1', dir=None, id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'],
                     only_mean=True, min_inv=0, filtering_source=None):
    """ Pipeline for preparing manual review results from multiple investigators"""
    raw_df = raw_to_df(path, site, method, sheet_name, dir)
    df = stnd_names(raw_df, metadata.alias_map)
    df = calc_diff(df, metadata, additional_cells="WBC-like")
    df = pivot_long(df, id_vars=id_vars)

    if filtering_source is not None:
        df = filter_by_df(df, filtering_source=filtering_source)

    df = add_mean_investigator(df, method, min_inv)
    if only_mean:
        df = df.query("Investigator=='Mean Investigator'")
    df = add_grade_column(df, metadata)
    df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID
    df = create_derived_variables_long(df, metadata)
    df = add_pos_column(df, metadata)
    return df


""" Old pipeline functions - kept for compatibility"""
def short_pipe(df, metadata, id_vars=["SampleID", "Site", "Method", "FileName"]):
    # standardize column names
    df = stnd_names(df, metadata.alias_map)

    # print warning if WBCs in differential don't add up to ~100
    # need to reach decision if this should be df = check_diff_sum instead when organizing functions
    df = check_diff_sum(df, metadata, tolerance=5)

    df = pivot_long(df, id_vars=id_vars)
    df = df.dropna(subset=["SampleID"])  # drop when no readable SampleID

    # prepare graded and boolean values
    grd_params = metadata.variable_groups.get("binary") + metadata.variable_groups.get("grade")
    df = add_grade_column(df, metadata,
                          raw_grade_cond=lambda d: (d["Method"].isin(["OMR", "Manual"]) & d["Variable"].isin(grd_params))
                          )
    df = df.dropna(subset=["Value", "Grade"], how='all')  # drop when neither value or grade in row

    # calculate derived variables (e.g., Variant Lymphocytes)
    df = create_derived_variables_long(df, metadata)
    df = add_pos_column(df, metadata)
    # reconsider performing only after concatenation of all long dataframes

    return df


def medium_pipe(file_name, site, method, metadata, sheet_name='Sheet1', dir=None,
                id_vars=["SampleID", "Site", "Method", "FileName"], stnrd_id=True, **kwargs):
    if stnrd_id:
        df = raw_to_df(file_name, site, method, sheet_name, dir)
    else:
        df = raw_bma_to_df(file_name, site, method, sheet_name, dir)
    df = short_pipe(df, metadata, id_vars=id_vars)

    return df

def bma_prep_pipeline(file_name, site, method, metadata, sheet_name='Sheet1', dir=None,
                      id_vars=["SampleID", "Site", "Method", "FileName", 'Investigator'], **kwargs):
    def stnd_bma_id(name):
        match = re.match(r"^([A-Za-z]*\d+)", str(name))
        return match.group(1) if match else name

    df = raw_bma_to_df(file_name, site, method, sheet_name, dir)
    # standardize column names
    df = stnd_names(df, metadata.alias_map)
    df['SampleID'] = df['SampleID'].apply(stnd_bma_id)

    if method == 'TEST':  # in future represent this as site rules
        df = calc_diff(df, metadata, diff_cells="NDC")
        # df = calc_diff(df, metadata, diff_cells="NDC", additional_cells="NDC-like")
        # df = calc_diff(df, metadata, diff_cells="NDC lineage")
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

