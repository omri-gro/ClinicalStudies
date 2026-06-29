# clinstudtools/preprocessing/graded.py
""" Qualitative Logic Preprocessing Functions """

import pandas as pd
import numpy as np
from clinstudtools.core.metadata import MetadataBundle

# to do: add step/function before/after this function where other shapes of same grade are converted to uniform look (e.g., "x" changed to "Negative")
# to do: integrate this somehow with checking positive cases (will behave like grade with single threshold, or checking where not nan/0 when pregraded)
def add_grade_column(df_long: pd.DataFrame, meta: "MetadataBundle", raw_grade_cond=None):
    """
    Expects df_long columns: SampleID, Variable, Value, Method, Site (case-sensitive).
    Adds column for grades based on number in value column.
    Uses MetaDataBundle, which is expected to include grading_specs dict with thresholds and grades for variable names.
    raw_grade_cond is callable condition defining when numeric grades were provided in df_long instead of values
    """
    # currently data is converted to grades but not to pd.CategoricalDtype - if want to change add the _coerce_grade_categorical function

    df = df_long.copy()
    req = {"Variable", "Value", "Method", "Site"}
    if not req.issubset(df.columns):
        raise ValueError(f"Expected columns: {req}")

    df["Grade"] = pd.NA
    df["Grade_from"] = pd.NA

    # Copy-as-is for non-numeric rows or those filling given_grade_cond
    value_num = pd.to_numeric(df["Value"], errors="coerce")
    non_numeric = value_num.isna() & df["Value"].notna()  # strings like "N/A", "Not evaluable", etc.

    # for values provided as numeric grades and previously read as values
    if raw_grade_cond is not None:
        non_numeric = non_numeric | raw_grade_cond(df)
        if not isinstance(non_numeric, pd.Series) or non_numeric.dtype != bool:
            raise TypeError("raw_grade_cond must return a boolean pandas Series")

    if non_numeric.any():
        df.loc[non_numeric, "Grade"] = df.loc[non_numeric, "Value"]
        df.loc[non_numeric, "Grade_from"] = "provided"
        # in these cases, the "values" were really just grades, so no need to include them in continuous values analysis
        df.loc[non_numeric, "Value"] = pd.NA


    # mask indices where MetaDataBundle claims raw was given as grade
    # currently not in use at pregraded_index not defined for either studies' context configuration
    if meta.pregraded_index is not None:
        keys = pd.MultiIndex.from_frame(df[["Site", "Method", "Variable"]])
        pregraded = meta.pregraded_index.reindex(keys, fill_value=False).to_numpy()
    else:
        pregraded = np.zeros(len(df), dtype=bool)

    # copy when source says raw values are already grades  (notice that currently nan of some types are also copied if belonging to these variables)
    mask_pregraded_numeric = pregraded & ~non_numeric
    if mask_pregraded_numeric.any():
        df.loc[mask_pregraded_numeric, "Grade"] = df.loc[mask_pregraded_numeric, "Value"]
        df.loc[mask_pregraded_numeric, "Grade_from"] = "provided"
        # in these cases, the "values" were really just grades, so no need to include them in continuous values analysis
        df.loc[mask_pregraded_numeric, "Value"] = pd.NA


    # derive grade from gradable numeric values (which were not already provided as grades)
    gradable_vars = set(meta.grading_specs.keys())
    need_convert = (~pregraded) & (~non_numeric) & df["Variable"].isin(gradable_vars)
    if need_convert.any():
        for var in sorted(gradable_vars):
            mv = need_convert & df["Variable"].eq(var)
            if not mv.any():
                continue
            spec = meta.grading_specs[var]
            vals = pd.to_numeric(df.loc[mv, "Value"], errors="coerce")  # -> float dtype w/ np.nan
            try:
                df.loc[mv, "Grade"] = cut_series_to_categorical(
                    vals,
                    thresholds=spec["thresholds"],
                    grades=spec["grades"],
                    right_closed=spec.get("right_closed", True),
                    clamp_out_of_range=spec.get("clamp_out_of_range", True))
                df.loc[mv, "Grade_from"] = "derived"
            except ValueError or ValueError as e:
                print(f"\033[93mGrading specs error for {var}: {e}. Values not converted to grades.\033[0m")
                df.loc[mv, "Grade_from"] = "conversion error"
    return df


def add_pos_column(df_long: pd.DataFrame, meta: "MetadataBundle",
                   normal_grades=[0, "0", "Normal", "Negative", "normal", "negative", "Not present", "not present"]):
    """
    Adds boolean column for positivity based on normal ranges in MetaDataBundle.
    If grade already exists, treat values in normal_vals as False (negative) and rest as True,
    Else use normal ranges and the "Value" column if possible.
    """
    # create empty boolean column
    df = df_long.copy()
    df["Positive"] = np.nan
    df["Positive"] = df["Positive"].astype('boolean')

    # where grade exists, use it for positivity  -  this section might need changing if grade ever not related to positivity
    grade_negative = df["Grade"].isin(normal_grades)
    grade_positive = df["Grade"].notna() & ~df["Grade"].isin(normal_grades)
    df.loc[grade_negative, "Positive"] = False
    df.loc[grade_positive, "Positive"] = True

    # find values where positivity will be based on normal ranges
    value_num = pd.to_numeric(df["Value"], errors="coerce")
    numeric = value_num.notna()
    normal_ranges = getattr(meta, "normal_ranges", {})
    norm_range_vars = set(normal_ranges.keys())
    need_convert = df["Positive"].isna() & numeric & df["Variable"].isin(norm_range_vars)

    # convert based on normal ranges
    if need_convert.any():
        for var in sorted(norm_range_vars):
            mv = need_convert & df["Variable"].eq(var)
            if not mv.any():
                continue
            norm_range = normal_ranges[var]
            if len(norm_range) == 2:
                df.loc[mv, "Positive"] = ~value_num[mv].between(norm_range[0], norm_range[1], inclusive="both")
            else:
                print(f'{norm_range} is not an appropriate normal range for {var}')
    return df


def cut_series_to_categorical(x: pd.Series,
                              thresholds,
                              grades,
                              *,
                              right_closed: bool = True,
                              clamp_out_of_range: bool = False) -> pd.Categorical:
    """
    Bin continuous values into grades.

    Args:
        x: numeric Series.
        thresholds: list-like of bin edges (e.g., [0,5,10,20,101]).
        grades: list-like of labels (e.g., [0,1,2,3] or ["none","mild","mod","sev"]).
        right_closed: if True, bins are right-closed ( (a,b] ).
        clamp_out_of_range: clip x to [min(thresholds), max(thresholds)] before cutting. Only used when len(thresholds) == len(grades) + 1.

    - If len(thresholds) == len(grades) - 1: treated as *interior cut points*; we pad with -inf/+inf
    - If len(thresholds) == len(grades) + 1: treated as full bin edges; must already bracket all bins.

    Returns:
        pandas Series with categories exactly as in `grades`.
    """
    y = pd.to_numeric(x, errors="coerce")  # robust to pd.NA / strings / etc.
    if len(thresholds) == len(grades) - 1:
        edges = np.concatenate(([-np.inf], thresholds, [np.inf]))
    elif len(thresholds) == len(grades) + 1:
        edges = thresholds
        if clamp_out_of_range:
            y = y.clip(lower=edges[0], upper=edges[-1])
    else:
        raise ValueError(f"\033[93mGrades must have one more or one less member than thresholds\033[0m")

    # Validate monotonicity (required by pandas.cut)
    if not np.all(np.diff(edges) > 0):
        raise ValueError("thresholds/edges must be strictly increasing.")

    out = pd.cut(y, bins=edges, labels=grades, right=right_closed, include_lowest=True)
    return out

