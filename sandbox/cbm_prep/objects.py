""" consider splitting into multiple files later """

import pandas as pd
from typing import List, Tuple, Optional, Dict, Any, Union, Sequence
from pathlib import Path
from itertools import product
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import warnings
from dataclasses import dataclass, asdict
from matplotlib.backends.backend_pdf import PdfPages

import sandbox as sb
from pipelines import medium_pipe, bma_prep_pipeline, mean_manual_pipe
import regressions as reg  # placeholder until integrating regression functions
from stats_sandbox import binary_classification_metrics, binary_classification_metrics_bootstrap
import plotting
import sys

sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies')
from clinstudtools.transforms import filter_by_reference, filter_by_condition
from clinstudtools.utils import ensure_list
from clinstudtools.table_integrity import safe_pivot, filter_to_ids
from clinstudtools.comparison_labels import normalize_filter_mapping, format_filter_label

@dataclass
class RegressionResult:
    slope: float
    intercept: float
    slope_ci: Optional[Tuple[float, float]] = None
    intercept_ci: Optional[Tuple[float, float]] = None
    r: Optional[float] = None  # correlation coefficient
    r2: Optional[float] = None    # coefficient of determination
    iterations: Optional[int] = None
    ci_method: str = "Bootstrap"
    reg_method: str = "Deming"
    r_type: str = "Pearson Correlation"
    n: Optional[int] = None  # number of datapoints

    def __post_init__(self):
        # auto-compute missing r2 from r, or r from r2
        if self.r is None and self.r2 is not None:
            self.r = (self.r2 ** 0.5) if self.r2 >= 0 else None
        elif self.r2 is None and self.r is not None:
            self.r2 = self.r ** 2

    def __repr__(self) -> str:  # to do - create function or attribute for string representations with CIs
        """Compact string representation for quick debugging/printing."""
        return (
            f"<RegressionResult slope={self.slope:.3f}, intercept={self.intercept:.3f}, "
            f"r={self.r if self.r is not None else 'NA'}, "
            f"method={self.reg_method}, ci_method={self.ci_method}>"
        )

    def as_dict(self) -> dict[str, Any]:
        """Convert result to a dict (safe for missing CI)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RegressionResult":
        """
        Factory method to normalize old regression result dicts into a RegressionResult.
        Handles naming differences, tuples, missing values.
        """

        # Map old-style keys to new attributes
        slope_ci = (
            data.get("slope_ci")
            or (data.get("slope_ci_bottom"), data.get("slope_ci_top"))
            if "slope_ci_bottom" in data and "slope_ci_top" in data
            else None
        )
        intercept_ci = (
            data.get("intercept_ci")
            or (data.get("intercept_ci_bottom"), data.get("intercept_ci_top"))
            if "intercept_ci_bottom" in data and "intercept_ci_top" in data
            else None
        )

        return cls(
            slope=data.get("slope"),
            intercept=data.get("intercept"),
            slope_ci=slope_ci,
            intercept_ci=intercept_ci,
            r=data.get("r") or data.get("correlation_coefficient"),
            r2=data.get("r2") or data.get("R2"),
            iterations=data.get("iterations"),
            ci_method=data.get("ci_method") or data.get("CI method", "Bootstrap"),
            reg_method=data.get("reg_method") or data.get("regression method", "Deming"),
            n=data.get("N") or data.get("n")
        )


class MethodComparator:
    def __init__(self, df: pd.DataFrame, measurement_col='Value'):
        """
        Initialize comparator with long-format dataframe.
        Expected columns: SampleID, Site, Method, Variable, Value
        """
        self.df = df.copy()
        self.measurement_col = measurement_col  # currently not used
        # stores regression results
        self.results = []
        # stores qualitative and semi-quantitative metrics (sensitivity, kappa, etc.)
        self.metrics = []

    @classmethod
    def from_paths_dict(cls, paths: dict, metadata: sb.MetadataBundle, measurement_col='Value',
                        more_id_vars=None, bma=False, inv=False, **kwargs):
        # possible kwargs - dir, stnrd_id, stnrd_id, min_inv, only_mean, filtering_source etc.
        # maybe don't use, might be better to create something like that every time due to all edge cases

        # function that gets dict with paths as values and descriptors of method/site/etc as keys
        # will read, prepare, concatenate and convert to MethodComparator
        # assumes all files are in same directory

        # to do: clean and add error messages
        # to do: use metadata to pass info like stnrd_id?
        # to do: make it possible for the values in paths dict to be a list of multiple paths, instead of a single path

        # list of columns that will not be pivoted
        id_vars = ["SampleID", "Site", "Method", "FileName", 'Investigator'] if (bma or inv) else ["SampleID", "Site", "Method", "FileName"]
        if isinstance(more_id_vars, list):
            id_vars += more_id_vars

        df_srcs_list = []
        for key, path in paths.items():
            try:
                if isinstance(key, tuple) or isinstance(key, list):
                    site, method, *add_inf = key
                    if add_inf:
                        reviewer = add_inf[0]
                    else:
                        reviewer = None
                else:
                    site = method = reviewer = None
                if bma:
                    df = bma_prep_pipeline(path, site, method, metadata, id_vars=id_vars, **kwargs)
                elif inv:
                    df = mean_manual_pipe(path, site, metadata=metadata, method=method, id_vars=id_vars, **kwargs)
                else:
                    df = medium_pipe(path, site, method, metadata, id_vars=id_vars, **kwargs)
                df_srcs_list.append(df)
            except Exception as e:
                print(f'\033[91mError when importing {key}: {e}\033[0m')
        all_dfs = pd.concat(df_srcs_list)
        return MethodComparator(all_dfs, measurement_col)


    @classmethod
    def from_excel_sheets(cls, excel_path):
        # to do: goes through all sheets in an excel, read each into a df, prepare, concatenate and convert to MethodComparator
        # need to decide on how to specify site/method/investigator
        pass

    """ Methods creating new MethodComparators from existing ones"""
    def filter_by_df(self, filtering_source, **kwargs):
        warnings.warn(
            "MethodComparator.filter_by_df is deprecated. Use clinstudtools.filtering.filter_by_reference instead on the raw DataFrame.",
            DeprecationWarning, stacklevel=2
        )
        # Update self.df in place and return self for chaining
        self.df = filter_by_reference(self.df, filtering_source, **kwargs)
        return self


    def only_when_cond(self, condition: str, filtering_cols=None, non_filtered_vars=None):
        warnings.warn("MethodComparator.only_when_cond is deprecated. Use clinstudtools.transforms.filter_by_condition.", DeprecationWarning, stacklevel=2)
        needed_cases = filter_by_condition(self.df, condition)
        out_df = filter_by_reference(self.df, needed_cases, filtering_cols, include_rows=True)
        return MethodComparator(out_df, self.measurement_col)


    def filter_by_unclass(self, unclass='Unclassified WBC', threshold=10, filtering_cols=None):
        """
        Create a new MethodComparator with only cases for which unclass < threshold.
        """
        warnings.warn(
            "MethodComparator.filter_by_unclass is deprecated. Perform this filtering "
            "on the raw DataFrame before instantiating MethodComparator.",
            DeprecationWarning, stacklevel=2
        )
        cond = f"Variable == '{unclass}' and Value < {threshold}"
        needed_cases = filter_by_condition(self.df, cond)
        out_df = filter_by_reference(self.df, needed_cases, filtering_cols, include_rows=True)
        return MethodComparator(out_df, self.measurement_col)


    def apply_to_df(self, function: str, *args, inplace=False, **kwargs):
        # still needs testing, use 'Pythonic method delegation' for improvements
        """
        Apply a DataFrame method to self.df by name.

        Parameters
        ----------
        function : str
            Name of the DataFrame method to call (e.g., 'query', 'dropna', 'assign').
        *args, **kwargs :
            Arguments to pass to the DataFrame method.
        inplace : bool, optional (default=False)
            If True, modify the internal DataFrame and return self.
            If False, return a new MethodComparator with the modified DataFrame.

        Returns
        -------
        MethodComparator or other
            If the called method returns a DataFrame and inplace=False,
            returns a new MethodComparator. Otherwise returns the raw result.
        """
        # If the method relies on caller environment variables, inject the caller's frame.
        if function in ('query', 'eval'):
            # Grab the caller's frame (1 level up the stack)
            caller_frame = sys._getframe(1)

            # Pass the caller's locals and globals to pandas
            kwargs.setdefault('local_dict', caller_frame.f_locals)
            kwargs.setdefault('global_dict', caller_frame.f_globals)

        # Get the DataFrame method dynamically
        method = getattr(self.df, function, None)
        if method is None or not callable(method):
            raise AttributeError(f"DataFrame has no method '{function}'")

        # Apply the method
        result = method(*args, **kwargs)

        # Handle return type
        if isinstance(result, pd.DataFrame):
            if inplace:
                self.df = result
                return self
            else:
                return MethodComparator(result)
        else:
            # For Series, scalars, or other objects — just return as-is
            return result


    def export_comparison_matrix(self, out_path=None, **kwargs) -> pd.DataFrame:
        # user can include needed_vals both and needed_grades (both need to be Iterable[str]),
        # to ask for values for some variables and grades for others
        if 'needed_vals' in kwargs and 'needed_grades' in kwargs:
            # create kwargs for needed values
            val_kwargs = kwargs
            val_kwargs.update({'needed_vars': kwargs["needed_vals"], 'value_col': "Value"})
            val_df = sb.to_comparison_matrix(self, metadata=getattr(self, "metadata", None), **val_kwargs)

            # create kwargs for needed values
            grade_kwargs = kwargs
            grade_kwargs.update({'needed_vars': kwargs["needed_grades"], 'value_col': "Grade"})
            grade_df = sb.to_comparison_matrix(self, metadata=getattr(self, "metadata", None), **grade_kwargs)

            row_identifiers = kwargs.get('row_identifiers', ("SampleID", "Site"))
            wide_df = pd.merge(val_df, grade_df, on=row_identifiers, how='outer')
        elif 'needed_grades' in kwargs:
            # create kwargs for needed values
            grade_kwargs = kwargs
            grade_kwargs.update({'needed_vars': kwargs["needed_grades"], 'value_col': "Grade"})
            wide_df = sb.to_comparison_matrix(self, metadata=getattr(self, "metadata", None), **grade_kwargs)
        else:
            wide_df = sb.to_comparison_matrix(self, metadata=getattr(self, "metadata", None), **kwargs)

        if out_path:
            sb.write_df_to_file(wide_df, out_path)
        return wide_df

    def clean_calculations(self):
        self.results = {}
        self.metrics = {}

    def _col_from_ids(
            self,
            ids: Union[pd.Index, pd.MultiIndex],
            *,
            variable: str,
            dim_col: str,
            level: str,
            column: str,
            id_cols: Sequence[str] = ("SampleID", "Site"),
            row_filters: Optional[dict] = None,
            on_duplicates: str = "raise",
            sort_by: Optional[Sequence[str]] = None,
    ) -> np.ndarray:
        """
        Return an array of `column` values for a single `level` of `dim_col`,
        aligned to the order of `ids` produced by _prepare_pairwise_arrays.

        Parameters
        ----------
        ids:
            Index/MultiIndex from _prepare_pairwise_arrays (the wide index).
            Must correspond to id_cols (e.g. SampleID+Site).
        variable:
            Variable name (same as used in _prepare_pairwise_arrays).
        dim_col:
            Comparison dimension column (e.g. "Method", "Investigator").
        level:
            Level within dim_col to extract for (e.g. "OMR" or "R1").
        column:
            Column to retrieve values from (e.g. "Positive").
        id_cols:
            The columns used to build ids (must match how ids were created).
        row_filters:
            Optional restrictions (e.g. {"Site": ["A","B"], "Method": "manual"}).
            Note: if row_filters includes dim_col, it should be compatible with `level`.
        on_duplicates, sort_by:
            Same semantics as in safe_pivot/robust_dup.

        Returns
        -------
        np.ndarray aligned to `ids` order.
        """

        df = self.df

        # Filter to the variable and other row filters (but do NOT require the comparison pair)
        subset = df[df["Variable"] == variable].copy()
        subset = self._apply_row_filters(subset, row_filters)

        # Only the one level we want
        subset = subset[subset[dim_col] == level]

        # Restrict to the ids we care about
        # Build a temporary index matching id_cols and filter by ids
        subset = filter_to_ids(subset, ids, id_cols)

        # Pivot to align the requested column (no real "columns" dimension needed since only one level,
        # but we can keep it consistent by pivoting on dim_col)
        wide = safe_pivot(
            subset,
            index=id_cols,
            columns=dim_col,
            values=column,
            on_duplicates=on_duplicates,
            sort_by=sort_by,
        )

        # Reindex to match the order of ids from the original comparison
        wide = wide.reindex(ids)

        # If the column is missing entirely (no rows for that level), return all NaN
        if level not in wide.columns:
            return np.full(shape=(len(ids),), fill_value=np.nan)

        return wide[level].to_numpy()


    def _apply_row_filters(
            self,
            df: pd.DataFrame,
            row_filters: Optional[dict] = None,
    ) -> pd.DataFrame:
        """
            Apply column-based row filters.

            row_filters example:
                {
                    "Site": ["A", "B"],
                    "Method": "OMR",
                    "Investigator": ["R1", "R2"],  # must use names of column in df as keys
                    "MCV": lambda s:s > 80         # can use callables, list-like or a single str/bool/float as values
                }
        """
        if not row_filters:
            return df

        for col, allowed in row_filters.items():
            if col not in df.columns:
                raise ValueError(f'Dataframe does not have {col} as column.')
            if allowed is None:
                continue
            if callable(allowed):
                df = df[allowed(df[col])]
            else:
                allowed = ensure_list(allowed)
                df = df[df[col].isin(allowed)]
        return df


    def _prepare_pairwise_arrays(
            self,
            level_a: str,
            level_b: str,
            variable: str,
            dim_col: str,
            measurement_col: str = "Value",
            row_filters: Optional[dict] = None,
            id_cols: Tuple[str, ...] = ("SampleID", "Site"),
            on_duplicates: str = "raise",
            sort_by: Optional[Sequence[str]] = None,
    ) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Internal helper: return matched (x, y, ids).

        Parameters:
            level_a: name of reference review method/first reviewer (e.g., 'OMR', 'Rev1', etc.)
            level_b: name of test review method/second reviewer (e.g., 'CBM', 'Rev2', etc.)
            variable: single variable to filter on (e.g., 'Total Neutrophil')
            dim_col: basis of comparison (usually either "Method" or "Investigator")
            measurement_col: column name holding the measurement values
            row_filter: dictionary to be sent to _apply_row_filters, with keys being names of columns in self.df
            id_cols: columns that will be used for identifying x/y pairs
            on_duplicates: behaviour when duplicates exist in self.df on column names in id_cols
            sort_by: column names to be used when on_duplicates=='first' or 'last'

        Returns:
            x: array of level_a measurements
            y: array of level_b measurements
            ids: list of tuples of id_cols values corresponding to each x/y pair
        """
        subset = self.df[self.df["Variable"] == variable].copy()
        subset = self._apply_row_filters(subset, row_filters)
        subset = subset[subset[dim_col].isin([level_a, level_b])]
        if measurement_col == "Value":   # assumes that functions on the "Value" column are always numeric
            subset[measurement_col] = pd.to_numeric(subset[measurement_col], errors="coerce")

        subset.dropna(subset=[measurement_col])

        wide = safe_pivot(subset, index=id_cols, columns=dim_col, values=measurement_col,
                          on_duplicates=on_duplicates, sort_by=sort_by,)


        # ensure both levels are present
        missing_cols = {level_a, level_b} - set(wide.columns)
        if missing_cols:
            missing_cols_str = ", ".join(f"{item}" for item in missing_cols)
            raise ValueError(f"Missing {variable} variable for {missing_cols_str}")

        wide = wide.dropna(subset=[level_a, level_b])
        if len(wide) == 0:
            raise ValueError(f"Missing {variable} variable in either {level_a} or {level_b}")

        x = wide[level_a].to_numpy()
        y = wide[level_b].to_numpy()
        ids = wide.index

        return x, y, ids


    def _prepare_arrays(
        self,
        ref_method: str,
        test_method: str,
        variable: str,
        measurement_col: Optional[str] = 'Value',
        site_filter: Optional[List[str]] = None,
        on_duplicates: str = "raise",  # "raise" (default) or "first"
    ) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Internal helper: return matched (x, y, ids).

        Parameters:
        ref_method: name of reference review method (e.g., 'OMR')
        test_method: name of test review method (e.g., 'DSS')
        variable: variable to filter on (e.g., 'Total Neutrophil')
        site_filter: list or set of sites to include (optional)
        measurement_col: column name holding the measurement values
        on_duplicates: either raise error if duplicates of same observation are found, or just take first one

    Returns:
        x: array of reference method measurements
        y: array of test method measurements
        ids: list of (Site, SampleID) tuples corresponding to each x/y pair
        """
        # see _prepare_arrays_strict if ever need to improve
        if not (isinstance(ref_method, str) and isinstance(test_method, str) and isinstance(variable, str) and isinstance(measurement_col, str)):
            raise TypeError('ref_method, test_method, variable and measurement_col should all be strings')

        subset = self.df[self.df["Variable"] == variable].copy()
        if site_filter is not None and site_filter is not [None]:
            if isinstance(site_filter, str):
                site_filter = [site_filter]
            subset = subset[subset["Site"].isin(site_filter)]
        subset = subset[subset["Method"].isin([ref_method, test_method])]

        # Ensure uniqueness per (SampleID, Site, Method) first (raise if not) - consider moving to sandbox or own method
        dup_key = ["SampleID", "Site", "Method"]
        dup_mask = subset.duplicated(dup_key, keep=False)
        if dup_mask.any():
            dup_summary = (
                subset.loc[dup_mask, dup_key]
                .value_counts()
                .reset_index(name="count")
                .sort_values("count", ascending=False)
            )

            if on_duplicates == "raise":
                preview = dup_summary.head(10).to_string(index=False)
                raise ValueError(
                    "Duplicate rows per (SampleID, Site, Method) before pivot.\n"
                    f"Examples (up to 10):\n{preview}"
                )

            elif on_duplicates == "first":
                # Keep the first occurrence per key, drop the rest
                # Using sort_index() first to make 'first' deterministic wrt input order
                before = len(subset)
                subset = subset.sort_index().drop_duplicates(dup_key, keep="first")
                after = len(subset)
                kept = len(dup_summary)
                dropped = before - after
                print(
                    f"[dedupe] on_duplicates='first': kept first row for {kept} "
                    f"duplicate key(s); dropped {dropped} rows."
                )

            else:
                raise ValueError("on_duplicates must be 'raise' or 'first'.")

        # Pivot: SampleID+Site as index, Methods as columns
        pivoted = subset.pivot_table(
            index=["SampleID", "Site"],
            columns="Method",
            values=measurement_col,
            aggfunc="first",
            dropna=False
        )

        # Drop missing pairs (cases where at least one of the methods has NaN)
        pivoted = pivoted.dropna(subset=[ref_method, test_method])

        x = pivoted[ref_method].values
        y = pivoted[test_method].values
        ids = pivoted.index.values  # add .get_level_values("SampleID") before .value to take only SampleIDs
        # To do: add an option for ids to include the site as well

        return x, y, ids

    def get_pairwise_data(self, level_a: str, level_b: str, variable: str, dim_col: str = "Method",
                          on_duplicates: str = "first"):
        """
        Public method to safely extract aligned arrays (x, y, ids) for custom external analysis.
        """
        return self._prepare_pairwise_arrays(
            level_a, level_b, variable, dim_col, on_duplicates=on_duplicates
        )

    def confusion_matrix(self,
                         level_a: str,
                         level_b: str,
                         variable: str,
                         row_filters: Optional[List[str]] = None,
                         site_filter: Optional[List[str]] = None,
                         dim_col="Method",
                         measurement_col="Grade",  # Default to Grade to avoid numeric coercion
                         id_cols: Tuple[str, ...] = ("SampleID", "Site"),
                         **kwargs):
        """Generate a confusion matrix (crosstab) for qualitative/ordinal comparisons."""
        if row_filters is None and site_filter is not None:
            row_filters = {"Site": site_filter}

        # Extract pairwise categorical/ordinal arrays
        x, y, ids = self._prepare_pairwise_arrays(level_a, level_b, variable, dim_col,
                                                  row_filters=row_filters,
                                                  measurement_col=measurement_col,
                                                  id_cols=id_cols,
                                                  on_duplicates="first")

        # Generate the confusion matrix
        # Note: Using categorical dtypes in your DataFrame earlier in the pipeline
        # ensures crosstab includes empty classes (0 counts).
        xtab = pd.crosstab(x, y,
                           rownames=[f"Ref ({level_a})"],
                           colnames=[f"Test ({level_b})"],
                           dropna=False)

        record = {
            "comparison_dim": dim_col,
            "level_a": level_a,
            "level_b": level_b,
            "variable": variable,
            "stratum": dict(row_filters or {}),
            "x": x,
            "y": y,
            "ids": ids,
            "xtab": xtab  # Store the raw dataframe matrix
        }

        # Store in metrics (or create a dedicated self.qualitative list)
        self.metrics.append(record)
        return record


    def sen_spe(self,
                level_a: str,  # e.g. reference method
                level_b: str,  # e.g. test method
                variable: str,
                row_filters: Optional[List[str]] = None,
                site_filter: Optional[List[str]] = None,  # for backwards compatibility
                dim_col="Method",
                measurement_col="Positive",
                id_cols: Tuple[str, ...] = ("SampleID", "Site"),  # identifiers of datapoints
                cis=True,
                **kwargs):
        """
        Calculate sensitivity/specificity for one variable/site/method pair.
        Possible kwargs: m_boot, alpha, random_state.
        # consider merging this function and fit function
        """
        if row_filters is None and site_filter is not None:
            row_filters = {"Site": site_filter}

        x, y, ids = self._prepare_pairwise_arrays(level_a, level_b, variable, dim_col, row_filters=row_filters,
                                                  measurement_col=measurement_col, id_cols=id_cols,
                                                  on_duplicates="first")

        if cis:
            result = binary_classification_metrics_bootstrap(x, y, **kwargs)
        else:
            result = binary_classification_metrics(x, y)

        record = {
            # identity / metadata
            "comparison_dim": dim_col,
            "level_a": level_a,
            "level_b": level_b,
            "variable": variable,
            "stratum": dict(row_filters or {}),  # includes split_by values

            # data
            "x": x,
            "y": y,
            "ids": ids,

            # result payload
            "mtr": result,
        }

        # save results
        self.metrics.append(record)
        return record

    def fit(self,
            level_a: str,  # e.g. reference method
            level_b: str,  # e.g. test method
            variable: str,
            row_filters: Optional[List[str]] = None,
            site_filter: Optional[List[str]] = None,  # for backwards compatibility
            dim_col="Method",
            measurement_col="Value",
            id_cols: Tuple[str, ...] = ("SampleID", "Site"),  # identifiers of datapoint
            **kwargs  # e.g., reg_method, ci, lambda_, etc.
            ):
        """Run regression for one variable/site/method pair."""
        if row_filters is None and site_filter is not None:
            row_filters = {"Site": site_filter}

        x, y, ids = self._prepare_pairwise_arrays(level_a, level_b, variable, dim_col, row_filters=row_filters,
                                                  measurement_col=measurement_col, id_cols=id_cols,
                                                  on_duplicates="first")

        # placeholder regression_func → replace with your real implementation
        stats = reg.regression_comp(x, y, **kwargs)  # placeholder - replace later with different regression functions
        if type(stats) == RegressionResult:
            result = stats
            if result.n is None:
                result.n = len(x)
        elif type(stats) == dict:  # backwards compatibility - for old regression functions returning dicts
            result = RegressionResult.from_dict(stats)
        else:
            result = RegressionResult(**stats)

        record = {
            # identity / metadata
            "comparison_dim": dim_col,
            "level_a": level_a,
            "level_b": level_b,
            "variable": variable,
            "stratum": dict(row_filters or {}),  # includes split_by values

            # data
            "x": x,
            "y": y,
            "ids": ids,

            # result payload
            "reg": result,
        }

        # add count of positives if possible - need to make this its own method later
        if 'Positive' in self.df.columns and pd.api.types.is_bool_dtype(self.df["Positive"]):
            try:
                ref_pos = self._col_from_ids(
                    ids,
                    variable=variable,
                    dim_col=dim_col,
                    level=level_a,
                    column="Positive",
                    row_filters=row_filters,)
                record["pos_count"] = int(ref_pos.sum())
            except TypeError:
                pass

        # save results
        self.results.append(record)
        return record

    def batch_fit(self, ref_methods, test_methods, variables, site_filters=None, model="deming", measurement_col="Value"):
        """Run regression across many combinations of methods, variables and sites in one call."""
        # currently ref_methods, test_methods, variables should be lists, consider option to allow strings as well

        ref_methods = ensure_list(ref_methods)
        test_methods = ensure_list(test_methods)
        variables = ensure_list(variables)
        site_filters = ensure_list(site_filters)

        for ref, test, var, sites in product(ref_methods, test_methods, variables, site_filters):
            if ref == test:
                continue
            try:
                self.fit(ref, test, var, site_filter=sites, reg_method=model, measurement_col=measurement_col)
            except Exception as e:
                print(f"Skipping {ref} vs {test} ({var}, {sites}): {e}")

    def batch_compare(self,
                      *,
                      levels_a,
                      levels_b,
                      variables,
                      comp_func='Regression',
                      dim_col='Method',
                      row_filters=None,
                      split_by=None,
                      **kwargs):
        """
        Run pairwise comparisons across many combinations of levels, variables, and strata.

        Parameters
        ----------
        levels_a, levels_b : str or sequence of str
            Levels of dim_col to compare pairwise (levels_a being the reference and levels_b being the test).

        variables : str or sequence of str
            Variables to iterate over (one comparison per variable).

        comp_func : str callable
            Function performing a single comparison (e.g. self.fit, self.sen_spe),
            or name of such function 'deming', 'binary', etc.
            Must accept:
                level_a, level_b, variable, dim_col, row_filters, **kwargs

        dim_col : str
            Column defining the comparison dimension (e.g. "Method", "Investigator").

        row_filters : dict, optional
            Row restriction filters applied *before* stratification.
            Example: {"Method": "OMR", "Investigator": ["MeanInvestigator", "CBM"]}

        split_by : str or sequence of str, optional
            Columns whose unique values define separate comparisons.
            Example: "Site" or ["Site", "Method"].

        kwargs :
            Passed through to compare_func.
        """
        # ---- find callable for wanted function ---
        if not callable(comp_func):
            if comp_func in ['Regression', 'regression', 'deming']:
                comp_func = self.fit
            elif comp_func in ['binary', 'sen_spe']:
                comp_func = self.sen_spe
            elif comp_func in ['confusion_matrix', 'crosstab', 'qualitative']:
                comp_func = self.confusion_matrix
            else:
                raise ValueError(f'{comp_func} not appropriate comparison function')

        # ---- normalize inputs ----
        # if iterated arguments are str or None, convert to list
        levels_a = ensure_list(levels_a)
        levels_b = ensure_list(levels_b)
        variables = ensure_list(variables)
        split_by = ensure_list(split_by)

        # ---- base dataframe after row restriction ----
        df_base = self.df
        if row_filters:
            df_base = self._apply_row_filters(df_base, row_filters)

        # ---- generate strata ----
        if split_by and split_by != [None]:
            grouped = df_base.groupby(split_by, dropna=False)
            strata = list(grouped)
        else:
            strata = [(None, df_base)]

        # ---- main loop ----
        for stratum_key, df_stratum in strata:
            # build stratum-specific filters
            if stratum_key is None:
                stratum_filters = {}
                stratum_label = "All"
            else:
                if len(split_by) == 1:
                    stratum_key = (stratum_key,)
                stratum_filters = dict(zip(split_by, stratum_key))
                stratum_label = format_filter_label(stratum_filters)

            # merge row_filters + stratum filters
            combined_filters = dict(row_filters or {})
            combined_filters.update(stratum_filters)

            for level_a, level_b, var in product(levels_a, levels_b, variables):
                if level_a == level_b:
                    continue
                try:
                    comp_func(level_a, level_b, var, dim_col=dim_col, row_filters=combined_filters, **kwargs)
                except Exception as e:
                    # print(f"[batch_compare] Skipping comparison: "
                    #       f"{dim_col} {level_a} vs {level_b}, "
                    #       f"variable='{var}', "
                    #       f"{stratum_label}. "
                    #       f"Reason: {e}")
                    continue


    def calc_all_biases(self, crit_points_dict):
        """
        Go through all results and calculate and save biases for all critical points
        Args:
            crit_points_dict: dict of variable names as keys and lists of numbers as values (like the one in a MetadataBundle)
        Returns:
            in each member of the results attribute which has critical points, add a dictionary where the points are
            keys and bias details are values
        """
        if self.results == {}:
            print('No results yet - calculate regressions first and then calculate biases')
            return
        for rec in self.results:
            reg = rec.get("reg")
            if reg is None:
                continue

            var_name = rec.get("variable")
            crit_points = crit_points_dict.get(var_name, [])
            if not crit_points:
                continue
            biases = sb.calc_bias_at_points(reg, crit_points)
            rec["biases"] = biases


    # def regressions_to_dataframe(self) -> pd.DataFrame:
    #     """
    #     Convert regressions in MethodComparator.results to a DataFrame.
    #     Expects results argument to be a dict with {key: {..., "reg": RegressionResult}}.
    #     """
    #     rows = []
    #     for key, val in self.results.items():
    #         ref_method, test_method, variable, site = key
    #         stats: RegressionResult = val["reg"]
    #
    #         # create printable version of regression results (to do: add this as attribute for RegressionResult)
    #         if math.isnan(stats.slope) or stats.slope == float('inf') or stats.slope == float('-inf'):
    #             reg_strngs = ["NA", "NA", "NA"]
    #         else:
    #             reg_strngs = []
    #             for param in ('slope', 'intercept'):
    #                 pnt_est = getattr(stats, param)
    #                 bottom = getattr(stats, f'{param}_ci')[0]
    #                 top = getattr(stats, f'{param}_ci')[1]
    #                 reg_strngs.append(f"{pnt_est:.2f}\n({bottom:.2f}-{top:.2f})")
    #             reg_strngs.append(f"{stats.r2:.2f}")
    #
    #         row = {
    #             "Ref Method": ref_method,
    #             "Test Method": test_method,
    #             "Variable": variable,
    #             "Site": site,
    #             "N": stats.n,
    #             "Slope": reg_strngs[0],
    #             "Intercept": reg_strngs[1],
    #             "R^2": reg_strngs[2],
    #             "Positives": val.get('pos_count', None)
    #         }
    #
    #         rows.append(row)
    #
    #     return pd.DataFrame(rows)

    def regressions_to_dataframe(self) -> pd.DataFrame:
        """
        Convert regressions in MethodComparator.results to a DataFrame.
        Expects results attribute to be a dict with {key: {..., "reg": RegressionResult}}.
        """
        rows = []
        for rec in self.results:
            stats: RegressionResult = rec.get("reg")
            if stats is None:
                continue
            if math.isnan(stats.slope) or not np.isfinite(stats.slope):
                slope_str = intercept_str = r_str = "NA"
            else:
                slope_str = f"{stats.slope:.2f}\n({stats.slope_ci[0]:.2f}-{stats.slope_ci[1]:.2f})"
                intercept_str = f"{stats.intercept:.2f}\n({stats.intercept_ci[0]:.2f}-{stats.intercept_ci[1]:.2f})"
                r_str = f"{stats.r:.2f}"

            row = self._base_result_row(rec)
            reg_results = {
                "N": stats.n,
                "Positives": rec.get("pos_count"),
                "Slope": slope_str,
                "Intercept": intercept_str,
                "r": r_str,
            }
            row.update(reg_results)
            rows.append(row)
        return pd.DataFrame(rows)


    def biases_to_dataframe(self) -> pd.DataFrame:
        """
        Convert biases in MethodComparator.results to a DataFrame.
        Expects results attribute to be a dict with {key: {..., "biases": dict}}.
        """
        rows = []
        for rec in self.results:
            biases = rec.get("biases", {})
            if not biases:
                continue
            for point, bias in biases.items():
                bias_strngs = []
                for bias_type in ('abs', 'rel'):
                    pnt_est = bias[f'{bias_type}_bias']
                    bottom = bias[f'{bias_type}_bias_ci'][0]
                    top = bias[f'{bias_type}_bias_ci'][1]
                    if math.isnan(pnt_est):
                        bias_str = "NA"
                    elif math.isnan(bottom) or math.isnan(top):
                        bias_str = f"{pnt_est:.2f}"
                    else:
                        bias_str = f"{pnt_est:.2f}\n({bottom:.2f}-{top:.2f})"
                    bias_strngs.append(bias_str)

                row = self._base_result_row(rec)
                results = {
                    "Critical Point": point,
                    "Bias": bias_strngs[0],
                    "%Bias": bias_strngs[1]
                }
                row.update(results)
                rows.append(row)
        return pd.DataFrame(rows)

    def metrics_to_dataframe(self) -> pd.DataFrame:
        """
        Convert MethodComparator.metrics to a DataFrame.
        Expects metrics attribute to be a dict with
        {key: {...,
               "mtr":{
                    "sensitivity": {"value": ..., "ci": (x1, x2)}, ...
                    "tp": ...}}}.
        """
        rows = []
        for rec in self.metrics:
            mtrcs = rec["mtr"]

            mtrcs_strngs = []
            for mtr_name in ['sensitivity', 'specificity', 'agreement']:
                mtr_val = mtrcs.get(mtr_name, np.nan)

                # Check if metric is the bootstrapped output format: {"value": X, "ci": (L, U)}

                if isinstance(mtr_val, dict):
                    pnt_est = mtr_val.get('value', np.nan)
                    if pd.isna(pnt_est):
                        mtrcs_strngs.append("NA")
                    else:
                        ci = mtr_val.get('ci', (np.nan, np.nan))
                        pnt_est_pct = pnt_est * 100
                        bottom_pct = ci[0] * 100
                        top_pct = ci[1] * 100
                        mtrcs_strngs.append(f"{pnt_est_pct:.1f}\n({bottom_pct:.1f}-{top_pct:.1f})")

                # Otherwise, treat it as the standard binary output (scalar)
                else:
                    if pd.isna(mtr_val):
                        mtrcs_strngs.append("NA")
                    else:
                        pnt_est_pct = mtr_val * 100
                        mtrcs_strngs.append(f"{pnt_est_pct:.1f}")

            row = self._base_result_row(rec)
            results = {
                "N": len(rec["x"]),
                "Sensitivity": mtrcs_strngs[0],
                "Specificity": mtrcs_strngs[1],
                "Agreement": mtrcs_strngs[2],
                "TP": mtrcs["tp"], "TN": mtrcs["tn"], "FN": mtrcs["fn"], "FP": mtrcs["fp"]}
            row.update(results)
            rows.append(row)
        return pd.DataFrame(rows)

    def crosstabs_to_dataframe(self) -> pd.DataFrame:
        """
        Convert confusion matrices in MethodComparator.metrics to a flat DataFrame.
        Expects metrics attribute to contain records with an "xtab" key.
        """
        rows = []
        for rec in self.metrics:
            xtab = rec.get("xtab")
            if xtab is None:
                continue

            base_row = self._base_result_row(rec)

            # Melt the crosstab to extract Ref Class, Test Class, and Count
            # This turns a 3x3 matrix into 9 rows of data
            melted = xtab.reset_index().melt(id_vars=xtab.index.name,
                                             var_name="Test Grade",
                                             value_name="Count")

            # Rename the reference column for clarity
            melted.rename(columns={xtab.index.name: "Ref Grade"}, inplace=True)

            for _, match_row in melted.iterrows():
                row = base_row.copy()
                row.update({
                    "Ref Grade": match_row["Ref Grade"],
                    "Test Grade": match_row["Test Grade"],
                    "Count": match_row["Count"],
                    "Agreement": "Match" if match_row["Ref Grade"] == match_row["Test Grade"] else "Mismatch"
                })
                rows.append(row)

        return pd.DataFrame(rows)

    def export_confusion_matrices_to_excel(self, filepath: str):
        """
        Export qualitative confusion matrices to a single Excel sheet, stacked vertically.
        As requested, Reference methods are columns, Test methods are rows.
        """
        # Ensure the output is an Excel file, as CSV doesn't support this layout cleanly
        if not filepath.lower().endswith(('.xlsx', '.xls')):
            filepath += '.xlsx'

        # Use ExcelWriter to append multiple matrices to the same sheet
        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            start_row = 0
            sheet_name = 'Confusion Matrices'

            for rec in self.metrics:
                xtab = rec.get("xtab")
                if xtab is None:
                    continue

                # 1. Format the metadata title for this matrix
                var_name = rec["variable"]
                level_a = rec["level_a"]
                level_b = rec["level_b"]

                # Use your existing formatting utility
                stratum_str = format_filter_label(rec.get("stratum", {}))

                title = f"Variable: {var_name} | Ref: {level_a} vs Test: {level_b}"
                if stratum_str:
                    title += f" | {stratum_str}"

                # 2. Transpose the matrix to make Ref=Columns and Test=Rows
                matrix_to_export = xtab.T
                matrix_to_export.index.name = f"Test ({level_b})"
                matrix_to_export.columns.name = f"Ref ({level_a})"

                # 3. Write the title string (converted to a 1x1 DataFrame for easy writing)
                pd.DataFrame([title]).to_excel(writer, sheet_name=sheet_name,
                                               startrow=start_row, startcol=0,
                                               header=False, index=False)

                # 4. Write the transposed matrix right below the title
                matrix_to_export.to_excel(writer, sheet_name=sheet_name,
                                          startrow=start_row + 1, startcol=0)

                # 5. Increment start_row for the next matrix
                # (1 row for title + matrix height + 1 row for columns + 2 rows of blank padding)
                start_row += matrix_to_export.shape[0] + 4

        print(f"Successfully exported confusion matrices to {filepath}")


    @staticmethod
    def _base_result_row(rec: dict) -> dict:
        """
        Build the common identifying columns for a result record.
        """
        comp_dim = rec["comparison_dim"]
        row = {
            f"Ref {comp_dim}": rec["level_a"],
            f"Test {comp_dim}": rec["level_b"],
            "Variable": rec["variable"],
        }
        stratum = normalize_filter_mapping(rec.get("stratum", {}))
        row.update(stratum)
        return row

    def save_results(self, filepath: str, result_type: str = "reg", fill_all=True) -> None:
        """
        Save MethodComparator.results to CSV or Excel.
        """
        if result_type in ["matrix_visual", "stacked_xtab"]:
            self.export_confusion_matrices_to_excel(filepath)
            return

        if result_type in ["reg", "Regression"]:
            df = self.regressions_to_dataframe()
        elif result_type in ["bias", "Bias"]:
            df = self.biases_to_dataframe()
        elif result_type in ["senspe", "SenSpe", "binary", "Binary", "sen_spe"]:
            df = self.metrics_to_dataframe()
        elif result_type in ['confusion_matrix', 'Confusion Matrix']:
            df = self.crosstabs_to_dataframe()
        else:
            raise ValueError(f"{result_type} result_type not supported")

        # fill in 'All' value for places where Site or Investigator field is blank
        for col in ['Site', 'Investigator']:
            if col in df.columns:
                df[col].fillna('All', inplace=True)

        if isinstance(filepath, (str, os.PathLike)):
            sb.write_df_to_file(df, filepath)
        else:
            return df

    def plot_all_regressions(self, pdf_path: str, *, style: Dict = None, scatter_kwargs: Dict = None, overlay_kwargs: Dict = None):
        # to do: add option where pdf_path=None, in which case plots are shown directly
        """
        Create one page per result: scatter + regression overlay; save to a single PDF.

        Parameters
        ----------
        pdf_path : str, optional
            Output path for the PDF (e.g., 'all_regressions.pdf').
            to do - have option to leave empty so plots only shown.
        style : dict, optional
            Base style applied to all pages. Can include template strings like:
              title="Regression for {name}", xlabel="{x_label}", ylabel="{y_label}"
            Common useful keys:
              - legend (bool), legend_loc, legend_title
              - grid (bool), tight_layout (bool)
              - equal_limits (bool), pad_limits (bool), xpad, ypad, pad_mode
              - palette, line_color, scatter_color, ci (bool), ci_mode ("lines"/"shade")
        scatter_kwargs : dict, optional
            Extra kwargs to pass to plot_scatter_basic (rarely needed; style usually suffices).
        overlay_kwargs : dict, optional
            Extra kwargs to pass to overlay_regression_line (e.g., n_points=400).
        """
        base_style = {'gird': True,
                      'equal_limits': True,
                      # regression overlay defaults:
                      'ci': True,
                      'ci_mode': 'shade'}
        user_style = style or {}
        scatter_kwargs = scatter_kwargs or {}
        overlay_kwargs = overlay_kwargs or {}

        with PdfPages(pdf_path) as pdf:
            for rec in self.results:
                # ---- labels ----
                comp_dim = rec["comparison_dim"]
                xname = rec["level_a"]
                yname = rec["level_b"]
                varname = rec["variable"]

                label = format_filter_label(rec.get("stratum", {}))
                plot_title = f"{varname}\n{label}" if label else varname
                scpc_style = {'title': plot_title, 'xlabel': xname, 'ylabel': yname}
                final_style = plotting.merge_styles(base_style, user_style, scpc_style)

                # ---- plotting ----
                fig, ax = plt.subplots(figsize=final_style.get("figsize", (6, 5)))
                fig, ax = plotting.plot_scatter_basic(rec, style=final_style, fig=fig, ax=ax, **scatter_kwargs)
                plotting.overlay_regression_line(data=rec, fig=fig, ax=ax, style=final_style,
                                                 **overlay_kwargs)
                pdf.savefig(fig)
                plt.close()


            # for data, rslt in self.results.items():
            #     xname, yname, varname, site = data
            #     if rslt['reg'].r is not None:
            #         plot_title = varname if site == 'All' else f'{varname}\n{site}'
            #         scpc_style = {'title': plot_title, 'xlabel': xname, 'ylabel': yname}
            #         final_style = plotting.merge_styles(base_style, user_style, scpc_style)
            #
            #         fig, ax = plt.subplots(figsize=final_style.get("figsize", (6, 5)))
            #         fig, ax = plotting.plot_scatter_basic(rslt, style=final_style, fig=fig, ax=ax, **scatter_kwargs)
            #         plotting.overlay_regression_line(fig=fig, ax=ax, result=rslt['reg'], style=final_style, **overlay_kwargs)
            #         pdf.savefig(fig)
            #         plt.close()
