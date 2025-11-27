""" consider splitting into multiple files later """

import sandbox as sb
import pandas as pd
from typing import List, Tuple, Optional, Dict, Any, Union, Sequence
from pathlib import Path
from itertools import product
import os
import numpy as np
import math
import matplotlib.pyplot as plt

# placeholder nutil integrating regression functions
import sys
sys.path.append(r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification')
import reg_types as reg
from dataclasses import dataclass, asdict
from matplotlib.backends.backend_pdf import PdfPages
import plotting
from pipelines import medium_pipe, bma_prep_pipeline


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
            r=data.get("r") or data.get("R2"),
            r2=data.get("r2") or data.get("correlation_coefficient"),
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
        # stores regression results by (ref, test, variable, site)
        self.results = {}
        # stores qualitative and semi-quantitative metrics (sensitivity, kappa, etc.) by (ref, test, variable, site)
        self.metrics = {}

    @classmethod
    def from_paths_dict(cls, paths: dict, metadata: sb.MetadataBundle, dir=None, measurement_col='Value',
                        more_id_vars=None, stnrd_id=True, bma=False):
        # maybe don't use, might be better to create something like that every time due to all edge cases

        # function that gets dict with paths as values and descriptors of method/site/etc as keys
        # will read, prepare, concatenate and convert to MethodComparator
        # assumes all files are in same directory

        # to do: clean and add error messages
        # to do: use metadata to pass info like stnrd_id?
        # to do: make it possible for the values in paths dict to be a list of multiple paths, instead of a single path

        # list of columns that will not be pivoted
        id_vars = ["SampleID", "Site", "Method", "FileName", 'Investigator'] if bma else ["SampleID", "Site", "Method", "FileName"]
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
                    df = bma_prep_pipeline(path, site, method, metadata, dir=dir, id_vars=id_vars)
                else:
                    df = medium_pipe(path, site, method, metadata, dir=dir, id_vars=id_vars, stnrd_id=stnrd_id)
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
    def filter_by_df(self,
                     filtering_source: Union[str, Path, pd.DataFrame],
                     filtering_cols: Optional[Sequence[str]] = None,
                     include_rows=False):
        """
        Using a dataframe of samples to exclude/include (or a file containing one), create new MethodComparator only without/with those rows.
        filtering_cols - Sequence of column names to use for filtering (e.g., ['Site', 'SampleID', 'Investigator', 'Method']),
                         else infer from the filtering df.
        include_rows - Change to True so resulting df will have only samples included in the filtering_source, instead of those that don't.
        """
        df2 = filtering_source if isinstance(filtering_source, pd.DataFrame) else sb.raw_to_df(filtering_source)
        if not filtering_cols:
            filtering_cols = set(df2.columns)
            # if want to use "FileName" as one of filtering columns, define filtering_cols directly and use DataFrame for filtering_source
            filtering_cols.discard("FileName")

        df1 = self.df.copy()
        common = list(set(df1.columns) & filtering_cols)  # columns existing in both

        df1_com = df1[common]
        df2_com = df2[common]

        # Build a set-like MultiIndex of unique df2 keys, test membership for df1
        keys_df2 = pd.MultiIndex.from_frame(df2_com.drop_duplicates())
        mask = pd.MultiIndex.from_frame(df1_com).isin(keys_df2)  # samples appearing in both are in True
        mask = ~mask if not include_rows else mask
        out_df = df1.loc[mask].copy()

        # return a new MethodComparator with filtered rows
        return MethodComparator(out_df, self.measurement_col)

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
        else:
            wide_df = sb.to_comparison_matrix(self, metadata=getattr(self, "metadata", None), **kwargs)

        if out_path:
            sb.write_df_to_file(wide_df, out_path)
        return wide_df

    def clean_calculations(self):
        self.results = {}
        self.metrics = {}

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

        # note this function actually does same thing as _prepare_arrays method from MethodComparator
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

    def fit(self,
            ref_method: str,
            test_method: str,
            variable: str,
            site_filter: Optional[List[str]] = None,
            model: Optional[str] = "deming",
            measurement_col="Value"):
        """Run regression for one variable/site/method pair."""

        x, y, ids = self._prepare_arrays(ref_method, test_method, variable, site_filter=site_filter, measurement_col=measurement_col, on_duplicates="first")

        if site_filter:
            site_filter = sb._ensure_list(site_filter)

        # placeholder regression_func → replace with your real implementation
        stats = reg.regression_comp(x, y, reg_method=model)  # placeholder - replace later with different regression functions
        if type(stats) == RegressionResult:
            result = stats
            if result.n is None:
                result.n = len(x)
        elif type(stats) == dict:  # backwards compatibility - for old regression functions returning dicts
            result = RegressionResult.from_dict(stats)
        else:
            result = RegressionResult(**stats)

        key = (ref_method, test_method, variable, ",".join(site_filter) if site_filter else "All")
        self.results[key] = {
            "x": x,
            "y": y,
            "ids": ids,
            "reg": result
        }
        return self.results[key]

    def batch_fit(self, ref_methods, test_methods, variables, site_filters=None, model="deming", measurement_col="Value"):
        """Run regression across many combinations of methods, variables and sites in one call."""
        # currently ref_methods, test_methods, variables should be lists, consider option to allow strings as well

        ref_methods = sb._ensure_list(ref_methods)
        test_methods = sb._ensure_list(test_methods)
        variables = sb._ensure_list(variables)
        site_filters = sb._ensure_list(site_filters)

        for ref, test, var, sites in product(ref_methods, test_methods, variables, site_filters):
            if ref == test:
                continue
            try:
                self.fit(ref, test, var, site_filter=sites, model=model, measurement_col=measurement_col)
            except Exception as e:
                print(f"Skipping {ref} vs {test} ({var}, {sites}): {e}")

    def batch_compare(self, ref_methods, test_methods, variables,
                      function='Regression', site_filters=None, **kwargs):
        """Run comparison calculation across many combinations of methods, variables and sites in one call.
        ref_methods, test_methods, variables, site_filters - str, list or None
        function - 'deming',
        kwargs - sent to calculating function (model, measurement, etc.)
        """
        # if iterated arguments are str or None, convert to list
        ref_methods = sb._ensure_list(ref_methods)
        test_methods = sb._ensure_list(test_methods)
        variables = sb._ensure_list(variables)
        site_filters = sb._ensure_list(site_filters)

        for ref, test, var, sites in product(ref_methods, test_methods, variables, site_filters):
            if ref == test:
                continue
            try:
                if function == 'regression':   # consider using register decorator if many functions
                    self.fit(ref, test, var, site_filter=sites, **kwargs)
                elif function == 'binary':
                    continue
            except Exception as e:
                print(f"Skipping {function} calculation for {var} ({ref} vs {test}, {sites}): {e}")


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
        for key, val in self.results.items():
            var_name = key[2]
            crit_points = crit_points_dict.get(var_name, [])
            if not crit_points:
                continue
            biases = sb.calc_bias_at_points(val["reg"], crit_points)
            self.results[key]["biases"] = biases


    def regressions_to_dataframe(self) -> pd.DataFrame:
        """
        Convert regressions in MethodComparator.results to a DataFrame.
        Expects results argument to be a dict with {key: {..., "reg": RegressionResult}}.
        """
        rows = []
        for key, val in self.results.items():
            ref_method, test_method, variable, site = key
            stats: RegressionResult = val["reg"]

            # create printable version of regression results (to do: add this as attribute for RegressionResult
            if math.isnan(stats.slope) or stats.slope == float('inf') or stats.slope == float('-inf'):
                reg_strngs = ["NA", "NA", "NA"]
            else:
                reg_strngs = []
                for param in ('slope', 'intercept'):
                    pnt_est = getattr(stats, param)
                    bottom = getattr(stats, f'{param}_ci')[0]
                    top = getattr(stats, f'{param}_ci')[1]
                    reg_strngs.append(f"{pnt_est:.2f}\n({bottom:.2f}-{top:.2f})")
                reg_strngs.append(f"{stats.r2:.2f}")

            row = {
                "ref_method": ref_method,
                "test_method": test_method,
                "Variable": variable,
                "Site": site,
                "N": stats.n,
                "Slope": reg_strngs[0],
                "Intercept": reg_strngs[1],
                "R^2": reg_strngs[2]
            }
            rows.append(row)

        return pd.DataFrame(rows)

    def biases_to_dataframe(self) -> pd.DataFrame:
        # to do - add argument here and to save_results about whether to save variables like "Sites" and "Test Method"
        # to do (long term) - merge this and results_to_dataframe functions
        """
        Convert MethodComparator.results to a DataFrame.
        Expects results argument to be a dict with {key: {..., "biases": dict}}.
        """
        results = self.results
        rows = []
        for (ref, test, var, sites), res in results.items():
            biases = res.get("biases", {})
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
                row = {
                        "Reference Method": ref,
                        "Test Method": test,
                        "Variable": var,
                        "Sites": sites,
                        "Critical Point": point,
                        "Bias": bias_strngs[0],
                        "%Bias": bias_strngs[1]
                        }
                rows.append(row)

        return pd.DataFrame(rows)

    def save_results(self, filepath: str, txt_form: bool = True, result_type: str = "reg") -> None:
        """
        Save MethodComparator.results to CSV or Excel.
        """
        if result_type in ["reg", "Regression"]:
            df = self.regressions_to_dataframe()
        elif result_type in ["bias", "Bias"]:
            df = self.biases_to_dataframe()
        else:
            raise ValueError(f"{result_type} result_type not supported")
        sb.write_df_to_file(df, filepath)

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
            for data, rslt in self.results.items():
                xname, yname, varname, site = data
                if rslt['reg'].r is not None:
                    plot_title = varname if site == 'All' else f'{varname}\n{site}'
                    scpc_style = {'title': plot_title, 'xlabel': xname, 'ylabel': yname}
                    final_style = plotting.merge_styles(base_style, user_style, scpc_style)

                    fig, ax = plt.subplots(figsize=final_style.get("figsize", (6, 5)))
                    fig, ax = plotting.plot_scatter_basic(rslt, style=final_style, fig=fig, ax=ax, **scatter_kwargs)
                    plotting.overlay_regression_line(fig=fig, ax=ax, result=rslt['reg'], style=final_style, **overlay_kwargs)
                    pdf.savefig(fig)
                    plt.close()

if __name__ == "__main__":
    analysis_name = 'bwh100_investigator_training'
    suffix = '_old_version'
    # suffix = ''
    save_name = f"{analysis_name}{suffix}"
    site = 'BWH'

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    clv_df_raw = sb.raw_to_df('BWH_ClV%.csv', site, 'Manual', dir=raw_dir)
    mnl_df_raw = sb.raw_to_df('BWH_manual%.csv', site, 'Manual', dir=raw_dir)
    cbm_df_raw = sb.raw_to_df(f'BWH_CBM%{suffix}.csv', site, 'CBM', dir=raw_dir)

    clv_df = sb.short_pipe(clv_df_raw, metadata)  # short_pipe moved to pipelines.py
    mnl_df = sb.short_pipe(mnl_df_raw, metadata)
    cbm_df = sb.short_pipe(cbm_df_raw, metadata)

    bwh_df = pd.concat([clv_df, mnl_df, cbm_df])

    bwh_df = bwh_df[bwh_df["Value"].notna()]  # remember to include in larger pipeline later

    methd_comp = MethodComparator(bwh_df)

    methd_comp.fit('Manual', 'CBM', 'Segmented Neutrophil')
    methd_comp.calc_all_biases(metadata.crit_points)
    methd_comp.save_results(r'results/seg_reg.csv')
    methd_comp.save_results(r'results/seg_reg.xlsx', result_type='bias')


    print('')
