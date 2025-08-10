""" consider splitting into multiple files later """

import sandbox as sb
import pandas as pd
from typing import List, Tuple, Optional, Dict, Any
from itertools import product
import os
import numpy as np

# placeholder nutil integrating regression functions
import sys
sys.path.append(r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification')
import reg_types as reg
from dataclasses import dataclass, asdict


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
    def __init__(self, df: pd.DataFrame, measurement_col='Measurement'):
        """
        Initialize comparator with long-format dataframe.
        Expected columns: SampleID, Site, Method, Variable, Value
        """
        self.df = df.copy()
        self.measurement_col = measurement_col
        self.results = {}  # stores regression results by (ref, test, variable, site)

    def _prepare_arrays(
        self,
        ref_method: str,
        test_method: str,
        variable: str,
        measurement_col: Optional[str] = 'Value',
        site_filter: Optional[List[str]] = None
    ) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Internal helper: return matched (x, y, ids).

        Parameters:
        ref_method: name of reference review method (e.g., 'OMR')
        test_method: name of test review method (e.g., 'DSS')
        variable: variable to filter on (e.g., 'Total Neutrophil')
        site_filter: list or set of sites to include (optional)
        measurement_col: column name holding the measurement values

    Returns:
        x: array of reference method measurements
        y: array of test method measurements
        ids: list of (Site, SampleID) tuples corresponding to each x/y pair
        """
        # note this function actually does same thing as _prepare_arrays method from MethodComparator
        subset = self.df[self.df["Variable"] == variable].copy()
        if site_filter is not None:
            subset = subset[subset["Site"].isin(site_filter)]

        # Pivot: SampleID+Site as index, Methods as columns
        pivoted = subset.pivot_table(
            index=["SampleID", "Site"],
            columns="Method",
            values=measurement_col
        )

        # Drop missing pairs (cases where at least one of the methods has NaN)
        pivoted = pivoted.dropna(subset=[ref_method, test_method])

        x = pivoted[ref_method].values
        y = pivoted[test_method].values
        ids = pivoted.index.get_level_values("SampleID").values
        # To do: add an option for ids to include the site as well

        return x, y, ids

    def fit(self,
            ref_method: str,
            test_method: str,
            variable: str,
            site_filter: Optional[List[str]] = None,
            model: Optional[str] = "deming"):
        """Run regression for one variable/site/method pair."""

        x, y, ids = self._prepare_arrays(ref_method, test_method, variable, site_filter=site_filter)

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

        key = (ref_method, test_method, variable, tuple(site_filter) if site_filter else "All")
        self.results[key] = {
            "x": x,
            "y": y,
            "ids": ids,
            "reg": result
        }
        return self.results[key]

    def batch_fit(self, ref_methods, test_methods, variables, site_filters=None, model="deming"):
        """Run regression across many combinations in one call."""
        if site_filters is None:
            site_filters = [None]

        for ref, test, var, sites in product(ref_methods, test_methods, variables, site_filters):
            if ref == test:
                continue
            try:
                self.fit(ref, test, var, site_filter=sites, model=model)
            except Exception as e:
                print(f"Skipping {ref} vs {test} ({var}, {sites}): {e}")

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


    def regressions_to_dataframe(self) -> pd.DataFrame:  # really referring only to regression results
        """
        Convert regressions in MethodComparator.results to a DataFrame.
        Expects results argument to be a dict with {key: {..., "reg": RegressionResult}}.
        """
        rows = []
        for key, val in self.results.items():
            ref_method, test_method, variable, site = key
            stats: RegressionResult = val["reg"]

            row = {
                "ref_method": ref_method,
                "test_method": test_method,
                "variable": variable,
                "site": site,
                **stats.as_dict()
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
                    if pnt_est == np.nan:
                        bias_str = "NA"
                    elif bottom == np.nan or top == np.nan:
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
        format = filepath.split(".")[-1]
        if format.lower() == "csv":
            df.to_csv(filepath, index=False)
        elif format.lower() in ("xlsx", "excel"):
            df.to_excel(filepath, index=False)
        else:
            raise ValueError("Format must be 'csv' or 'excel'.")




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

    clv_df = sb.short_pipe(clv_df_raw, metadata)
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
