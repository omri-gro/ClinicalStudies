# core - project classes (except MethodComparator)
import pandas as pd
from .utils import load_yaml

class MetadataBundle:  # to do: add attribute for thresholds & grades?
    #  Holds all metadata and configuration for the pipeline.
    def __init__(self, metadata_path):
        # consider a simpler init without metadata_path, and move current content into from_yaml method
        context = load_yaml(metadata_path)
        self.variables = context.get("variables", {})
        self.alias_map = self.build_alias_map(self.variables)  # not the standard way of creating variables, could have just defined self.alias_map within the function
        self.variable_groups = self.build_variable_groups(self.variables)
        self.crit_points = self.build_lists_map(self.variables, keyword="crit_points")  # to do: add functionality where critical points are defined by grading thresholds if one is provided and the other isn't.
        self.normal_ranges = self.build_lists_map(self.variables, keyword="normal_range")  # notice this currently doesn't directly require exactly 2 values in normal_range
        self.grading_specs = self.building_grading_specs(self.variables)

        self.src_fixes = self.build_src_fixes(context)

        self.pregraded_index = self.build_pregraded_index()

    @staticmethod
    def build_alias_map(variables):
        alias_map = {}
        for canonical, props in variables.items():
            alias_map[canonical] = canonical  # map canonical to itself
            for alias in props.get("aliases", []):
                alias_map[alias] = canonical
        return alias_map

    @staticmethod
    def build_variable_groups(variables):
        """
        Precompute all variable groups as a dict: group_name → list of variable names.
        """
        group_map = {}
        for var, props in variables.items():
            for group in props.get("groups", []):
                group_map.setdefault(group, []).append(var)
        return group_map

    @staticmethod
    def build_lists_map(variables, keyword):
        """
        Create dictionary of all variables where the keyword holds a list (of critical points, ranges, etc.)
        """
        points_map = {}
        for var, props in variables.items():
            crit_points = props.get(keyword)
            if isinstance(crit_points, list):
                points_map[var] = crit_points
        return points_map

    @staticmethod
    def building_grading_specs(variables):
        """
        Create a dictionary of variable_name -> dict(thresholds, grades, right_closed, clamp_out_of_range),
        where only variables that have thresholds+grades appear.
        """
        grading_specs = {
            name: spec for name, spec in variables.items()
            if isinstance(spec, dict)
            and spec.get("thresholds") is not None
            and spec.get("grades") is not None
        }
        return grading_specs

    def build_src_fixes(self, context):
        """
        Create a dictionary of (site, method) -> {"rule name": ...}
        """
        fix_dict = context.get("Source fixes", {})
        src_fixes = {}
        for k, v in fix_dict.items():
            site, method = self.parse_site_method_key(k)
            src_fixes[(site, method)] = v or {}
        return src_fixes

    def build_pregraded_index(self):  # to do: consider applying dimensionless approach for indices as well
        """
        Create a Series where True if index (Site, Method, Variable) is provided as grade already
        Within the src_fixes attribute, this can be specified for variables or groups of variables
        """
        rows = []
        for (site, method), rule in self.src_fixes.items():
            for var in (rule.get("raw_given_as_grade") or []):
                if var in self.variables.keys():
                    rows.append((site, method, var))
                else:
                    for var_name in self.variable_groups.get(var, []):
                        rows.append((site, method, var_name))
        if rows:
            idx = pd.MultiIndex.from_tuples(rows, names=["Site", "Method", "Variable"])
            return pd.Series(True, index=idx, dtype=bool)
        else:
            return None

    @staticmethod
    def parse_site_method_key(key: str):
        if isinstance(key, tuple):
            site, method = key
        else:
            k = key.strip()
            if k.startswith("(") and k.endswith(")"):
                k = k[1:-1]
            site, method = [p.strip() for p in k.split(",", 1)]
        return site, method
