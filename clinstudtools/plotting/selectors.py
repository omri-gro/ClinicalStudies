""" result -> data selection logic """
import sys
sys.path.append(r'C:\Users\omrig\DataAnalysisProjects\ClinicalStudies\clinstudtools')
from utils import expect_single


def select_records(results, **criteria):
    """
    Filter a list of result objects by simple equality criteria.

    Parameters
    ----------
    results : iterable
        Collection of result-like objects (usually list of dicts or list of objects that have attributes similar to the kwargs).
    **criteria :
        Key-value pairs that must match exactly on each result (kwargs).

    Returns
    -------
    list
        All results matching the given criteria.
    """
    out = []
    for r in results:
        if all(r.get(k) == v for k, v in criteria.items()):
            out.append(r)
    return out


""" selector factories """
def select_method_comparison_xy(ref_method, test_method):
    def resolver(results):
        recs = select_records(results, kind="method_comparison", ref_method=ref_method, test_method=test_method)
        rec = expect_single(recs, "method comparison")
        df = rec["data"]
        return {"x": df["ref_mean"], "y": df["test"]}

    return resolver

def select_regression(ref_method, test_method):
    def resolver(results):
        recs = select_records(results, kind="method_comparison", ref_method=ref_method, test_method=test_method)
        rec = expect_single(recs, "regression")
        return rec["reg"]
    return resolver

def select_reference_agreement_bars(method):
    def resolver(results):
        recs = select_records(results, kind="inter_reviewer_agreement", method=method)
        rec = expect_single(recs, "reference agreement")
        df = rec["data"]
        return {"x": df["ref_mean"], "y_min": df["investigator_1"], "y_max": df["investigator_2"]}
    return resolver
