""" Human-readable representations of comparison metadata """
""" Turning internal comparison descriptors into readable output, formatting consistency, etc. """
""" pure metadata (dicts, lists, etc.) into string transformations (which can be inside a dict or list) """
""" no matplotlib, figures or plotting decisions """

def normalize_filter_mapping(filters: dict) -> dict:
    """
    Return a copy of a filter/stratum mapping with
    single-element containers unwrapped for readability.
    """
    formatted = {}
    for k, v in filters.items():
        if isinstance(v, (list, tuple, set)):
            if len(v) == 1:
                formatted[k] = next(iter(v))
            else:
                formatted[k] = tuple(v)
        else:
            formatted[k] = v
    return formatted

def format_filter_label(filters: dict) -> str:
    if not filters:
        return ""

    formatted = normalize_filter_mapping(filters)
    return ", ".join(f"{k}={v}" for k, v in formatted.items())

