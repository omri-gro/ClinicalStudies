import numpy as np

def normalize_data(data, *, required):
    """
    Normalize plotting input into named numpy arrays.

    Parameters
    ----------
    data :
        dict-like, DataFrame, or tuple
    required : tuple of str
        Names of required fields, e.g. ('x', 'y') or ('x', 'y_min', 'y_max')

    Returns
    -------
    dict[str, np.ndarray]
    """
    if isinstance(data, dict):
        out = {k: data[k] for k in required}
    elif hasattr(data, "__dataframe__"):
        out = {k: data[k] for k in required}
    else:
        out = dict(zip(required, data))

    return {k: np.asarray(v) for k, v in out.items()}

