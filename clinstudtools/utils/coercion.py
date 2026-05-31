# coercion.py
# input normalize: coerce inputs into a predictable shape

# can later add: as_numpy_1d(x), as_numpy_2d(x), as_series(x, name=None), as_index(x),
# ensure_list_or_none (returns None instead of []), is_iterable_but_not_str(x), is_scalar(x), maybe_scalar(x),
# def ensure_dict(x, *, value_if_scalar=True), ensure_mapping(x, keys=None), broadcast_dict(value, keys),
# ensure_columns(x), ensure_str(x)

import pandas as pd
from typing import Any, List, Iterable, Mapping, Optional, Sequence, Union, Tuple
import numpy as np
import warnings


def as_df(
    obj: Any,
    *,
    attr: Optional[str] = "df",
    copy: bool = True,
) -> pd.DataFrame:
    """
    Coerce input into a pandas DataFrame.

    Accepts:
      - pandas.DataFrame
      - pandas.Series (converted via .to_frame())
      - dict (converted via pd.DataFrame(dict))
      - list[dict] (converted via pd.DataFrame(list_of_dicts))
      - an object with attribute `attr` that resolves to one of the above (default: 'df')

    Parameters
    ----------
    obj:
        Input object.
    attr:
        If not None, try to unwrap this attribute first (e.g., 'df').
        Set to None to disable unwrapping.
    copy:
        If True, returns a copy when the result is already a DataFrame.

    Returns
    -------
    pandas.DataFrame
    """
    if attr is not None and not isinstance(obj, (pd.DataFrame, pd.Series)):
        unwrapped = getattr(obj, attr, None)
        if unwrapped is not None:
            obj = unwrapped
    if isinstance(obj, pd.DataFrame):
        return obj.copy() if copy else obj
    if isinstance(obj, pd.Series):
        df = obj.to_frame()
        return df.copy() if copy else df
    # pd.DataFrame(...) already supports dict, list-of-dicts, numpy arrays, etc.
    try:
        df = pd.DataFrame(obj)
    except Exception as e:
        raise TypeError(
            f"Expected DataFrame/Series/dict/list-of-dicts or an object with '{attr}'. "
            f"Got type={type(obj)!r}."
        ) from e
    return df.copy() if copy else df

def ensure_list(x):
    """
    Return x as a list.
    None -> []
    scalar -> [scalar]
    iterable (except str) -> list(x)
    """
    if isinstance(x, list):
        return x
    elif isinstance(x, tuple):
        return list(x)
    elif isinstance(x, pd.Series):
        return x.to_list()
    else:
        return [x]


# Backward-compatible alias
def _ensure_list(*args, **kwargs):
    warnings.warn(
        "utils._ensure_list is deprecated; use utils.coercion.ensure_list",
        DeprecationWarning,
        stacklevel=2,
    )
    return ensure_list(*args, **kwargs)

def _as_df(*args, **kwargs):
    warnings.warn(
        "utils._as_df is deprecated; use utils.coercion.as_df",
        DeprecationWarning,
        stacklevel=2,
    )
    return as_df(*args, **kwargs)




