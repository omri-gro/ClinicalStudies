# utils/validation.py

from typing import Union

def expect_single(
    items,
    *,
    what: str = "item",
    context: Union[str, None] = None
):
    """
    Assert that `items` contains exactly one element and return it.
    """
    n = len(items)
    if n != 1:
        msg = f"Expected exactly one {what}, found {n}"
        if context:
            msg += f" ({context})"
        raise ValueError(msg)
    return items[0]
