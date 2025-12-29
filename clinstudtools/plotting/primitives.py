""" low-level, stateless plotting functions drawing single element onto given axis """
import matplotlib.pyplot as plt
import numpy as np
from primitives_helpers import normalize_data

def plot_scatter_basic(data, *, style=None, fig=None, ax=None, **_):
    """
       Draw a basic x/y scatter plot.

       Parameters
       ----------
       data :
           One of:
           - dict-like with keys 'x' and 'y'
           - pandas DataFrame with columns 'x' and 'y'
           - tuple (x, y)
       fig, ax :
           Optional matplotlib Figure and Axes to draw on.
           If not provided, a new figure and axis are created.
       style : dict, optional
           Visual styling options (marker, alpha, color, grid, etc.):
            - marker: marker style (default: 'o')
            - markersize: float
            - alpha: float (transparency)

       Returns
       -------
       fig, ax :
           The matplotlib Figure and Axes used for drawing.
       """
    style = style or {}

    # normalize axis
    if ax is None or fig is None:
        fig, ax = plt.subplots()

    xy = normalize_data(data, required=("x", "y"))
    x, y = xy["x"], xy["y"]
    ax.scatter(x, y,
               color=style.get("scatter_color"),
               marker=style.get("marker", "o"),
               alpha=style.get("alpha", 0.7),
               label=style.get("label"))
    # apply_axes_style(ax, style)
    # apply_figure_style(fig, style)
    return fig, ax


