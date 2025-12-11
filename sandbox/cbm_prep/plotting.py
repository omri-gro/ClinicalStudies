from __future__ import annotations
import json
from pathlib import Path
from typing import Callable, Dict, List, Tuple, Optional

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from plot_styles import DEFAULT_STYLE, PLOT_STYLES

PLOT_FUNCTIONS = {}

"""Core functions"""
def register(name):
    def deco(f):
        PLOT_FUNCTIONS[name] = f
        return f
    return deco

def get_plot(name: str) -> Callable:
    try:
        return PLOT_FUNCTIONS[name]
    except KeyError:
        raise ValueError(f"Unknown plot '{name}'. Available: {list(PLOT_FUNCTIONS)}")

def available_plots() -> List[str]:
    return sorted(PLOT_FUNCTIONS)

def _get_xy(data, xkey="x", ykey="y"):
    # Accept dict-like, pandas DataFrame/Series, or array-likes directly
    # dict-like (such as the ones in a MethodComparator results attr) will have an xkey and ykey
    try:
        # Mapping-like (dict, pandas DataFrame)
        x = data[xkey]
        y = data[ykey]
    except Exception:
        # Assume already passed as a 2-tuple/list: (x, y)
        x, y = data
    return np.asarray(x), np.asarray(y)


"""Helpers"""
def set_equal_limits_and_scale(fig=None, ax=None):
    """
    Sets the x and y axes to have:
    - The same scale (1 unit in x = 1 unit in y)
    - The same limits (xmin == ymin, xmax == ymax)

    Parameters:
    - fig: matplotlib.figure.Figure (optional)
    - ax: matplotlib.axes.Axes (optional)

    If both are None, uses the current active figure.
    """
    import matplotlib.pyplot as plt

    if ax is not None:
        axes = [ax]
    elif fig is not None:
        axes = fig.get_axes()
    else:
        axes = plt.gcf().get_axes()

    for a in axes:
        x_min, x_max = a.get_xlim()
        y_min, y_max = a.get_ylim()

        data_min = min(x_min, y_min)
        data_max = max(x_max, y_max)

        a.set_xlim(data_min, data_max)
        a.set_ylim(data_min, data_max)
        a.set_aspect('equal')


def _deep_merge(a: dict, b: dict) -> dict:
    """Shallow+nested dict merge: values in b override a.
       Combines two dictionaries, but if a key’s value is itself a dictionary, it merges those nested"""
    out = dict(a)
    for k, v in b.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def merge_styles(*layers: dict) -> dict:
    """Allows merging multiple layers of styles in order of increasing priority."""
    out = {}
    for layer in layers:
        out = _deep_merge(out, layer or {})
    return out


def _fmt_template_values(style_or_str, context: dict):
    """If style values are str with {placeholders}, format them with context."""
    if isinstance(style_or_str, str):
        return style_or_str.format(**context)
    if isinstance(style_or_str, dict):
        out = {}
        for k, v in style_or_str.items():
            out[k] = _fmt_template_values(v, context)
        return out
    return style_or_str


def apply_axes_style(ax, style: dict | None = None):
    """Apply common axis styling from a style dict."""
    style = style or {}

    # Labels & title
    if style.get("title") is not None:
        ax.set_title(style["title"])
    if style.get("xlabel") is not None:
        ax.set_xlabel(style["xlabel"])
    if style.get("ylabel") is not None:
        ax.set_ylabel(style["ylabel"])

    # Grid (default True unless explicitly False)
    if style.get("grid", True):
        ax.grid(True)

    # Optional extras if present (no-ops if missing)
    if "xlim" in style: ax.set_xlim(*style["xlim"])
    if "ylim" in style: ax.set_ylim(*style["ylim"])
    if "xscale" in style: ax.set_xscale(style["xscale"])
    if "yscale" in style: ax.set_yscale(style["yscale"])
    if "aspect" in style: ax.set_aspect(style["aspect"])

    # Legend control
    if style.get("legend", False):
        legend_kwargs = {}
        if style.get("legend_loc") is not None:
            legend_kwargs["loc"] = style["legend_loc"]
        if style.get("legend_title") is not None:
            legend_kwargs["title"] = style["legend_title"]
        if style.get("legend_fontsize") is not None:
            legend_kwargs["fontsize"] = style["legend_fontsize"]


def apply_figure_style(fig, style: dict | None = None):
    style = style or {}
    if style.get("tight_layout", True):
        fig.tight_layout()
    if style.get("equal_limits"):
        set_equal_limits_and_scale(fig=fig)


"""Registered plots"""
"""
All of these receive:
 * a data argument which is either:
    a dictionary-like (with 'x' and 'y' as keys) OR
    a pandas DataFrame/Series with 'x' and 'y' columns OR
    a tuple of (x, y), with x and y being array-like or lists of numbers
 * possible style dictionary like the ones in 'plot_styles.py
 * possible fig or ax that were created previously, on which the plot will be drawn (otherwise will create new ones)
Functions return fig and ax
"""
@register("scatter_basic")
def plot_scatter_basic(data, style=None, fig=None, ax=None, **_):
    style = style or {}
    if ax is None or fig is None:
        fig, ax = plt.subplots()
    x, y = _get_xy(data, style.get("xkey", "x"), style.get("ykey", "y"))  # to do: consider finding another way to deal with different xkey instead of 'style'
    ax.scatter(x, y, color=style.get("scatter_color"))
    apply_axes_style(ax, style)
    apply_figure_style(fig, style)
    return fig, ax

@register("scatter_groups")
def plot_scatter_groups(
    data,
    style=None,
    fig=None,
    ax=None,
    *,
    x=None,            # optional 1D array of x values
    y=None,            # optional 1D array of y values
    groups=None,       # optional 1D array of group labels
    xkey="x",          # keys for extracting from `data` if arrays not passed directly
    ykey="y",
    group_key="group",
    legend=True,       # whether to add a legend
    **_,
):
    """
    Scatter plot with points colored by group.

    Parameters
    ----------
    data : dict-like, optional
        Should contain xkey, ykey, and optionally group_key.
    style : dict, optional
        Keys:
          - palette: dict mapping group -> color OR list of colors
          - marker: marker style (default: 'o')
          - markersize: float
          - alpha: float transparency
          - xlabel, ylabel, title, grid: standard styling keys
    fig, ax : existing figure and axes or None to create new
    x, y, groups : arrays, optional
        Directly specify coordinates and grouping. Overrides keys.
    xkey, ykey, group_key : str
        Keys to extract from `data` if x/y/groups not given directly.
    legend : bool
        If True, adds a legend mapping colors to groups.
    """
    style = style or {}
    if ax is None or fig is None:
        fig, ax = plt.subplots()

    # Extract x, y, groups
    if x is None or y is None:
        x = np.asarray(data[xkey])
        y = np.asarray(data[ykey])
    else:
        x = np.asarray(x)
        y = np.asarray(y)

    if groups is None:
        if group_key in data:
            groups = np.asarray(data[group_key])
        else:
            groups = np.full_like(x, fill_value="default", dtype=object)
    else:
        groups = np.asarray(groups)

    # Determine unique groups and colors
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    if n_groups > 6:
        raise ValueError(f"Too many groups ({n_groups}); maximum supported is 6.")

    default_colors = ["C0", "C1", "C2", "C3", "C4", "C5"]
    palette = style.get("palette", default_colors)
    if isinstance(palette, dict):
        color_map = {g: palette.get(g, default_colors[i % 6]) for i, g in enumerate(unique_groups)}
    else:
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(unique_groups)}

    # Plot each group
    for g in unique_groups:
        mask = groups == g
        ax.scatter(
            x[mask],
            y[mask],
            label=str(g) if legend else None,
            color=color_map[g],
            marker=style.get("marker", "o"),
            s=style.get("markersize", 40),
            alpha=style.get("alpha", 0.8),
        )

    # Legend
    if legend:
        ax.legend(title=style.get("legend_title", group_key))

    # Apply common axes styling if available
    if "apply_axes_style" in globals():
        apply_axes_style(ax, style)

    return fig, ax


@register("inter_pairs")
def inter_pairs(
        data,
        style=None,
        fig=None,
        ax=None,
        *,
        x1=None,       # optional 1D array of reviewer 1's values
        x2=None,       # optional 1D array of reviewer 2's values
        y=None,        # optional 1D array of y value
        x1key="Rev1",  # keys for extracting from `data` if arrays not passed directly
        x2key="Rev2",
        ykey="y",
        **_,
):
    """
    Horizontal line plot for each pair of points from same sample
    (inter-reviewer, two measurements/reviews of single slide, both compared to same y-value).
    # to add later - option for markers at edges in addition to the central marker,
                     differentiation between reviewers by color

    Parameters
    ----------
    data : dict-like, optional
        Should contain x1key, x2key and ykey.
    style : dict, optional
        Keys:
          - color: color for central point
          - marker: marker style (default: 'o')
          - markersize: float
          - xlabel, ylabel, title, grid: standard styling keys
    fig, ax : existing figure and axes or None to create new
    x1, x2, y : arrays, optional
        Directly specify coordinates. Overrides keys.
    x1key, x2key, ykey : str
        Keys to extract from `data` if x1/x2/y not given directly.
    """
    style = style or {}
    if ax is None or fig is None:
        fig, ax = plt.subplots()

    # Extract x1, x2, y
    if x1 is None or x2 is None or y is None:
        x1 = np.asarray(data[x1key])
        x2 = np.asarray(data[x2key])
        y = np.asarray(data[ykey])
    else:
        x1 = np.asarray(x1)
        x2 = np.asarray(x2)
        y = np.asarray(y)
        if not (len(x1) == len(x2) == len(y)):
            raise ValueError("x1, x2, y must have the same length.")

    segments = np.column_stack([
        np.stack([x1, y], axis=1),
        np.stack([x2, y], axis=1)
    ]).reshape(-1, 2, 2)

    lc = LineCollection(segments)  # to do: add style kwargs for lines
    ax.add_collection(lc)

    # Apply common axes styling if available
    if "apply_axes_style" in globals():
        apply_axes_style(ax, style)

    return fig, ax

@register("overlay_regression_line")
def overlay_regression_line(
    data=None,
    style: dict | None = None,
    fig=None,
    ax=None,
    *,
    result,                 # RegressionResult: has slope, slope_ci, intercept, intercept_ci
    x=None,                 # optional 1D array; if None uses current x-limits
    n_points: int = 200,    # resolution for drawing the line(s)
    **_,
):
    """
    Overlay a regression line (and optional CI) onto an existing Axes.

    Parameters
    ----------
    data : unused (kept for uniform registry signature)
    style : dict, optional
        Keys (all optional):
          - line_color       : color for regression line (default: palette['fit'] or 'C2')
          - line_width       : float (default: 2.0)
          - line_style       : e.g., '-', '--' (default: '-')
          - label            : legend label for the line
          - ci               : bool, draw confidence interval if True (default: False)
          - ci_mode          : 'lines' or 'shade' (default: 'lines')
          - ci_color         : color (default: line_color)
          - ci_alpha         : alpha for shade (default: 0.15)
          - ci_linestyle     : linestyle for CI lines (default: '--')
          - ci_linewidth     : linewidth for CI lines (default: 1.2)
          - zorder           : zorder for the line (default: None -> matplotlib default)
          - xkey / ykey      : ignored here; present for API symmetry
    fig, ax : existing Matplotlib Figure and Axes (REQUIRED)
    result : object with attributes:
        - slope (float), intercept (float)
        - slope_ci (tuple[low, high]), intercept_ci (tuple[low, high])  # optional for CI
    x : 1D array-like
        Domain to draw across. If None, uses ax.get_xlim().
    n_points : int
        Number of points for the line/CI sampling.
    Returns
    -------
    fig, ax

    Notes
    -----
    CI rendering here combines (slope_ci, intercept_ci) as bounding lines:
        y_low  = slope_ci[0] * x + intercept_ci[0]
        y_high = slope_ci[1] * x + intercept_ci[1]
    This is a visual approximation; a true pointwise CI requires variance/covariance.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    if fig is None or ax is None:
        raise ValueError("overlay_regression_line requires an existing fig and ax.")

    style = style or {}
    # colors/palette fallback
    palette = (DEFAULT_STYLE.get("palette", {}) if "DEFAULT_STYLE" in globals() else {})
    line_color = style.get("line_color", palette.get("fit", "C2"))

    # Domain
    if x is None:
        x_min, x_max = ax.get_xlim()
        xline = np.linspace(x_min, x_max, n_points)
    else:
        x_arr = np.asarray(x)
        x_min, x_max = np.nanmin(x_arr), np.nanmax(x_arr)
        xline = np.linspace(x_min, x_max, n_points)

    # Main regression line
    yline = result.slope * xline + result.intercept
    ax.plot(
        xline,
        yline,
        color=line_color,
        linewidth=style.get("line_width", 2.0),
        linestyle=style.get("line_style", "-"),
        label=style.get("label", "Regression"),
        zorder=style.get("zorder", None),
    )

    # Confidence interval (optional)
    if style.get("ci", False):
        if not (hasattr(result, "slope_ci") and hasattr(result, "intercept_ci") and
                result.slope_ci is not None and result.intercept_ci is not None):
            # If CI requested but missing, don't crash—just skip.
            # You could also raise if you'd prefer strict behavior.
            pass
        else:
            m_lo, m_hi = result.slope_ci
            b_lo, b_hi = result.intercept_ci
            y_low = m_lo * xline + b_lo
            y_high = m_hi * xline + b_hi

            ci_mode = style.get("ci_mode", "lines")
            ci_color = style.get("ci_color", line_color)

            if ci_mode == "shade":
                ax.fill_between(
                    xline,
                    y_low,
                    y_high,
                    alpha=style.get("ci_alpha", 0.15),
                    color=ci_color,
                    linewidth=0,
                    zorder=style.get("zorder", None),
                    label=style.get("ci_label", None),
                )
            else:
                # default 'lines'
                ax.plot(
                    xline, y_low,
                    linestyle=style.get("ci_linestyle", "--"),
                    linewidth=style.get("ci_linewidth", 1.2),
                    color=ci_color,
                    zorder=style.get("zorder", None),
                    label=style.get("ci_label_low", None),
                )
                ax.plot(
                    xline, y_high,
                    linestyle=style.get("ci_linestyle", "--"),
                    linewidth=style.get("ci_linewidth", 1.2),
                    color=ci_color,
                    zorder=style.get("zorder", None),
                    label=style.get("ci_label_high", None),
                )

    # Optional: apply common axes styling if you use a shared helper
    if "apply_axes_style" in globals():
        apply_axes_style(ax, style)

    return fig, ax




