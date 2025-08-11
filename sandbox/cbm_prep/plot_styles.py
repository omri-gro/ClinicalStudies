DEFAULT_STYLE = {
    "grid": True,
    "rc": {               # global rcParams
        "font.size": 11,
    },
    "palette": {          # named colors
        "line": "C0",
        "scatter": "C1",
        "fit": "C2",
    },
    "equal_limits": False,   # if True -> square limits via helper
    "tight_layout": True,
    "legend_loc": "best",
    "legend_title": None,
    "legend_fontsize": None
}

# Per-plot style suggestions (override by runtime kwargs)
PLOT_STYLES = {
    "scatter_basic": {
        "title": "Scatter",
        "xlabel": "X",
        "ylabel": "Y",
        "scatter_color": None  # falls back to palette.scatter
    },
    "scatter_with_fit": {
        "title": "Scatter + Regression",
        "xlabel": "X",
        "ylabel": "Y",
        "line_width": 2.0
    },
    "side_by_side_demo": {
        "figsize": (10, 4),
        "left_title": "Left",
        "right_title": "Right"
    },
    "scatter_groups": {
        "title": "Scatter",
        "xlabel": "X",
        "ylabel": "Y",
        "scatter_color": None,
        "legend": True,
        "legend_title": "Site",
        "legend_loc": "upper right"
    }
}
